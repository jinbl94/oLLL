#include <iostream>
#include <string>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>

using namespace NTL;

void sort(vec_long& P, const long m, const vec_RR& D, const vec_ZZ& V2);

void updateM(const long i, const mat_ZZ& B, const long m, mat_RR& M, vec_ZZ& V2, vec_RR& D);

void init_param(const mat_ZZ& B, const long m, mat_RR& M, vec_ZZ& V2, vec_RR& D);

void reduce(mat_ZZ& B, vec_long P, const long m, mat_RR& M, vec_ZZ& V2, vec_RR& D, long& counter, bool& verbose);

static long oLLL(mat_ZZ& B, bool verbose = false);

int main(int argc, char* argv[])
{

    if(argc > 1){
        std::cout << "Parameters: ";
        for(int i = 1; i < argc; i++){
            std::cout << argv[i] << ' ';
        }
        std::cout << std::endl;
    }

    int loop;
    const long m = 10;
    ZZ BOUND;
    BOUND = 100;
    mat_ZZ B;
    B.SetDims(m, m);

    for(long i = 0; i < m; i++){
        for(long j = 0; j < m; j++){
            RandomBnd(B[i][j], BOUND);
        }
    }

    loop = oLLL(B, true);

    std::cout << "\nloop: " << loop << std::endl;

    return loop;
}

void sort(vec_long& P, const long m, const vec_RR& D, const vec_ZZ& V2)
// exchange P[i] and P[j] (i > j) if D[P[i]] < D[P[j]]
{
    long p = 0, t;
    for(long i = 1; i < m; i++){
        if(V2[P[i]] < V2[P[p]]){
            p = i;
        }
        RR t_d = D[P[i]];
        long t_p = P[i], j = i - 1;
        for(; j >= 0; j--){
            if(t_d < D[P[j]]){
                P[j + 1] = P[j];
            }else{
                break;
            }
        }
        P[j + 1] = t_p;
    }
    t = P[0];
    P[0] = P[p];
    P[p] = t;
}

void updateM(const long i, const mat_ZZ& B, const long m, mat_RR& M, vec_ZZ& V2, vec_RR& D)
{
    ZZ t_zz;
    RR t_rr;

    // V2[i] = B[i]^2
    InnerProduct(t_zz, B[i], B[i]);
    conv(V2[i], t_zz);
    
    // subtract i related parts
    for(long k = 0; k < m; k++){
        sqr(t_rr, M[k][i]);
        D[k] -= t_rr;
    }
    D[i] = 0;

    long j;
    for(j = 0; j < i; j++){
        InnerProduct(t_zz, B[i], B[j]);
        conv(M[i][j], t_zz);
        M[j][i] = M[i][j];

        conv(t_rr, V2[j]);
        div(M[i][j], M[i][j], t_rr);
        sqr(t_rr, M[i][j]);
        D[i] += t_rr;

        conv(t_rr, V2[i]);
        div(M[j][i], M[j][i], t_rr);
        sqr(t_rr, M[j][i]);
        D[j] += t_rr;
    }
    for(j = j + 1; j < m; j++){
        InnerProduct(t_zz, B[i], B[j]);
        conv(M[i][j], t_zz);
        M[j][i] = M[i][j];

        conv(t_rr, V2[j]);
        div(M[i][j], M[i][j], t_rr);
        sqr(t_rr, M[i][j]);
        D[i] += t_rr;

        conv(t_rr, V2[i]);
        div(M[j][i], M[j][i], t_rr);
        sqr(t_rr, M[j][i]);
        D[j] += t_rr;
    }
}

void init_param(const mat_ZZ& B, const long m, mat_RR& M, vec_ZZ& V2, vec_RR& D)
// refresh all parameters
{
    ZZ t_zz;
    RR t_rr;

    for(long i = 0; i < m; i++){
        InnerProduct(t_zz, B[i], B[i]);
        conv(V2[i], t_zz);
        D[i] = 0;
    }

    for(long i = 1; i < m; i++){
        for(long j = 0; j < i; j++){
            InnerProduct(t_zz, B[i], B[j]);
            conv(M[i][j], t_zz);
            M[j][i] = M[i][j];

            // M[i][j] = <B[i], B[j]> / <B[j], B[j]>
            // D[i] = D[i] + M[i][j]^2
            conv(t_rr, V2[j]);
            div(M[i][j], M[i][j], t_rr);
            sqr(t_rr, M[i][j]);
            D[i] += t_rr;

            // M[j][i] = <B[i], B[j]> / <B[i], B[i]>
            // D[j] = D[j] + M[j][i]^2
            conv(t_rr, V2[i]);
            div(M[j][i], M[j][i], t_rr);
            sqr(t_rr, M[j][i]);
            D[j] += t_rr;
        }
    }
}

void reduce(mat_ZZ& B, const vec_long P, const long m, const long n, 
        mat_RR& M, vec_ZZ& V2, vec_RR& D, long& counter, bool& verbose)
// round jth basis vector and insert it to 'correct' position, return the index of next basis vector
{
    counter = 0;
    ZZ factor;
    RR t, t_z;
    for(long i = 1; i < m; i++){
        for(long j = 0; j < i; j++){
            // calculate roundoff factor according to the 'integer' part of mu_{i,j}
            t = M[P[i]][P[j]];
            round(t_z, t);
            if(verbose){
                std::cout << "mu(" << P[i] << "," << P[j] << "): " << t
                    << " t_z: " << t_z
                    << "\n";
            }
            if(t_z == 0){
                continue;
            }
            conv(factor, t_z);

            // B_{P_i} = B_{P_i} - B_{P_j} * factor
            for(long k = 0; k < n; k++){
                B[P[i]][k] -= (B[P[j]][k] * factor);
            }
            counter++;
        }
        //refresh B[P[i]] related parameters
        updateM(P[i], B, m, M, V2, D);
    }
}

static long oLLL(mat_ZZ& B, bool verbose)
{
    const long m = B.NumRows();
    const long n = B.NumCols();

    // permutation of B
    vec_long P;
    P.SetLength(m);
    for(long i = 0; i < m; i++){
        P[i] = i;
    }

    // {(b_i)^2}_{i \in [m]}
    vec_ZZ V2;
    V2.SetLength(m);

    // {sum_j \mu_{i,j}^2}_{i \in [m]}
    vec_RR D;
    D.SetLength(m);
    vec_RR D_1;
    D.SetLength(m);

    // {<b_i, b_j>}_{i,j \in [m]}
    mat_RR M;
    M.SetDims(m, m);

    init_param(B, m, M, V2, D);
    if(verbose){
        std::cout << "Before reduction:\n" << "B:\n" << B
            << "\nV2:\n" << V2 << " V2[P[0]]: " << V2[P[0]]
            << "\nM:\n" << M
            << "\nD:\n" << D
            << "\n\n";
    }

    long loopcount = 0;
    long counter = 1;
    while(counter){
        sort(P, m, D, V2);
        reduce(B, P, m, n, M, V2, D, counter, verbose);
        loopcount++;
        if(loopcount > 3)
            break;
    }

    if(verbose){
        std::cout << "After reduction:\n" << "B:\n" << B
            << "\nV2:\n" << V2 << " V2[P[0]]: " << V2[P[0]]
            << "\nM:\n" << M
            << "\nD:\n" << D
            << "\n\n";
    }

    return loopcount;
}
