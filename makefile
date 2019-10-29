TARGET=olll
OBJS=olll.o
CC=g++
FLAGS= -lntl -g

default:$(TARGET)
$(TARGET):$(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(FLAGS)
$(OBJS):%.o:%.cpp
	$(CC) -c $< -o $@

.PHONY:clean
clean:
	@rm -f $(OBJS) $(TARGET)
