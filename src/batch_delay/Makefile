CC = g++
CFLAGS = -Wall -O3 -std=c++17
INC = -I../../include

batch_delay : batch_delay.o
	$(CC) $(CFLAGS) $(LIB) -o batch_delay batch_delay.o
	rm batch_delay.o

batch_delay.o : batch_delay.cpp
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $<

clean :
	rm -f batch_delay.o batch_delay
