CC = g++
CFLAGS = -Wall -O3 -std=c++17
INC = -I../../include

streamtube_concentration : streamtube_concentration.o
	$(CC) $(CFLAGS) $(LIB) -o streamtube_concentration streamtube_concentration.o
	rm streamtube_concentration.o

streamtube_concentration.o : streamtube_concentration.cpp
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $<

clean :
	rm -f streamtube_concentration.o streamtube_concentration
