CC = g++
CFLAGS = -Wall -O3 -std=c++17
INC = -I../../include

streamtube_gillespie : streamtube_gillespie.o
	$(CC) $(CFLAGS) $(LIB) -o streamtube_gillespie streamtube_gillespie.o
	rm streamtube_gillespie.o

streamtube_gillespie.o : streamtube_gillespie.cpp
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $<

clean :
	rm -f streamtube_gillespie.o streamtube_gillespie
