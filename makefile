.PHONY: clean
CC = g++ 
CFLAGS = -g -Wall -std=c++11 -fopenmp 

all: channelFlow
	./channelFlow 0.05

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

channelFlow: channelFlow.o
	$(CC) $(CFLAGS) -o $@ $+

clean:
	rm -f *.o core.*
	rm -f channelFlow

