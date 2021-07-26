.PHONY: clean
CC = clang++ 
CFLAGS = -g -Wall -std=c++11 -fopenmp -L /usr/local/opt/llvm/lib -I /usr/local/opt/llvm/include

all: channelFlow
	./channelFlow 0.05

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

channelFlow: channelFlow.o
	$(CC) $(CFLAGS) -o $@ $+

clean:
	rm -f *.o core.*
	rm -f channelFlow

