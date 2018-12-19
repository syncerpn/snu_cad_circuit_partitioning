CC = gcc -pg -g -Wall
C++ = g++ -g -Wall

.c.o: 
	$(CC) -c $<

.cc.o:
	$(C++) -c $<

PARSER = getmodule.o
PARSER += writemodule.o
PARSER += inputparse.o
PARSER += hash.o
PARSER += global.o
PARSER += util.o
PARSER += partition_alg.o

PART = pmain.o

install:
	make part

inputparse.o: 	standard.h inputparse.c
global.o: 		standard.h global.c
pmain.o:		standard.h pmain.c

part: $(PARSER) $(PART) 
	$(CC) -o part $(PARSER) $(PART) -lm
	cp part ../test

clean:
	rm -rf $(PARSER) $(PART) part
