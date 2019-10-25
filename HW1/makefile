CC=g++
CFLAGS=-O2 -ansi -pedantic
## for Mac OS X
LIBS=-framework Accelerate
## For Linux
#LIBS=-llapack -lblas -lg2c -lm
LAPACKFLAGS=`./lapackflags.bash`

packet: packet.cpp
	$(CC) -o $@ $@.cpp $(CFLAGS) $(LAPACKFLAGS) $(LIBS) 

clean:
	rm packet *.dat 2>&-
