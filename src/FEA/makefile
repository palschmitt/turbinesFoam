CC=g++
CFLAGS=-I/home/pal/OpenFOAM/OpenFOAM-7/src/OpenFOAM/lnInclude/ \
		-I ../lnInclude

all: main.O FrameAnalysis.o
	$(CC) $(CFLAGS) -o FEATest main.o FrameAnalysis.o

main.O: main.C
	$(CC) $(CFLAGS) -c main.C
	
clean:
	rm *.o FEATest	


