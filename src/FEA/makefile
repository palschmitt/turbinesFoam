CC=g++

LIB_SRC=/home/pal/OpenFOAM/OpenFOAM-7/src

CFLAGS=-I ../lnInclude \
		-I$(LIB_SRC)/OpenFOAM/lnInclude \
		-I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
		-I$(LIB_SRC)/sampling/lnInclude \
	    -I$(LIB_SRC)/thermophysicalModels/solidThermo/lnInclude \
	    -I$(LIB_SRC)/transportModels/compressible/lnInclude \
	    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
	    -I$(LIB_SRC)/thermophysicalModels/specie/lnInclude \
	    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
	    -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
	    -I$(LIB_SRC)/fvOptions/lnInclude 
	    
LCFLAGS=-larmadillo -llapack -lblas\
    -lfiniteVolume \
    -lsampling \
    -lmeshTools \
    -lturbulenceModels \
    -lcompressibleTurbulenceModels \
    -lfvOptions

all: main.o FrameAnalysis.o
	$(CC) $(CFLAGS) $(LCFLAGS) -o FEATest main.o FrameAnalysis.o
	
FrameAnalysis.o: 	FrameAnalysis.C FrameAnalysis.H
	$(CC) $(CFLAGS)  -c FrameAnalysis.C 

main.O: main.C
	$(CC) $(CFLAGS)  -c main.C
	
clean:
	rm *.o FEATest	


