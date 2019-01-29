# list of executable files to produce
SHIF = shif
SAMP = samp
DIST = dist

# .o files necessary to build the executables
OBJ_SHIF = SimulTree.o  TreeExtras.o NLOpt.o Uncertainty.o  Model.o Utils.o Tree.o TestShiftDetection.o
OBJ_SAMP = Model.o TreeExtras.o Uncertainty.o Utils.o Tree.o Distribution.o Sampling.o Sample.o
OBJ_DIST = DrawDensity.o DrawTreeCairo.o DrawTreePSTricks.o DrawTreeTikz.o DrawTreeGeneric.o Model.o TreeExtras.o TreeBounds.o Uncertainty.o Utils.o Tree.o Distribution.o ComputeDistribution.o
########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT -g

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(SHIF) $(SAMP) $(DIST)

# build the executable
$(SHIF): $(OBJ_SHIF)
	$(LD) $^ -o $@ -lm -lgsl -lgslcblas -lnlopt -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  


# build the executable
$(SAMP): $(OBJ_SAMP)
	$(LD) $^ -o $@ -lm -lgsl -lgslcblas

# build the executable
$(DIST): $(OBJ_DIST)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lnlopt -lcairo

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release
