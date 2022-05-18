MODULENAME=cons1d
F90= gfortran
CC = gcc
CXX= g++
AR = ar -rvs
CFLAGS = -g #-fpe0
FFLAGS =  -fbacktrace -fbounds-check
INCLUDES = 
OBJF90 = shape.o amodules.o apde.o clearMem.o computeMoments.o computeRHS.o init.o \
         initvar.o matrixsolve.o output.o overset.o solvedq.o timestep.o slopelimiter.o  xdriver.o

LDFLAGS=

default: $(OBJECTS) $(OBJF90) $(INCLUDES)
	$(F90) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -r *.o *.mod $(MODULENAME).exe

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o
