FC = gfortran

CC = g++

CXX = $(CC)

# If you link using gfortran, then use the FORTRAN_LIB_ATLAS. If you link 
# with g++, then use the CPP_LIB_ATLAS. Basically, I'm using a mixture of 
# libraries, so whichever linker I use must know which libraries of the 
# other language to include. 
FORTRAN_LIB_ATLAS = -llapack -lf77blas -lcblas -latlas -lstdc++

CPP_LIB_ATLAS = -llapack -lgfortranbegin -lgfortran -lblas

LDFLAGS = $(CPP_LIB_ATLAS)

SOURCES = ComputationalElement.cpp Fitter.cpp CircularPistonTransducerElement.cpp \
CircularPistonTransducer.cpp OptimizeTransducer.cpp Runner.cpp main.cpp 

OBJ = ComputationalElement.o Fitter.o CircularPistonTransducerElement.o \
CircularPistonTransducer.o OptimizeTransducer.o Runner.o main.o 

CFLAGS = -O2 -Wall
CPPFLAGS = $(CFLAGS) 

EXECUTABLE_NAME = go

PREPROCESSOR = -DDEBUG 

all : main.cpp Runner.h OptimizeTransducer.h CircularPistonTransducer.h \
CircularPistonTransducerElement.h TransducerElement.h ComputationalElement.h \
Fitter.h array3d.h array3d.cpp $(OBJ)
	$(CC) $(OBJ) $(CPP_LIB_ATLAS) $(PREPROCESSOR) -o $(EXECUTABLE_NAME)
	
main.o : main.cpp 
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c main.cpp

runner : runner.cpp runner.h
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c runner.cpp

ComputationalElement.o : ComputationalElement.h ComputationalElement.cpp 
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c ComputationalElement.cpp

Fitter.o : Fitter.cpp Fitter.h
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c Fitter.cpp

CircularPistonTransducerElement.o : CircularPistonTransducerElement.cpp \
CircularPistonTransducerElement.h TransducerElement.h ComputationalElement.h
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c CircularPistonTransducerElement.cpp
 
CircularPistonTransducer.o : CircularPistonTransducer.cpp CircularPistonTransducer.h \
CircularPistonTransducerElement.h TransducerElement.h ComputationalElement.h \
Fitter.h array3d.h array3d.cpp
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c CircularPistonTransducer.cpp
 
OptimizeTransducer.o : OptimizeTransducer.cpp OptimizeTransducer.h TransducerElement.h \
ComputationalElement.h Fitter.h    
	$(CC) $(CFLAGS) $(PREPROCESSOR) -c OptimizeTransducer.cpp
	
clean:
	rm -f $(OBJ) $(EXECUTABLE_NAME) 
