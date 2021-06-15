CC = /opt/local/bin/x86_64-apple-darwin11-gcc-4.5.4
CFLAGS = -Wall -g #-O3 -ftree-vectorize 
INCLUDE = -I/opt/local/include
FLAGS = -DMVTIME -DPWD=\"`pwd`\"
LDFLAGS = $(CFLAGS)
LIBDIR =
SUPERLU_INC = 
SUPERLU_LIB = 
LIBS = -L/opt/local/lib -lsuperlu_4.3 -lfftw3 -lfftw3f -llapack -latlas -lcblas -lf77blas -lgfortran -lm
