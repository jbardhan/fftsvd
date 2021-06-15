CC = gcc
CFLAGS = -g -D_MALLOC
INCLUDE = -I$(SUPERLU_INC)
FLAGS = #-DMVTIME
LDFLAGS = $(CFLAGS)
LIBDIR =
SUPERLU_INC = /usr/include/superlu
SUPERLU_LIB = /usr/lib
LIBS = -L$(SUPERLU_LIB) -lsuperlu_5.2.2 -llapack -lblas -lfftw3 -lfftw3f -lm


#CFLAGS = -g -O3 -march=opteron
#CFLAGS =  -g -O3 -mtune=pentium4 -mmmx -msse -msse2 -mfpmath=sse,387
