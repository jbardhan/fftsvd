CC = /usr/vac/bin/xlc_r
CFLAGS = -O3 -qstrict -qarch=pwr4 -qtune=pwr4 -q64 -qsmp=omp -qlanglvl=extc99
FLAGS = -DMVTIME -DOMP
LDFLAGS = $(CFLAGS)
INCLUDE = -I/home/ai0/maltman/SuperLU_3.0/SRC -I/home/ai0/maltman/fftw-3.0.1/api
LIBS = /home/ai0/maltman/fftw-3.0.1/.libs/libfftw3.a /home/ai0/maltman/SuperLU_3.0/superlu.a  /home/ai0/maltman/SuperLU_3.0/blas.a /home/ai0/maltman/LAPACK/lapack.a /home/ai0/maltman/LAPACK/blas.a -lm -lf
