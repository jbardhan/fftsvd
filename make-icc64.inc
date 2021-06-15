CC = icc
CFLAGS = -O3 -xW -unroll -ip -mp -openmp
INCLUDE = -I/usr/local/include -I/usr/local/include/SuperLU -I/programs/i386/lib/intel/ipp/5.0/ia32/include
FLAGS = -DINTEL_C -DMVTIME -DOMP
LDFLAGS = $(CFLAGS) -static-libcxa -i-static
LIBDIR = -L/usr/local/lib64 -L/programs/i386/lib/intel/ipp/5.0/em64t/lib -L/programs/i386/lib/intel/mkl/8.0.1/lib/em64t
LIBS = -lsuperlu /usr/lib64/liblapack.a /usr/lib64/libblas.a -lg2c -lfftw3 -lfftw3f -lippsemergedem64t -lippsmergedem64t -lippcoreem64t
