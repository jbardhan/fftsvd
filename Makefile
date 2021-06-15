include make-$(HOSTCOMPUTER).inc

OBJS = Vector.o ComplexVector.o Matrix.o Vector3D.o Panel.o Integration.o \
       VertFace.o QUI.o Charge.o GMRES.o Cube.o Tree.o \
       FFT.o Preconditioner.o EquivDensity.o GreensFunction.o \
       QuadratureRule.o LJparameters.o SurfaceOperator.o QualocationOperator.o \
       FlatPanel.o GST.o TOR.o FlatIntegration.o calcpc_GST.o calcpc_TOR.o \
       Polynomial.o SMatrix.o SVector.o ComplexSVector.o FFTSVDpbeAPI.o

all: FFTSVDcap FFTSVDsolv FFTSVDlj FFTSVDstern FFTSVDsolvcav FFTSVDpbe \
     FFTSVDsalt FFTSVDforce FFTSVDsolvecfqual FFTSVDsolvcurved \
     FFTSVDcapcurved FFTSVDsolvecfqualcurved FFTSVDcoul FFTSVDsolvecfqualcav \
     FFTSVDljcurved FFTSVD_GBall FFTSVDcurvedhessian curvedarea meshmaker testdlp triangulatecurved

clean:
	-rm FFTSVDcap FFTSVDsolv FFTSVDlj FFTSVDstern FFTSVDsolvcav FFTSVDpbe \
            FFTSVDsalt FFTSVDforce FFTSVDsolvecfqual meshmaker \
            FFTSVDsolvcurved FFTSVDcapcurved FFTSVDsolvecfqualcurved \
            FFTSVDcoul FFTSVDsolvecfqualcav FFTSVDljcurved FFTSVD_SGB \
            FFTSVD_GBall FFTSVDcurvedhessian curvedarea *.o *.il testm3

fresh: clean all

obj: $(OBJS)

Matrixdriver: Matrixdriver.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ Matrixdriver.o $(OBJS) $(LIBDIR) $(LIBS)

dqdriver: dqdriver.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ dqdriver.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVD_SGB: FFTSVD_SGB.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVD_SGB.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVD_GBall: FFTSVD_GBall.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVD_GBall.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDpot: FFTSVDpot.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDpot.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDcurvedhessian: FFTSVDcurvedhessian.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDcurvedhessian.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDhessian: FFTSVDhessian.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDhessian.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDforce: FFTSVDforce.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDforce.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDcap: FFTSVDcap.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDcap.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDcapcurved: FFTSVDcapcurved.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDcapcurved.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolv: FFTSVDsolv.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolv.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDbcadiag: FFTSVDbcadiag.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDbcadiag.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvcurved: FFTSVDsolvcurved.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvcurved.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvecf: FFTSVDsolvecf.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvecf.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvecfqual: FFTSVDsolvecfqual.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvecfqual.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvecfqualcav: FFTSVDsolvecfqualcav.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvecfqualcav.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvecfqualcurved: FFTSVDsolvecfqualcurved.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvecfqualcurved.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDlj: FFTSVDlj.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDlj.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDljcurved: FFTSVDljcurved.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDljcurved.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDstern: FFTSVDstern.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDstern.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsalt: FFTSVDsalt.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsalt.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDsolvcav: FFTSVDsolvcav.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDsolvcav.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDpbe: FFTSVDpbe.o FFTSVDpbeAPI.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDpbe.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDcoul: FFTSVDcoul.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ FFTSVDcoul.o $(OBJS) $(LIBDIR) $(LIBS)

selfint: selfint.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ selfint.o $(OBJS) $(LIBDIR) $(LIBS)

curvedarea: curvedarea.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ curvedarea.o $(OBJS) $(LIBDIR) $(LIBS)

meshmaker: meshmaker.o FFTSVDpbeAPI.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ meshmaker.o $(OBJS) $(LIBDIR) $(LIBS)

meshmaker2: meshmaker2.o FFTSVDpbeAPI.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ meshmaker2.o $(OBJS) $(LIBDIR) $(LIBS)

integrationtiming: integrationtiming.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ integrationtiming.o $(OBJS) $(LIBDIR) $(LIBS)

accuracy: accuracy.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ accuracy.o $(OBJS) $(LIBDIR) $(LIBS)

testdlp: testdlp.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ testdlp.o $(OBJS) $(LIBDIR) $(LIBS)

testM3: testM3.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ testM3.o $(OBJS) $(LIBDIR) $(LIBS)

checkgst: checkgst.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ checkgst.o $(OBJS) $(LIBDIR) $(LIBS)

triangulatecurved: triangulatecurved.o $(OBJS);
	$(CC) $(LDFLAGS) -o $@ triangulatecurved.o $(OBJS) $(LIBDIR) $(LIBS)

FFTSVDpbeAPI.o: FFTSVDpbeAPI.c
	$(CC) $(CFLAGS) $(INCLUDE) -DPWD=\"`pwd`\" $(FLAGS) -c FFTSVDpbeAPI.c

meshmaker.o:
	$(CC) $(CFLAGS) $(INCLUDE) -DPWD=\"`pwd`\" $(FLAGS) -c meshmaker.c

meshmaker2.o:
	$(CC) $(CFLAGS) $(INCLUDE) -DPWD=\"`pwd`\" $(FLAGS) -c meshmaker2.c

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c $<

