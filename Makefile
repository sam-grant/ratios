Target  = wiggle_copy.exe IdealWiggle.exe RatioMethod.exe RatioMethodFit.exe Fourier.exe integralTest.exe FFT.exe #wiggle_fit.exe
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
all:$(Target)

wiggle_copy.exe: wiggle_copy.cc
	g++ -o $@ wiggle_copy.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

IdealWiggle.exe: IdealWiggle.cc
	g++ -o $@ IdealWiggle.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

RatioMethod.exe: RatioMethod.cc
	g++ -o $@ RatioMethod.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

RatioMethodFit.exe: RatioMethodFit.cc
	g++ -o $@ RatioMethodFit.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Fourier.exe: Fourier.cc
	g++ -o $@ Fourier.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

integralTest.exe: integralTest.cc
	g++ -o $@ integralTest.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

FFT.exe: FFT.cc
	g++ -o $@ FFT.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
clean:
	rm *exe
