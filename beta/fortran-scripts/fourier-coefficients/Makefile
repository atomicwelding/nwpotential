SRCS = fourier.f90 mod_math.f90
LDFLAGS = -framework Accelerate

.PHONY : fourier


fourier: $(SRCS)
	gfortran -c mod_math.f90
	gfortran -c fourier.f90
	gfortran $(SRCS) -o FOURIER $(LDFLAGS)

clean:
	rm *.o *.mod
	rm FOURIER