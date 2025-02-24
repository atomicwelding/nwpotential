.PHONY : all 3dpes 2dcut 1dcut clean

all: 3dpes 2dcut 1dcut

3dpes: 3dpes.f90
	gfortran 3dpes.f90 -o 3DPES

2dcut: 2dcut.f90
	gfortran 2dcut.f90 -o 2DCUT

1dcut: 1dcut.f90
	gfortran 1dcut.f90 -o 1DCUT

clean:
	rm *.o *.mod
	rm 3DPES 1DCUT 2DCUT
