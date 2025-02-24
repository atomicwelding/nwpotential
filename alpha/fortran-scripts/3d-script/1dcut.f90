! Made by Erwan Le Doeuff
! X, Y are provided in units of delta, so there should belong to [0,1]
program onedcut
  implicit none
  double precision, parameter :: pi = 3.141592, Delta = 3.
  double precision, dimension(3) :: FD, Falpha, Freq

  integer :: Npts 
  double precision :: stepZ
  double precision, dimension(:), allocatable :: PES
  
  integer :: i
  double precision :: X, Y, Z, zmin, zmax

  Npts = 0
  do while ( Npts < 2 )
     print *, "Npts ( > 1)"
     read(*,*) Npts
  end do

  X = -1.
  do while ( X < 0. .OR. X > 1. )
     print *, "X (1 >= X >= 0)"
     read *, X
  end do

  Y = -1.
  do while ( Y < 0. .OR. Y > 1. )
     print *, "Y (1 >= Y >= 0)"
     read *, Y
  end do

  Zmin = -1.
  do while ( zmin < 0 )
     print *, "Zmin (>= 0, in angstroms)"
     read *, Zmin
  end do

  Zmax = -1.
  do while ( Zmax < Zmin )
     print *, "Zmax (>= Zmin, in angstroms)"
     read *, Zmax
  end do

  stepZ = (zmax-zmin)/(Npts - 1)

  ! Fourier coefficients
  FD = [5.9159767628, -0.6925077438, 0.1561723948] 
  Falpha = [1.2864552140, 0.1848397553, 0.1678981185] 
  Freq = [0.9682927579, 0.2927517444, 0.0914923772] 

  ! computing the map
  allocate(PES(Npts))
  PES = 0.
  do i=1, Npts
     Z = zmin + (i-1) * stepZ
     PES(i) = Morse(X,Y,Z)
  end do

  ! storing data
  open(unit=10, file="1dcut.txt", status="unknown")
  write(10, *) "X = ", X, "Y = ", Y, "Npts =", Npts
  do i=1, Npts
     Z = zmin + (i-1) * stepZ
     write(10, *) Z, PES(i)
  end do
  close(10)

  deallocate(PES)

contains
  function Fourier(Fp, X, Y)
    double precision, dimension(3), intent(in) :: Fp
    double precision, intent(in) :: X, Y

    double precision :: Fourier

    Fourier = 0.
    Fourier = Fourier + Fp(1)
    Fourier = Fourier + Fp(2) * ( cos(2 * pi * X) + cos(2 * pi * Y))
    Fourier = Fourier + Fp(3) * ( cos(2 * pi * (X + Y)) &
         + cos( 2 * pi * (X - Y)) )
  end function Fourier

  double precision function D(X, Y)
    double precision, intent(in) :: X, Y
    D = Fourier(FD, X, Y)
  end function D

  double precision function alpha(X,Y)
    double precision, intent(in) :: X, Y
    alpha = Fourier(Falpha, X, Y)
  end function alpha

  double precision function req(X,Y)
    double precision, intent(in) :: X, Y
    req = Fourier(Freq, X,Y)
  end function req

  double precision function Morse(X, Y, Z)
    double precision, intent(in) :: X, Y, Z
    Morse = exp(- 2 * alpha(X,Y) * (Z - req(X,Y)))
    Morse = Morse - 2 * exp( - alpha(X,Y) * (Z - req(X,Y)))
    Morse = D(X,Y) * Morse
  end function Morse
end program onedcut
