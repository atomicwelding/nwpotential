! Made by Erwan Le Doeuff
! X, Y are provided in units of delta, so there should belong to [0,1]
program LEPSPES
  implicit none
  double precision, parameter :: pi = 3.141592
  double precision, dimension(3) :: FD, Falpha, Freq

  integer :: Npts
  double precision :: stepXY, stepZ
  double precision, allocatable, dimension(:, :, :) :: PES
  
  integer :: i, j, k
  double precision :: X, Y, Z
  double precision :: Xmin, Ymin, Zmin, Vmin

  Npts = 0
  do while ( Npts < 2 )
     print *, "Npts ( > 1)"
     read(*,*) Npts
  end do


  ! setting things up
  stepXY = 1./Npts
  stepZ = 10./Npts
  allocate(PES(Npts, Npts, Npts))

  ! Fourier coefficients
  FD = [5.9159767628, -0.6925077438, 0.1561723948] 
  Falpha = [1.2864552140, 0.1848397553, 0.1678981185] 
  Freq = [0.9682927579, 0.2927517444, 0.0914923772]
  
  ! computing the map
  PES = 0.
  Xmin = 0.
  Ymin = 0.
  Zmin = 0.
  Vmin = 0.
  loopX: do i=1,Npts
     X = (i-1)*stepXY
     loopY: do j=1,Npts
        Y = (j-1)*stepXY
        loopZ: do k=1,Npts
           Z = (k-1)*stepZ
           PES(i,j,k) = Morse(X, Y, Z)
           if(Vmin > PES(i,j,k)) then
              Vmin = PES(i,j,k)
              Xmin = X
              Ymin = Y
              Zmin = Z
           end if
        end do loopZ
     end do loopY
  end do loopX

  print *, "minimum found at"
  print *, "X = ", Xmin, "(units of delta)"
  print *, "Y = ", Ymin, "(units of delta)"
  print *, "Z = ", Zmin, "(angstrom)"
  print *, "Vmin = ", Vmin, "(eV)"


  ! storing data
  open(unit=10, file="3dpes.txt", status="unknown")
  do i=1, Npts
     X = (i-1) * stepXY
     do j=1, Npts
        Y = (j-1) * stepXY
        do k=1, Npts
           Z = (k-1) * stepZ
           write(10, *) X, Y, Z, PES(i,j,k)
        end do
     end do
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
end program LEPSPES
