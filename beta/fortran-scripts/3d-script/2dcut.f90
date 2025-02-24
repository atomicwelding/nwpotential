! Made by Erwan Le Doeuff
! X, Y are provided in units of delta, so they should belong to [0,1]
program twodcut
  implicit none
  double precision, parameter :: pi = 3.141592
  double precision, dimension(6) :: FD, Freq
  double precision, dimension(6) :: FAZREP, FAOREP, FAZATT, FAOATT
  double precision, dimension(6) :: FTANHA, FTANHB

  integer :: Npts
  double precision :: stepXY
  double precision, dimension(:,:), allocatable :: PES
  
  integer :: i, j
  double precision :: X, Y, Z

  Npts = 0
  do while ( Npts < 2 )
     print *, "Npts ( > 1)"
     read(*,*) Npts
  end do

  Z = -1.
  do while ( Z < 0 )
     print *, "Z (>= 0, in angstroms)"
     read *, Z
  end do

  stepXY = 1./(Npts - 1)

  ! Fourier coefficients
  FD = [5.9009452164, -0.7354881167, -0.0673756003, 0.2119766772, 0.0491779149, -0.0411859602] 
  Freq = [1.2331064455, 0.3452006206, 0.0191448852, -0.0509045757, -0.0232676901, 0.0246932115] 
  FAZREP = [0.9241847545, 0.2189220041, 0.2550120056, 0.0662515014, 0.0057752654, -0.0026098713] 
  FAOREP = [0.3375306856, 0.1614828743, -0.1379101276, -0.0339673087, -0.0789225623, -0.0354611566] 
  FAZATT = [1.4732959345, -0.2558616251, 0.0309628695, 0.2177459374, 0.0866476968, -0.0376320444] 
  FAOATT = [0.1296798796, 0.1335917544, -0.0031712512, -0.1005041278, -0.0448206263, 0.0108314373] 
  FTANHA = [1.9463933408, 0.5496970117, 0.4992227554, -0.0972003639, 0.1634757370, 0.2211289257] 
  FTANHB = [-1.6245201491, 0.5997645035, 0.3850462511, -0.6448563896, -0.0297759883, 0.1674036849] 

  ! Computing the map
  allocate(PES(Npts, Npts))
  PES = 0.
  loopX: do i=1, Npts
     X = (i-1) * stepXY
     loopY: do j=1, Npts
        Y = (j-1) * stepXY
        PES(i,j) = Morse(X, Y, Z)
     end do loopY
  end do loopX

  ! Storing data
  open(unit=10, file="2dcut.txt", status="unknown")
  write(10, *) "Z = ", Z, "Npts = ", Npts
  do i=1, Npts
     do j=1, Npts
        X = (i-1) * stepXY   
        Y = (j-1) * stepXY 
        write(10, *) X, Y, PES(i,j)
     end do
  end do
  close(10)

  deallocate(PES)

contains
  function Fourier(Fp, X, Y)
    double precision, dimension(6), intent(in) :: Fp
    double precision, intent(in) :: X, Y
    double precision :: Fourier

    Fourier = 0.
    Fourier = Fourier + Fp(1)
    Fourier = Fourier + Fp(2) * ( cos(2 * pi * X) + cos(2 * pi * Y))
    Fourier = Fourier + Fp(3) * ( cos(2 * pi * (X + Y)) &
         + cos(2 * pi * (X - Y)) )
    Fourier = Fourier + Fp(4) * (cos(4 * pi * X) + cos(4 * pi * Y))
    Fourier = Fourier + Fp(5) * (cos(2 * pi * (2 * X + Y)) &
         + cos(2 * pi * (X + 2 * Y)) + cos(2 * pi * (2 * X - Y)) &
         + cos(2 * pi * (X - 2 * Y)))
    Fourier = Fourier + Fp(6) * (cos(4 * pi * (X + Y)) + cos(4 * pi * (X - Y)))
  end function Fourier

  double precision function D(X, Y)
    double precision, intent(in) :: X, Y
    D = Fourier(FD, X, Y)
  end function D
  
  double precision function req(X,Y)
    double precision, intent(in) :: X, Y
    req = Fourier(Freq, X,Y)
  end function req

  double precision function azrep(X,Y)
    double precision, intent(in) :: X, Y
    azrep = Fourier(FAZREP, X,Y)
  end function azrep

  double precision function aorep(X,Y)
    double precision, intent(in) :: X, Y
    aorep = Fourier(FAOREP, X,Y)
  end function aorep

  double precision function azatt(X,Y)
    double precision, intent(in) :: X, Y
    azatt = Fourier(FAZATT, X,Y)
  end function azatt

  double precision function aoatt(X,Y)
    double precision, intent(in) :: X, Y
    aoatt = Fourier(FAOATT, X,Y)
  end function aoatt

  double precision function tanha(X,Y)
    double precision, intent(in) :: X, Y
    tanha = Fourier(FTANHA, X,Y)
  end function tanha

  double precision function tanhb(X,Y)
    double precision, intent(in) :: X, Y
    tanhb = Fourier(FTANHB, X,Y)
  end function tanhb

  double precision function fa(X,Y,Z)
    double precision, intent(in) :: X, Y, Z
    fa = 0.5 * (1 + tanh(tanha(X,Y) * Z + tanhb(X,Y)))
  end function fa

  double precision function alpha_mod(X,Y,Z)
    double precision, intent(in) :: X, Y, Z
    alpha_mod = (1 - fa(X,Y,Z)) * (azrep(X,Y) + aorep(X,Y) * Z) &
         + (fa(X,Y,Z)) * (azatt(X,Y) + aoatt(X,Y) * Z)
  end function alpha_mod

  double precision function Morse(X, Y, Z)
    double precision, intent(in) :: X, Y, Z
    Morse = exp(- 2 * alpha_mod(X,Y,Z) * (Z - req(X,Y)))
    Morse = Morse - 2 * exp( - alpha_mod(X,Y,Z) * (Z - req(X,Y)))
    Morse = D(X,Y) * Morse
  end function Morse

end program twodcut
