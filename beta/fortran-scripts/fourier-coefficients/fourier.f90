program fourier_coefs
  use iso_fortran_env, only : stdout => output_unit
  use math, only : solve
  
  implicit none
  integer, parameter :: N = 6
  double precision, dimension(N,N) :: A

  double precision, dimension(N) :: BD, BREQ
  double precision, dimension(N) :: BAZATT, BAZREP, BAOATT, BAOREP ! BAZ->A0 , BAO->A1
  double precision, dimension(N) :: BTANHA, BTANHB

  ! Fitted values
  BD = [4.833511, 6.377278, 7.382040, 5.149473, 5.394620, 6.817161]
  BAZREP = [2.022437, 0.541444, 1.100547, 1.136776, 0.786462, 0.722033]
  BAOREP = [-0.069871, 0.474494, -0.084422, 0.727781, 0.334543, 0.089125]
  BAZATT = [1.730317, 1.771598, 2.060582, 1.119403, 0.962540, 1.977717]
  BAOATT = [0.031893, -0.043323, -0.143909, 0.331250, 0.352351, -0.115216]
  BREQ = [1.816304, 1.142394, 0.621643, 1.575456, 1.384302, 0.791984]
  BTANHA = [4.945993, 1.195805, 1.439399, 1.726881, 2.583052, 1.281390]
  BTANHB = [-0.728908, -3.349518, -2.889758, -1.300011, 0.000000, -2.618644]


  ! Solving
  call init(N, A)
  call solve(N, A, BD)
  print *, "D"
  call disp(N, BD)

  call init(N, A)
  call solve(N, A, BAZREP)
  print *, "A0 REP"
  call disp(N, BAZREP)

  call init(N, A)
  call solve(N, A, BAOREP)
  print *, "A1 REP"
  call disp(N, BAOREP)

  call init(N, A)
  call solve(N, A, BAZATT)
  print *, "A0 ATT"
  call disp(N, BAZATT)

  call init(N, A)
  call solve(N, A, BAOATT)
  print *, "A1 ATT"
  call disp(N, BAOATT)

  call init(N, A)
  call solve(N, A, BREQ)
  print *, "REQ"
  call disp(N, BREQ)

  call init(N, A)
  call solve(N, A, BTANHA)
  print *, "TANH A"
  call disp(N, BTANHA)

  call init(N, A)
  call solve(N, A, BTANHB)
  print *, "TANH B"
  call disp(N, BTANHB)
  
contains
  subroutine init(N, A)
    integer, intent(in) :: N
    double precision, dimension(N,N), intent(inout) :: A
    A = reshape([1., 1., 1., 1., 1., 1., &
         2., 0., -2., 1., 0., -1., &
         2., -2., 2., 0., 0., 0., &
         2., 2., 2., 0., -2., 0., &
         4., 0., -4., -2., 0., 2., &
         2., 2., 2., -2., 2., -2.], &
         [N,N])
  end subroutine init

  subroutine disp(N, x)
    integer, intent(in) :: N
    double precision, dimension(N), intent(in) :: x

    integer :: i
    do i=1,N
       print *, x(i)
    end do
  end subroutine disp
end program fourier_coefs
