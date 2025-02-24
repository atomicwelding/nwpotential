program fourier_coefs
  use iso_fortran_env, only : stdout => output_unit
  use math, only : solve
  
  implicit none
  integer, parameter :: N = 3
  double precision, dimension(N,N) :: A
  double precision, dimension(N) :: BD, BALPHA, BREQ

  BD = [4.843306, 5.603632, 7.613337]
  BALPHA = [1.991931, 0.950659, 1.252572]
  BREQ = [1.736781, 0.785308, 0.565774]

  call reset(N, A)
  call solve(N, A, BD)

  call reset(N, A)
  call solve(N, A, BALPHA)

  call reset(N, A)
  call solve(N, A, BREQ)

  print *, "D"
  print *, BD

  print *, "ALPHA"
  print *, BALPHA

  print *, "REQ"
  print *, BREQ
  
contains
  subroutine reset(N, A)
    integer, intent(in) :: N
    double precision, dimension(N,N), intent(inout) :: A
     A = reshape( [1., 1., 1., &
                2., 0., -2., &
                2., -2., 2.], [N, N] )
  end subroutine reset
end program fourier_coefs
