! Made by Erwan Le Doeuff
module math
  implicit none

  private
  public :: solve
contains
  subroutine solve(N, A, B)
    integer, intent(in) :: N
    double precision, dimension(N,N), intent(inout) :: A
    double precision, dimension(N), intent(inout) :: B

    integer, dimension(N) :: IPIV
    integer :: info
    
    ! compute the LU factorization
    call dgetrf(N, N, A, N, IPIV, info)
    call handle_error(info, "LU factorization misses")

    ! solve the system Ax=B
    ! stores the result in B
    call dgetrs('N', N, 1, A, N, IPIV, B, N, info)
    call handle_error(info, "Solving linear equations misses")
  end subroutine solve

  subroutine handle_error(info, msg)
    integer, intent(in) :: info
    character(len=*) , intent(in) :: msg

    if (info /= 0) then
       print *, "[ERROR], INFO = ", info
       print *, "=> ", msg
    end if
  end subroutine handle_error
end module math
