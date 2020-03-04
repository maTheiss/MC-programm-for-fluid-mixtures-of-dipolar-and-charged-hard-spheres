module util_mod
  use working_prec_mod, only: ik, rk
  implicit none
  private

  public :: det
  public :: inverse
  public :: trace
  public :: math_cross_product
  public :: math_euclidean_norm
contains

  ! det: returns determinant of 3x3 matrix (sarrus)
  pure real(rk) function det (matrix)
    real(rk), dimension(3,3), intent(in) :: matrix
    det  = matrix(1,1)*matrix(2,2)*matrix(3,3) &
      + matrix(1,2)*matrix(2,3)*matrix(3,1) &
      + matrix(1,3)*matrix(2,1)*matrix(3,2) &
      - matrix(1,3)*matrix(2,2)*matrix(3,1) &
      - matrix(1,1)*matrix(2,3)*matrix(3,2) &
      - matrix(1,2)*matrix(2,1)*matrix(3,3)
  end function det

  pure function inverse(matrix) result (invbox)
    real(rk), intent(in) :: matrix(3,3)
    real(rk)             :: invbox(3,3)
    real(rk)             :: adj(3,3) ! adjunkt

    adj(1,1) = matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)
    adj(2,1) = matrix(3,1)*matrix(2,3) - matrix(3,3)*matrix(2,1)
    adj(3,1) = matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)
    adj(1,2) = matrix(3,2)*matrix(1,3) - matrix(3,3)*matrix(1,2)
    adj(2,2) = matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1)
    adj(3,2) = matrix(3,1)*matrix(1,2) - matrix(3,2)*matrix(1,1)
    adj(1,3) = matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2)
    adj(2,3) = matrix(2,1)*matrix(1,3) - matrix(2,3)*matrix(1,1)
    adj(3,3) = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)

    invbox = 1._rk/det(matrix) * adj
  end function inverse

  pure real(rk) function trace (matrix)
    real(rk), intent(in) :: matrix(3,3)

    trace = ( matrix(1,1) + matrix(2,2) + matrix(3,3) )
  end function trace

  pure function math_cross_product (a, b) result (res)
    real(rk), dimension(3), intent(in) :: a, b
    real(rk), dimension(3) :: res

    res(1) = a(2) * b(3) - a(3) * b(2)
    res(2) = a(3) * b(1) - a(1) * b(3)
    res(3) = a(1) * b(2) - a(2) * b(1)
  end function math_cross_product

  pure function math_euclidean_norm (a) result (a_n)
    real(rk), dimension(:), intent(in) :: a
    real(rk)                           :: a_n

    a_n = sqrt (sum (a**2))
  end function math_euclidean_norm

  
end module util_mod

