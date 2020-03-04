module working_prec_mod
  implicit none
  
  ! working precision
  ! real kind: selected_real_kind(i,j) = i decimals and j exponent-range
  integer, parameter :: rk = selected_real_kind(15,307) !< (15, 307) is double precission 64 bit = 8 byte
  integer, parameter :: rks = selected_real_kind(6,37) !< (6, 37) is single precission 32 bit = 4 byte, 1 bit sign, 8 bits exponent, 23 bits decimals

  ! integer kind: selected_int_kind(i) = -10^i < n < 10^i
  integer, parameter :: ik = selected_int_kind(8) !< kind of integers
  integer, parameter :: ck20 = 20 !< max. number of characters for files

 
end module working_prec_mod
