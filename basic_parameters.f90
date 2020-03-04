module basic_parameters_mod
use working_prec_mod , only: rk, ik 
  implicit none

  
  real(rk), parameter       :: PI = 3.14159265359_rk
  real(rk), parameter       :: ec = 1.60217733d-19  ! elementary charge in C (SI-units)
  real(rk), parameter       :: eps0 = 8.854187817d-12   ! Permittivität des Vakuums  in SI: As/(Vm)  & As=C
  real(rk), parameter       :: kB = 1.380658d-23    ! in J/K
  real(rk)                  :: CoulCombo    
!      CoulCombo = ec * ec * 1.0e10_rk / ( 4.0_rk * PI * eps0 * kB )  ->Umrechnung auf CGS-units (in K)
  
  integer(ik), parameter    :: intervalls = 256_ik  ! number of intervalls for mean radial distribution function
 
  integer(ik), parameter    :: Kmax = 9_ik    ! Kmax is an Ewald sum parameter.
  real(rk), parameter       :: Alpha = 5.6_rk   ! Alpha is an Ewald sum parameter, Alpha = kappa * L, for kappa in A + T. 


end module basic_parameters_mod 
