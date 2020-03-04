!> \file main.f90 
!> \brief MC Code to compute the internal energy (with Ewald summation) for fluid mixtures of dipolar and charged hard spheres 
!> \author M.Theiss 

!> main programm 
program main
  use working_prec_mod, only: rk, ik
  use globals_mod  
  use ran_gen_mod     , only: init_random_seed
  use nvt_mod         , only: nvt
  use read_input_mod  , only: read_input 
  implicit none
 
  logical     :: read_conf, first_run
  integer(ik) :: nequil_in, nprod_in
  integer(ik) :: nhs_in, nions_in, ndipoles_in
  real(rk)    :: Temp_in, density_in, sigma_in  

  !> reading data from input file 
  call read_input (read_conf, first_run, nprod_in, nequil_in, nsamp, nhs_in, nions_in &
            , ndipoles_in, Temp_in, density_in  &
            , q_ions, mu_dipoles, sigma_in, lambda )
 
  !> initialize the random seed
  call init_random_seed ()

  !> initialize system:
  ! calculate volume and boxlength
  ! set particles on lattice
  call init_globals (read_conf, first_run, nhs_in, nions_in, ndipoles_in &
    , density_in, Temp_in, sigma_in )
  
  !> MC 
  call nvt (nequil_in, nprod_in)  
   
end program main
