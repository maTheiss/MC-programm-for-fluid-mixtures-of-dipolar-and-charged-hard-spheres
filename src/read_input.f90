module read_input_mod
use working_prec_mod, only: rk, ik 
use sim_parameters_mod   
  implicit none
  private

  public :: read_input    
  public :: identifier_files
  public :: restarted_files


contains
  
!--------------------------------------------------------------------------------------------------
! input data
!--------------------------------------------------------------------------------------------------
subroutine read_input(read_conf, first_run, nprod_in, nequil_in, nsamp_in, nhs_in, nions_in &
            ,ndipoles_in, Temp_in, density_in, q_ions_in, mu_dipoles_in, &
            sigma_in, lambda_in )

logical     :: read_conf, first_run 
integer(ik) :: nequil_in, nprod_in, nsamp_in 
integer(ik) :: nhs_in, nions_in, ndipoles_in
real(rk)    :: Temp_in, density_in
real(rk)    :: q_ions_in, mu_dipoles_in
real(rk)    :: sigma_in, lambda_in 
character (len=40) ::  dum_text
! character (len=40) :: in_file 

    call get_command_argument(1,in_file) 
     
    OPEN (30,FILE = in_file)   !'../input/input.inp')

    READ (30,*) dum_text 
    
    READ (30,*) dum_text, read_conf
    READ (30,*) dum_text, first_run 
    READ (30,*) dum_text, nprod_in
    READ (30,*) dum_text, nequil_in
    READ (30,*) dum_text, nsamp_in 
    READ (30,*) dum_text, nhs_in
    READ (30,*) dum_text, nions_in
    READ (30,*) dum_text, ndipoles_in
    READ (30,*) dum_text, q_ions_in
    READ (30,*) dum_text, mu_dipoles_in
    READ (30,*) dum_text, sigma_in  
    READ (30,*) dum_text, lambda_in 

    READ (30,*) dum_text
    
    READ (30,*) dum_text, Temp_in
    READ (30,*) dum_text, density_in  
    
    CLOSE (30)
   
    
end subroutine read_input


!--------------------------------------------------------------------------------------------------
! Names for input/output files 
!--------------------------------------------------------------------------------------------------
subroutine identifier_files 
integer(ik)        :: Dotpos

    Dotpos = SCAN(in_file,'.')

    if ( Dotpos .NE. 0 ) then
        SimName = in_file( 1 : DotPos-1 )
    else
        SimName = in_file
    end if


    ResultsFile = trim(SimName)//'.res'
    UavFile = trim(SimName)//'.plt'
    FinalConf_PartPos = trim(SimName)//'.ppos'
    FinalConf_ChargePos = trim(SimName)//'.chpos'
    RunnEnFile = trim(SimName)//'.ren'
    
!     PairCorrFile = trim(SimName)//'.rdf'
!     KirkwoodFile = trim(SimName)//'.kwf'

end subroutine identifier_files

!--------------------------------------------------------------------------------------------------
! Names for input/output files after restarted simulation 
!--------------------------------------------------------------------------------------------------
subroutine restarted_files
 
     
!     if ( read_conf ) then

    ResultsFile = trim(SimName)//'_rst.res'
    UavFile = trim(SimName)//'_rst.plt'
    FinalConf_PartPos = trim(SimName)//'_rst.ppos'
    FinalConf_ChargePos = trim(SimName)//'_rst.chpos' 
    RunnEnFile = trim(SimName)//'_rst.ren'
    
!     PairCorrFile = trim(SimName)//'_rst.rdf'
!     KirkwoodFile = trim(SimName)//'_rst.kwf'

!     end if
end subroutine restarted_files



end module read_input_mod