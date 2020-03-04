!--------------------------------------------------------------------------------------------------
! initialize quantities needed for sampling
!-------------------------------------------------------------------------------------------------- 
subroutine read_res_file(nhs, nions, ndipoles, npart, lambda, boxlength, density)
use basic_parameters_mod, only: alpha
use working_prec_mod
use sim_parameters_mod
integer, intent(in)     :: nhs
integer, intent(in)     :: nions
integer, intent(in)     :: ndipoles
integer, intent(in)     :: npart
real(rk), intent(in)    :: lambda 
real(rk), intent(in)    :: boxlength
real(rk), intent(in)    :: density
character (len=40)      :: dum_text
integer                 :: nhs_check, nions_check, ndipoles_check, npart_check
real(rk)                :: lambda_check, alpha_check, boxlength_check, density_check, fraction_check
integer                 :: i 
 
    OPEN (res_file_num,FILE = './output_file/'//ResultsFile) 
        
        do i = 1, 18 !(Z23-4 blanks)
            READ (res_file_num,*) dum_text
        end do 
        
        READ (res_file_num,*) max_disp(1), max_disp(2), max_disp(3) 
        READ (res_file_num,*) max_rot(1), max_rot(2), max_rot(3)
        
        do i = 1, 4
            READ (res_file_num,*) dum_text
        end do
        
        READ (res_file_num,*) nhs_check
        READ (res_file_num,*) nions_check
        READ (res_file_num,*) ndipoles_check
        READ (res_file_num,*) npart_check
        
        READ (res_file_num,*) lambda_check
        READ (res_file_num,*) alpha_check
        READ (res_file_num,*) !Temp! 
        READ (res_file_num,*) density_check
        READ (res_file_num,*) boxlength_check 
        READ (res_file_num,*) fraction_check   
   
    CLOSE (res_file_num)   

   if (nhs/=nhs_check .or. nions/=nions_check .or. ndipoles/=ndipoles_check & 
      .or. fraction_check/=real(nions,rk)/real(npart,rk)) then
        print*, 'Check number of particles and fraction'
        print*, nhs_check, nions_check, ndipoles_check, npart_check
        print*, fraction_check,  real(nions)/real(npart)
        stop
    else if (lambda/=lambda_check .or. alpha/=alpha_check .or. density/=density_check ) then 
        print*, 'Check values from res-file and input-file:'
        print*,  lambda_check, lambda        , 'lambda '
        print*,  alpha_check, alpha          , 'alpha '
        print*,  density_check, density      , 'density '
        print*,  boxlength_check, boxlength  , 'boxlength '
        stop
    else 
        print*, ResultsFile,' -file checked. all good. '
    end if  
        
end subroutine read_res_file 

