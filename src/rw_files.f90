module rw_files_mod
  use working_prec_mod, only: rk, ik  
  use sim_parameters_mod 
  use globals_mod
  implicit none
  private
  
  public :: sample
  public :: store_part_config
  public :: store_ch_config
  public :: store_results
  public :: initialize_sampling  
  public :: store_running_energy  
 
contains
 
 
!--------------------------------------------------------------------------------------------------
! this routine creates a (plot)file including total attempted moves and average energy 
!--------------------------------------------------------------------------------------------------
subroutine sample(i, en1, en2 ) !, vir, press)  
integer(ik), intent(inout)                   :: i
real(rk), intent(in)                         :: en1  
real(rk), intent(in)                         :: en2  
real(rk)                                     :: enpt1, enpt2
 
!     IF ( i == nsamp ) WRITE (66, *) real(npart,rk), Temp
    IF (npart /= 0) THEN
        enpt1   = en1/real(npart,rk) 
        enpt2   = en2/real(npart,rk) 
    ELSE
        enpt1   = 0.0D0 
        enpt2   = 0.0D0 
    END IF

    WRITE (plt_file_num, *) i, enpt1, enpt2  

  RETURN
end subroutine sample


!--------------------------------------------------------------------------------------------------
! this routine creates a (plot)file for running energy 
!--------------------------------------------------------------------------------------------------
subroutine store_running_energy(i, en) !, vir, press)  
integer(ik), intent(inout)                   :: i
real(rk), intent(in)                         :: en  
real(rk)                                     :: enpt
 
!     IF ( i == nsamp ) WRITE (66, *) real(npart,rk), Temp
    IF (npart /= 0) THEN
        enpt   = en/real(npart,rk)
    ELSE
        enpt   = 0.0D0 
    END IF

    WRITE (re_file_num, *) i, enpt 

  RETURN
end subroutine store_running_energy


!--------------------------------------------------------------------------------------------------
! writing U-average with standard deviation and input data
!--------------------------------------------------------------------------------------------------
subroutine store_results(nequil, nprod, nsamp_blocks, Usum, total_att, Utot) 
 integer(ik), intent(inout)                          :: nequil 
 integer(ik), intent(inout)                          :: nprod
 integer(ik), intent(in)                             :: nsamp_blocks
 real(rk), dimension(nsamp_blocks), intent(inout)    :: Usum
 integer(ik), dimension(nsamp_blocks), intent(inout) :: total_att
 real(rk), intent(inout)                             :: Utot
 real(rk)   :: U_average  = 0.0_rk
 real(rk)   :: stdev      = 0.0_rk
 integer(ik):: i 
 integer(ik):: sum_total_natt_disp
 integer(ik):: sum_total_nacc_disp
    
    sum_total_natt_disp = sum(total_natt(1,:))
    sum_total_nacc_disp = sum(total_nacc(1,:))
    
    Usum = Usum/real(npart,rk)  
    
    !average of energy of all sample blocks 
    do i = 1, nsamp_blocks  
        U_average = U_average + Usum(i) / real(total_att(i),rk) 
    end do 
    U_average = U_average / real(nsamp_blocks,rk)  
    
    !standard deviation of the mean 
    do i = 1, nsamp_blocks
        stdev = stdev + (Usum(i) / real(total_att(i),rk) - U_average)**2.0_rk
    end do 
    stdev = stdev / (nsamp_blocks - 1_ik) 
    stdev = sqrt(stdev)                                 ! standard deviation 
    stdev = stdev/sqrt(real(nsamp_blocks,rk))     ! standard deviation of the mean 
        
    OPEN (res_file_num,FILE = './output_file/'//ResultsFile) 
    
        WRITE (res_file_num,*) 'U_average*beta/N  &   st. dev. of the mean'
        WRITE (res_file_num,'(2f18.10)') U_average, stdev  
        
        WRITE (res_file_num,*) ''
        
        WRITE (res_file_num,*) 'block,   U_average*beta/N per block'
        DO i = 1,  nsamp_blocks
            WRITE (res_file_num,'(i6,f18.10)')  i,  Usum(i) / real(total_att(i),rk) 
        END DO 
        
        WRITE (res_file_num,*) ''
        WRITE (res_file_num,*) 'U*beta/N: last configuration '
        WRITE (res_file_num,"(F10.4)") Utot/real(npart,rk) 
        WRITE (res_file_num,*) ''
        
        WRITE (res_file_num,*) 'accepted and attempted moves from equil. and prod. phase'
        
        WRITE (res_file_num,"(2I10,A14,F10.2,A)")   total_nacc(2,3), total_natt(2,3)   , 'rotation   '  &
                                    , real(total_nacc(2,3))/real(total_natt(2,3))*100.0_rk,   ' % acceptance '
        WRITE (res_file_num,"(2I10,A14,F10.2,A)")  sum_total_nacc_disp, sum_total_natt_disp  , 'translation' &
                                    , real(sum_total_nacc_disp)/real(sum_total_natt_disp)*100.0_rk, ' % acceptance '
       
        WRITE (res_file_num,"(A10,A,F10.2,A)") '','            ---> hs     ',real(total_nacc(1,1))  &
                                                            / real(total_natt(1,1))*100.0_rk, ' % acceptance '
        WRITE (res_file_num,"(A10,A,F10.2,A)") '','            ---> ions   ',real(total_nacc(1,2))  &
                                                            /real(total_natt(1,2))*100.0_rk, ' % acceptance '
        WRITE (res_file_num,"(A10,A,F10.2,A)") '','            ---> dipoles',real(total_nacc(1,3))  &
                                                            /real(total_natt(1,3))*100.0_rk, ' % acceptance '
       
        WRITE (res_file_num,"(2I10,A14,F10.2,A)")   total_accepted_swap, total_attempted_swap     , 'particle swap'  &
                                    , real(total_accepted_swap)/real(total_attempted_swap)*100.0_rk,   ' % acceptance '
                                    
                                             
        WRITE (res_file_num,*) ''
        
        WRITE (res_file_num, *) 'max. translation and rotation (hs,ions,dipoles)'
        WRITE (res_file_num, "(4F10.6)")  max_disp(1), max_disp(2), max_disp(3)        
        WRITE (res_file_num, "(4F10.6)")  max_rot(1), max_rot(2), max_rot(3) 
!         "(1x,A, F10.3, A)" 
        WRITE (res_file_num,*) ''
    
        WRITE (res_file_num,*) 'input data' 
        WRITE (res_file_num,*) nequil,  '   number of equil. cycles'
        WRITE (res_file_num,*) nprod ,  '   number of prod. cycles'
        WRITE (res_file_num,*) nsamp ,  '   sampling number'
        WRITE (res_file_num,*) ''
        WRITE (res_file_num,*) nhs      ,  '   number of hs particles'
        WRITE (res_file_num,*) nions    ,  '   number of ions'
        WRITE (res_file_num,*) ndipoles ,  '   number of dipoles'
        WRITE (res_file_num,*) npart    ,  '   number of particles'
        WRITE (res_file_num,*) ''
        WRITE (res_file_num,"(F10.3,A)") lambda   ,  '   lambda '
        WRITE (res_file_num,"(F10.3,A)") alpha    ,  '   alpha'  
        WRITE (res_file_num,"(F10.3,A)") Temp     ,  '   temperature'
        WRITE (res_file_num,"(F10.5,A)") density  ,  '   density '  
        WRITE (res_file_num,"(F10.5,A)") boxlength,  '   boxlength ' 
        WRITE (res_file_num,"(F10.5,A)") real(nions,rk)/real(npart,rk)      ,  '   fraction'
!         WRITE (res_file_num,*) sigma_list(:)    ,  '   sigma'
         
        
    CLOSE (res_file_num)   
    
     
     
end subroutine store_results 

!--------------------------------------------------------------------------------------------------
! initialize quantities needed for sampling
!-------------------------------------------------------------------------------------------------- 
subroutine initialize_sampling (nsamp_blocks, Usum, total_att)
!  use displace_particle_mod, only: displace_particle 
!  use rotate_particle_mod, only: rotate_particle 
 
 integer(ik), intent(in)                         :: nsamp_blocks
 real(rk), dimension(nsamp_blocks), intent(inout):: Usum
 real(rk), dimension(nsamp_blocks), intent(inout):: total_att

    Usum(:)      = 0.0_rk 
    total_att(:) = 0_ik
       

end subroutine initialize_sampling 


!--------------------------------------------------------------------------------------------------
! wirting particle positions 
!--------------------------------------------------------------------------------------------------
subroutine store_part_config !(iout)  
integer :: i


    OPEN (pt_file_num,FILE = './output_file/'//FinalConf_PartPos) 
    
    WRITE (pt_file_num, *) 'positions (x,y,z) of', npart, 'particles' 

    !hs-particle configuration 
    DO i = 1, npart 
       WRITE (pt_file_num, *) particle_positions(1,i),particle_positions(2,i),particle_positions(3,i)
    END DO
    
    CLOSE (pt_file_num)  
    

  RETURN
end subroutine store_part_config
 
!--------------------------------------------------------------------------------------------------
! wirting charge positions 
!--------------------------------------------------------------------------------------------------
subroutine store_ch_config!(iout)  
integer :: i


    OPEN (ch_file_num,FILE = './output_file/'//FinalConf_ChargePos) 
    
    WRITE (ch_file_num, *) 'positions (x,y,z) of', ncharges, 'charges (ions and dipole charges)' 

    !hs-particle configuration 
    DO i = 1, ncharges
       WRITE (ch_file_num, *) charge_positions(1,i), charge_positions(2,i), charge_positions(3,i)
    END DO
    
    CLOSE (ch_file_num)  

  RETURN
end subroutine store_ch_config



end module rw_files_mod











