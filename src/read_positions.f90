!--------------------------------------------------------------------------------------------------
! reading hs and charge positions 
!--------------------------------------------------------------------------------------------------
subroutine read_positions
 use working_prec_mod 
 use globals_mod 
IMPLICIT NONE
character (len=40) ::  dum_text
integer(ik)        :: i 
 
    print*, ''
    print*, 'configurations taken from:'
    print*, FinalConf_PartPos
    print*, FinalConf_ChargePos
    print*,'' 
    
    !reading hs-positionsv
    OPEN (pt_file_num,FILE = './output_file/'//FinalConf_PartPos) 
    READ (pt_file_num,*) dum_text  !boxlength
     do i = 1, npart
        READ(pt_file_num,*) particle_positions(1,i), particle_positions(2,i), particle_positions(3,i)
     end do 
    CLOSE (pt_file_num)

    !reading ch-positions 
    OPEN (ch_file_num,FILE = './output_file/'//FinalConf_ChargePos) 
    READ (ch_file_num,*) dum_text  !boxlength
     do i = 1, ncharges
        READ(ch_file_num,*) charge_positions(1,i), charge_positions(2,i), charge_positions(3,i)
     end do
    CLOSE (ch_file_num)

  RETURN
end subroutine read_positions

 