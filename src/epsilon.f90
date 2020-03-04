module eps_mod 
  use working_prec_mod, only: rk, ik
  use globals_mod 
  use utilities_mod
  implicit none
  private
   
   public:: eps  
contains
 
  subroutine eps(moves)  
  integer(ik), intent(in) :: moves 
  
  integer(ik) :: i, ii 
  integer(ik) :: pid, stion, endion, lenion 
  integer(ik), dimension(2) :: chargeid  
  
  real(rk)    :: len_vec, totDM
  real(rk), dimension(3) :: eps_tmp, vec 



  eps_tmp(:) = 0.0_rk

  do pid = 1, npart 

    lenion = LENGTH(pid) 
    !considering only dipole-particles 
    if (lenion .ne. 2) cycle 
       
    chargeid = return_charge_index (pid) 
    stion  = chargeid(1) 
    endion = stion + lenion - 1_ik
    
    !get the charges back to the particles 
    call outfold ( pid, stion, endion ) 
    
    vec(:)  = charge_positions(:,stion) - charge_positions(:,endion)       
    len_vec = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2) 
    !unit-vector 
    vec(:)  = vec(:) / len_vec  
    
    !considering only pure components (otherwise: multiply with dipol-strength!) 
    eps_tmp(:) = eps_tmp(:) + vec(:) 

    ! periodic boundaries:
    do ii = stion, endion
      do i = 1,3  
        if (charge_positions(i,ii) > boxlength) charge_positions(i,ii)= charge_positions(i,ii) - boxlength
        if (charge_positions(i,ii) < 0.0_rk)    charge_positions(i,ii)= charge_positions(i,ii) + boxlength  
      end do  
    end do 

  end do 
  
  !total dipole moment 
  totDM = eps_tmp(1)**2 + eps_tmp(2)**2 + eps_tmp(3)**2 

  OPEN (re_file_num+1,FILE = './output_file/'//trim(SimName)//'.totSDM', position='append') 
    WRITE(re_file_num+1,*) moves, totDM, eps_tmp(1), eps_tmp(2), eps_tmp(3)  
    !WRITE(re_file_num+1,"(I20,4F20.5)") moves, totDM, eps_tmp(1), eps_tmp(2), eps_tmp(3)  
  CLOSE (re_file_num+1)
   
end subroutine eps  

end module eps_mod 









