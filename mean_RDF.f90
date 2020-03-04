!> module containing the mean radial distribution function 
module mean_RDF_mod
  use working_prec_mod, only: rk, ik
  use globals_mod 
  implicit none
  private
   
   public:: mean_RDF
contains
 
  subroutine mean_RDF ( g_RDF ) 
    integer(ik):: ref_part 
    integer(ik):: n_part_k = 0
    real(rk):: dis_1 
    real(rk):: dis_max
    real(rk), allocatable:: dis_k(:), dis_k2(:) 
    real(rk), intent(out):: g_RDF(intervalls) 
    real(rk):: dis_kmin
    real(rk):: d2
    integer(ik):: i, k
     
    allocate(dis_k(0:intervalls))
    allocate(dis_k2(0:intervalls))
    
    ref_part = 1_ik  
    
    dis_1 =  sigma_list(ref_part)
    dis_max = 4._rk * dis_1
    
    do k = 0, intervalls
        dis_k(k) = ( dis_1**2._rk + ( dis_max**2._rk - dis_1**2._rk ) * real(k,rk)/real(intervalls,rk) )**0.5_rk 
        dis_k2(k) = dis_k(k)*dis_k(k) 
    end do 
     
    do k = 1, intervalls 
      do i = 2, npart
        d2 = squared_distance (particle_positions(:,ref_part), particle_positions(:,i))
        if ( dis_k2(k-1_ik) .lt. d2 .and. d2 .lt. dis_k2(k) )    n_part_k = n_part_k + 1_ik
      end do 
      g_RDF(k) = 3._rk*boxlength**3._rk*real(n_part_k,rk) / ( 2._rk*PI* ( dis_k(k)**3._rk - dis_k(k-1_ik)**3._rk ) * real(npart,rk) * ( real(npart,rk) - 1._rk )  )
      n_part_k = 0_ik
    end do 
     
     
    deallocate(dis_k) 
    deallocate(dis_k2) 
    
  end subroutine mean_RDF

end module mean_RDF_mod









