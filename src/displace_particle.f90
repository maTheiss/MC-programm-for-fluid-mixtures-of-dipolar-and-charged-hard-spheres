!> \brief displacing a randomly choosen particle  
!> \author M.Theiss 
module displace_particle_mod
  use working_prec_mod, only: rk, ik
  use basic_parameters_mod
  use sim_parameters_mod
  use globals_mod
  use ES_Fourier_mod
  use Energy_mod
  use realmolecule_mod , only: realmolecule_move
  implicit none
  private

  public :: displace_particle !< public subroutine 
   
contains

!---------------------------------------------------------------------------------------------------------------------
!> displacing a randomly choosen particle 
!---------------------------------------------------------------------------------------------------------------------
  subroutine displace_particle ( dUREAL, dUFOURIER, ener_delt)
       
    real(rk) , intent(out) :: dUREAL, dUFOURIER !< delta in energy (real and fourier contribution)  
    real(rk) , intent(out) :: ener_delt !< delta in energy
    integer(ik) :: i  !< the infamous counting variables
    integer(ik) :: pid !< id of selected particle
    integer(ik) :: chargeid(2) !< id(s) of ion (dipole) charge(s)
    integer(ik) :: stion, endion, Nc
    real(rk) :: ran, acceptance, displace_num(3) !< random numbers
    real(rk) :: displace(3) 
    real(rk) :: ener_old, ener_new !< old/new energy values
    real(rk) :: position_old(3) !< dummy for old position of hard sphere
    real(rk) :: charge_position_old(3) !< dummy for old position of charge ion/dipole
    real(rk) :: charge_position_old2(3) !< dummy for old position of charge dipole
!     complex(rk), dimension(0:Kmax, ncharges):: EXPX_OLD
!     complex(rk), dimension(-Kmax:Kmax, ncharges):: EXPY_OLD, EXPZ_OLD
    complex(rk), allocatable, dimension(:,:):: EXPX_OLD, EXPY_OLD, EXPZ_OLD
    complex(rk), dimension(Nkvec) :: SUMQEXPV_OLD
! !     real(rk), intent(inout)     :: SUMQX, SUMQY, SUMQZ  
    real(rk):: UREAL, UREALo, USCH,  USDIP, UFOURIER, USURF  !< energy  in ewald sum (ES)
    
    dUREAL = 0.0_rk
    dUFOURIER = 0.0_rk
    ener_delt = 0.0_rk
    ener_new = 0.0_rk
    ener_old = 0.0_rk
    
    ! pick random particle
    call random_number (ran)
    pid = int(ran * real(npart, rk), ik) + 1_ik
     
    ! increase count 
    natt(1,species_list(pid)) = natt(1,species_list(pid)) + 1_ik 
      
    ! this function returns chargeid(dim=2),
    ! where chargeid(2) = 0 is an ion
    ! else it is a dipole.
    Nc = LENGTH(pid)
    chargeid = return_charge_index (pid) 
    stion = chargeid(1) 
    endion = stion + Nc - 1_ik
    
    allocate (EXPX_OLD(0:Kmax,1:Nc))
    allocate (EXPY_OLD(-Kmax:Kmax,1:Nc))
    allocate (EXPZ_OLD(-Kmax:Kmax,1:Nc))
    
    if (species_list(pid) /= 1_ik) then  
!       call DispRot_Energy ( pid, LENGTH(pid), stion, endion, UREAL, USCH,  USDIP, UFOURIER, USURF)  
      call RealMolecule_move(  stion, endion,  UREALo)
!       call RealInteract_move(  pid,  UREALo) !(etw. langsamer)
!       UREALo = UREALo !* CoulCombo 
      ener_old = UREALo !- USCH -  USDIP + UFOURIERo + USURF
      EXPX_OLD(:,1:Nc) = EXPX(:,stion:endion)
      EXPY_OLD(:,1:Nc) = EXPY(:,stion:endion)
      EXPZ_OLD(:,1:Nc) = EXPZ(:,stion:endion)  
      SUMQEXPV_OLD(:) = SUMQEXPV(:) 
    end if 
    
    ! give particle random displacement
    call random_number (displace_num) 
    displace = (displace_num-0.5_rk)*2.0_rk*max_disp(species_list(pid))*boxlength
    
    position_old(:) = particle_positions (:,pid) 
    particle_positions(:,pid) = position_old(:) + displace
    
    if (is_overlapping(pid)) then 
      particle_positions(:,pid) = position_old(:) 
      ener_delt = 0.0
      dUREAL = 0.0_rk
      dUFOURIER = 0.0_rk
      deallocate(EXPX_OLD)
      deallocate(EXPY_OLD)
      deallocate(EXPZ_OLD)
      RETURN  
    end if
         
    ! periodic boundaries:
    do i = 1,3 
      if (particle_positions(i,pid) > boxlength) particle_positions(i,pid) = particle_positions(i,pid) - boxlength   
      if (particle_positions(i,pid) < 0.0_rk )   particle_positions(i,pid) = particle_positions(i,pid) + boxlength 
    end do
    
       ! is the particle a ion or dipole?
    if (species_list(pid) /= 1_ik) then
 
      charge_position_old(:) = charge_positions(:,stion)
      charge_positions(:,stion) = charge_position_old(:) + displace
      if (chargeid(2) /= 0_ik) then
        charge_position_old2(:) = charge_positions(:,chargeid(2))
        charge_positions(:,chargeid(2)) = charge_position_old2(:) + displace
      end if
       
    ! periodic boundaries for dipole charges only  
    do i = 1,3  
     if (charge_positions(i,stion) > boxlength) charge_positions(i,stion)=charge_positions(i,stion) - boxlength
     if (charge_positions(i,stion) < 0.0_rk)    charge_positions(i,stion)=charge_positions(i,stion) + boxlength
    end do  
    if (chargeid(2) /= 0_ik) then
     do i = 1, 3 
      if (charge_positions(i,chargeid(2)) > boxlength) charge_positions(i,chargeid(2))=charge_positions(i,chargeid(2)) - boxlength 
      if (charge_positions(i,chargeid(2)) < 0.0_rk)    charge_positions(i,chargeid(2))=charge_positions(i,chargeid(2)) + boxlength  
     end do 
    end if 
       
     call DispRot_Energy ( pid, Nc, stion, endion, UREAL, USCH,  USDIP, UFOURIER, USURF)  
     ener_new = UREAL + UFOURIER !- USCH -  USDIP + USURF 
    end if  
    
    ener_delt = ener_new - ener_old
    dUREAL = UREAL-UREALo
    dUFOURIER = UFOURIER !- UFOURIERo 
    
    ! accept or reject move
    call random_number(acceptance)
    if (acceptance < exp( - lambda * ener_delt)) then  
        nacc(1,species_list(pid)) = nacc(1,species_list(pid)) + 1_ik 
        ! do nothing
!         print*, 'Accepted'!, ener_delt/real(npart,rk)/Temp   !; pause
    else
        ! set back 
        particle_positions(:,pid) = position_old(:)
        if (species_list(pid) /= 1_ik) then 
            charge_positions(:,chargeid(1)) = charge_position_old(:)
            if (chargeid(2) /= 0_ik) charge_positions(:,chargeid(2)) = charge_position_old2(:)
        end if 
        EXPX(:,stion:endion) = EXPX_OLD(:,1:Nc)
        EXPY(:,stion:endion) = EXPY_OLD(:,1:Nc)
        EXPZ(:,stion:endion) = EXPZ_OLD(:,1:Nc)  
        SUMQEXPV(:) = SUMQEXPV_OLD(:)
        ener_delt = 0.0_rk
        dUREAL = 0.0_rk
        dUFOURIER = 0.0_rk 
!        print*, 'Rejected'!; pause 
    end if
    
    deallocate(EXPX_OLD)
    deallocate(EXPY_OLD)
    deallocate(EXPZ_OLD)
      
  end subroutine displace_particle
   
end module displace_particle_mod
