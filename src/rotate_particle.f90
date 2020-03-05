!> \brief rotating a randomly choosen particle  
!> \author M.Theiss 
module rotate_particle_mod
  use working_prec_mod, only: rk, ik
  use basic_parameters_mod
  use sim_parameters_mod
  use globals_mod 
  use ES_Fourier_mod
  use Energy_mod
  use utilities_mod, only: outfold 
  use realmolecule_mod , only: realmolecule_move
  implicit none
  private

  public :: rotate_particle  
   
contains 

!---------------------------------------------------------------------------------------------------------------------
!> rotating a randomly choosen particle 
!---------------------------------------------------------------------------------------------------------------------
  subroutine rotate_particle ( dUREAL, dUFOURIER, ener_delt)
         
    real(rk) , intent(out) :: dUREAL, dUFOURIER !< delta in energy (real and fourier contribution)   
    real(rk) , intent(out) :: ener_delt !< delta in energy
    integer(ik) :: i, ii !< the infamous counting variables
    integer(ik) :: pid !< id of selected particle
    integer(ik) :: chargeid(2) !< id(s) of ion (dipole) charge(s)
    integer(ik) :: stion, endion, Nc
    real(rk) :: ran, acceptance  !< random numbers
    real(rk) :: ener_old, ener_new !< old/new energy values
    real(rk) :: charge_pos_old_x(2), charge_pos_old_y(2), charge_pos_old_z(2) !< dummies for old positions of charges ion/dipole 
    complex(rk), allocatable, dimension(:,:):: EXPX_OLD, EXPY_OLD, EXPZ_OLD
    complex(rk), dimension(Nkvec) :: SUMQEXPV_OLD
    real(rk):: UREAL, UREALo, USCH,  USDIP, UFOURIER, USURF  !< energy  in ewald sum (ES)
    
    !Rotation
    integer     :: axis 
    real(rk)    :: dtheta, angle, ax
    real(rk)    :: xcom, ycom, zcom !centre of mass <-> particle_positions 
    real(rk),dimension(2,2)    :: M
    real(rk),dimension(2)      :: T
    
    ener_delt = 0.0_rk
    ener_new = 0.0_rk
    ener_old = 0.0_rk
    dUREAL  = 0.0_rk
    dUFOURIER = 0.0_rk 
    
    ! pick random dipole
    call random_number (ran)
    pid = int(ran * real(npart, rk), ik) + 1_ik
     
    ! increase count
    natt(2,species_list(pid)) = natt(2,species_list(pid)) + 1_ik 
    
    !save old positions running_ener
    xcom = particle_positions(1,pid)
    ycom = particle_positions(2,pid)
    zcom = particle_positions(3,pid)
    
    Nc = LENGTH(pid)
    chargeid = return_charge_index (pid) 
    stion = chargeid(1) 
    endion = chargeid(2) 
    
    if (Nc /= 2) then 
      nacc(2,species_list(pid)) = nacc(2,species_list(pid)) + 1_ik
      return 
    end if 
    
    allocate (EXPX_OLD(0:Kmax,1:Nc))
    allocate (EXPY_OLD(-Kmax:Kmax,1:Nc))
    allocate (EXPZ_OLD(-Kmax:Kmax,1:Nc))
    
    call RealMolecule_move(  stion, endion,  UREALo)
    ener_old = UREALo !- USCH -  USDIP + UFOURIERo + USURF
      
    EXPX_OLD(:,1:Nc) = EXPX(:,stion:endion)
    EXPY_OLD(:,1:Nc) = EXPY(:,stion:endion)
    EXPZ_OLD(:,1:Nc) = EXPZ(:,stion:endion) 
    SUMQEXPV_OLD(:) = SUMQEXPV(:)
      
    charge_pos_old_x(1:Nc)  = charge_positions(1,stion:endion) 
    charge_pos_old_y(1:Nc)  = charge_positions(2,stion:endion)
    charge_pos_old_z(1:Nc)  = charge_positions(3,stion:endion)
    
    ! reconstructing position of dipole charges 
    call outfold ( pid, stion, endion ) 
     
    ! give particle random displacement
    call random_number(angle)
    dtheta = ( 2.0 * angle - 1.0 ) * Pi * max_rot(species_list(pid))
    
    M = reshape( (/ cos( dtheta ), -sin( dtheta ), sin( dtheta ), cos( dtheta ) /), &
                            (/ 2, 2 /) )

    call random_number(ax)
    axis = int( 3.0 * ax ) + 1
    
    Select Case ( axis )
            
        case ( 1 )              ! Rotation in y-z plane.
                
                ii = 1_ik
                do i = stion, endion 
                
                        T(1) = charge_positions(2,i) - ycom
                        T(2) = charge_positions(3,i) - zcom

                        T = matmul( M, T )
                            
                        charge_positions(2,i) = T(1) + ycom 
                        charge_positions(3,i) = T(2) + zcom 
                        charge_positions(1,i) = charge_pos_old_x(ii)
                        ii = ii + 1_ik
                end do
            

        case ( 2 )              ! Rotation in x-z plane.

                ii = 1_ik
                do i = stion, endion 
                        
                        T(1) = charge_positions(3,i) - zcom
                        T(2) = charge_positions(1,i) - xcom

                        T = matmul( M, T )
                            
                        charge_positions(3,i) = T(1) + zcom 
                        charge_positions(1,i) = T(2) + xcom 
                        charge_positions(2,i) = charge_pos_old_y(ii)
                        ii = ii + 1_ik
        
                end do
            
        case ( 3 )              ! Rotation in x-y plane.

                ii = 1_ik
                do i = stion, endion 
                
                        T(1) = charge_positions(1,i) - xcom
                        T(2) = charge_positions(2,i) - ycom
                        
                        T = matmul( M, T ) 
                        
                        charge_positions(1,i) = T(1) + xcom 
                        charge_positions(2,i) = T(2) + ycom 
                        charge_positions(3,i) = charge_pos_old_z(ii)
                        ii = ii + 1_ik
                        
                end do
    
    end select 
    
    ! periodic boundaries:
    do ii = stion, endion 
      do i = 1,3  
        if (charge_positions(i,ii) > boxlength) charge_positions(i,ii)= charge_positions(i,ii) - boxlength
        if (charge_positions(i,ii) < 0.0_rk)    charge_positions(i,ii)= charge_positions(i,ii) + boxlength  
      end do  
    end do 
     
    call DispRot_Energy ( pid, Nc, stion, endion, UREAL, USCH,  USDIP, UFOURIER, USURF)  
    ener_new = UREAL + UFOURIER !- USCH -  USDIP + USURF 
     
    ener_delt = ener_new - ener_old
    dUREAL = UREAL-UREALo
    dUFOURIER = UFOURIER !- UFOURIERo   
    
    ! accept or reject move
    call random_number(acceptance)
    if (acceptance < exp( - lambda * ener_delt)) then
        nacc(2,species_list(pid)) = nacc(2,species_list(pid)) + 1_ik    
!         print*, 'Accepted'; pause
    else
        ! set back  
        charge_positions(1,stion:endion) = charge_pos_old_x(1:Nc)  
        charge_positions(2,stion:endion) = charge_pos_old_y(1:Nc)  
        charge_positions(3,stion:endion) = charge_pos_old_z(1:Nc)  
        EXPX(:,stion:endion) = EXPX_OLD(:,1:Nc)
        EXPY(:,stion:endion) = EXPY_OLD(:,1:Nc)
        EXPZ(:,stion:endion) = EXPZ_OLD(:,1:Nc) 
        SUMQEXPV(:) = SUMQEXPV_OLD(:)
        ener_delt = 0.0_rk
        dUREAL  = 0.0_rk
        dUFOURIER = 0.0_rk
!        print*, 'Rejected'; pause
    end if 
    
    deallocate(EXPX_OLD)
    deallocate(EXPY_OLD)
    deallocate(EXPZ_OLD)
    
  end subroutine rotate_particle 
end module rotate_particle_mod
