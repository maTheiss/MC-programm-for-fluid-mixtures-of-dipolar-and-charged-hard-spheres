!> \brief NVT main routine  
!> \author M.Theiss 
module nvt_mod
  use working_prec_mod, only: rk, ik
  use sim_parameters_mod
  use globals_mod      
  use displace_particle_mod, only: displace_particle 
  use rotate_particle_mod, only: rotate_particle 
!  use mean_RDF_mod, only: mean_RDF
  use ES_Fourier_mod
  use Energy_mod
  use rw_files_mod    
  use adjust_displacement_mod, only: adjust_max_dispRot, adjust_disp, adjust_rot 
  use eps_mod  
  implicit none
  private

  public :: nvt

  logical     :: is_production 

contains

  !> NVT main routine  
  subroutine nvt (nequil, nprod)
    integer(ik) , intent(inout) :: nequil, nprod
    integer(ik) :: ncycles
    integer(ik) :: iii, icycl, i
    integer(ik) :: cyc = 1_ik
    real(rk)    :: which_move
    real(rk)    :: ener_delt
    real(rk)    :: running_ener = 0.0_rk
    real(rk)    :: energy = 0.0_rk 
   
!     real(rk),dimension(intervalls)   :: gr 
!     integer(ik):: mean = 0
    
    real(rk):: UREAL, USCH,  USDIP, UFOURIER, USURF  !< energy  in ewald sum (ES)
    
    real(rk):: UREAL1, USCH1,  USDIP1, UFOURIER1, USURF1   , Utot, Utot1, dUReal, dUFourier
    
    integer(ik)  :: block = 1_ik
    real(rk), allocatable, dimension(:):: Usum
    real(rk), allocatable, dimension(:):: U_cyc
    integer(ik), allocatable, dimension(:):: total_att 
    integer(ik) :: moves = 0_ik
    
    integer(ik):: sum_total_att
    real(rk)   :: U_average 
    
!      gr (:) = 0._rk 
          
     !opens file for continuous writing 
     OPEN (plt_file_num,FILE = './output_file/'//UavFile) 
     OPEN (re_file_num,FILE = './output_file/'//RunnEnFile) 
     OPEN (re_file_num+1,FILE = './output_file/'//trim(SimName)//'.totSDM', status='replace')
     CLOSE (re_file_num+1)

     call Fourier_Setup( ) 
     
     allocate( SUMQEXPV(Nkvec) ) 
     allocate( Usum(nsamp_blocks) ) 
     allocate( U_cyc(max(nequil,nprod)) )
     allocate( total_att(nsamp_blocks) ) 
     Usum(:)  = 0.0_rk 
     U_cyc(:) = 0.0_rk
     total_att(:) = 0_ik 
    
    total_attempted_swap = 0 
    total_accepted_swap  = 0
    
    ! production / equilibration
    do iii = 1,2
        
      print*, ''
      if (iii == 1) then
        ncycles = nequil
        is_production = .false.
        print*, 'Starting equilibration!'
      else
        ncycles = nprod
        is_production = .true.
        print*, 'Starting production!'
      end if
      
      ! calculate initial system energy 
      call TotalEnergy( UREAL1, USCH1,  USDIP1, UFOURIER1, USURF1, Utot1) 
      energy = (UREAL1 - USCH1 - USDIP1 + UFOURIER1 + USURF1 )  
      running_ener = Utot1 !energy     
      
      ! intialize the subroutine that adjusts the maximum displacement 
      call adjust_max_dispRot 
        
      ! initialize Usum and counter 
      if (is_production) then  
        Usum(:)     = 0.0_rk
        total_att(:)= 0_ik 
        block       = 1_ik
        U_cyc(:)    = 0.0_rk
        cyc         = 1_ik
      end if 
      
      ! cycle loop
      do icycl = 1, ncycles 
       
        ! move loop
        do i = 1, 1000
          ! decide which move? 
          call random_number(which_move) 
          if (which_move .le. 0.5)  then
            call displace_particle ( dUreal, dUFourier, ener_delt)  
          else 
            call rotate_particle ( dUreal, dUFourier, ener_delt)    
          end if              
          running_ener = running_ener + ener_delt
          UREAL1 = UREAL1 + dUreal 
          UFOURIER1 = UFOURIER1 + dUFourier
          Usum(block) = Usum(block) + running_ener 
          U_cyc(cyc)  = U_cyc(cyc) + running_ener 
          total_att(block) = total_att(block) + 1_ik    !NOT appropriate for npart == nions!
          moves = moves + 1_ik
        end do 
         
     
!        !calc. radial distribution function 
!         if (mod(icycl,4) == 0 .and. is_production == .true.  ) then
!          print*, 'calc. RDF:'
!          call mean_RDF( gr )
!          mean = mean + 1_ik
!          gr(:) = gr(:) + gr(:)
!         end if
!         gr(:) = gr(:) / real(mean,rk) 
                 
        if ( mod(icycl,nsamp) == 0) then
          print*, 'Total att:', moves 

          !saving average of running energy  
          sum_total_att = sum(total_att(1:block))
          U_average = sum(Usum(1:block))/real(sum_total_att,rk)
          U_cyc(cyc) = U_cyc(cyc)/real(nsamp,rk)/1000.0_rk   
          print*, 'U_av', U_average/real(npart,rk) 
          print*, 'U_cy', U_cyc(cyc)/real(npart,rk)
          call sample(moves, U_average, U_cyc(cyc))
          cyc = cyc + 1_ik 
          
          !saving running energy  
          call store_running_energy(moves, running_ener) 

          !computing and saving total dipole moment
          call eps(moves)

        end if
         
        !adjust maximum displacements  
        ! not in prod. phase (detailed balance)
        if ( mod(icycl,nsamp/10) == 0 .and. iii == 1 ) then   
          total_natt = total_natt + natt
          total_nacc = total_nacc + nacc 
          call adjust_max_dispRot     
        end if 
         
        if ( mod(icycl,ncycles/nsamp_blocks) == 0) then
          !counter for sampling blocks 
          block = block + 1_ik 
        end if 
        
      end do !cycle loop for equil. or prod. 
      
      total_nacc = total_nacc + nacc   
      total_natt = total_natt + natt 
      
    end do ! iii (equil/prod)
              
    ! save final configurations of particle and charge positions 
    call store_part_config 
    call store_ch_config 
    !save results: U_average, ..  
    call store_results(nequil, nprod, nsamp_blocks, Usum(1:nsamp_blocks), total_att(1:nsamp_blocks), Utot) 
 
!    !results for radial distribution function 
!    WRITE (72,*) ' '
!    write (*,*) 'RDF_mean.dat...' 
!    OPEN (72,FILE = './output_file/RDF_mean.dat') 
!     DO i = 1,  256      
!       WRITE (72,'(i6,10f18.10)')  i,  gr(i)
!     END DO 
!    WRITE (72,*) ' '
!    CLOSE (72)   

    CLOSE (plt_file_num)
    CLOSE (re_file_num)

    deallocate (particle_positions)
    deallocate (charge_positions)
    deallocate (sigma_list)
    deallocate (charge_list)
    deallocate (species_list)
    deallocate (LENGTH)
    deallocate ( SUMQEXPV )  
    deallocate ( CONST ) 
    deallocate ( KX ) 
    deallocate ( KY ) 
    deallocate ( KZ ) 
    deallocate ( EXPX )
    deallocate ( EXPY ) 
    deallocate ( EXPZ )  
    deallocate ( Usum ) 
    deallocate ( total_att ) 
    
    write(*,*)
    write(*,*) trim(SimName), ' complete'

  end subroutine nvt
end module nvt_mod



