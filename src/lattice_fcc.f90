SUBROUTINE lattice_fcc
 use globals_mod 
 use working_prec_mod

! Code converted using TO_F90 by Alan Miller
! Date: 2013-10-02  Time: 09:42:19

!     Place `Npart' Particles On A Simple Cubic
!     Lattice With Density 'Rho'
 
IMPLICIT NONE

integer(ik) ::i, ii, j, k, itel, n  
real(rk) ::dx, dy, dz, del
integer(ik) :: chargeid(2)  
real(rk):: ran 
real(rk):: pos(3,npart) 
real::    A(npart), B(npart)

pos(:,:) = 0_rk

! output
print*, 'Set startconfig -> cubic fcc!'
 
n = int((npart/4.0_rk)**(1.0_rk/3.0_rk)-1e-6)+1_ik 

if (n == 0) n = 1_ik
del = 1.0_rk/REAL(n,rk)
itel = 0_ik
dx = -del
do i = 1, n
      dx = dx + del
      dy = -del
      do j = 1, n
        dy = dy + del
        dz = -del
        do k = 1, n
          dz = dz + del
          if (itel < npart) then
            itel = itel + 1_ik
 
            pos(1,itel) = dx
            pos(2,itel) = dy 
            pos(3,itel) = dz 

            if (itel < npart .and. k<=n .and. j<=n) then
              itel = itel + 1
               
              pos(1,itel) = dx
              pos(2,itel) = dy+0.5_rk*del
              pos(3,itel) = dz+0.5_rk*del
              
            end if
            
            if (itel < npart .and. k<=n .and. i<=n) then
              itel = itel + 1
               
              pos(1,itel) = dx+0.5_rk*del
              pos(2,itel) = dy 
              pos(3,itel) = dz+0.5_rk*del

            end if
            
            if (itel < npart .and. j<=n .and. i<=n) then
              itel = itel + 1
               
              pos(1,itel) = dx+0.5_rk*del
              pos(2,itel) = dy+0.5_rk*del
              pos(3,itel) = dz 

            end if
            
          end if        
          
        end do
      end do
end do
  
 call sort(npart,A,B)
 
  do i = 1, npart
   do ii = 1, npart 
    if ( B(i) == A(ii) ) then 
        particle_positions(:,i) = pos(:,ii)   
    end if 
   end do 
  end do 
  
 particle_positions(:,:) = particle_positions(:,:) * boxlength 
  
  do i = 1, npart
    if (is_overlapping(i)) then 
        print*, 'overlapping in lattice_fcc'
        stop
    end if 
  end do 
  
  if (nions /= 0_ik) then 
    do i = 1, nions
        charge_positions(:,i) = particle_positions(:,i+nhs) 
    end do 
   end if  

   if (ndipoles /= 0_ik) then 
    do ii = 1, npart 
      if ( ii .gt. nions + nhs ) then 
        chargeid = return_charge_index (ii)
        call random_number(ran)

        charge_positions(1,chargeid(1)) =  particle_positions(1,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * cos(ran*2._rk*PI-PI)
        charge_positions(2,chargeid(1)) =  particle_positions(2,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * sin(ran*2._rk*PI-PI)
        charge_positions(3,chargeid(1)) =  particle_positions(3,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * cos(ran*PI)  
        charge_positions(1,chargeid(2)) =  particle_positions(1,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * cos(ran*2._rk*PI-PI)
        charge_positions(2,chargeid(2)) =  particle_positions(2,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * sin(ran*2._rk*PI-PI)
        charge_positions(3,chargeid(2)) =  particle_positions(3,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * cos(ran*PI)  
                                                
     end if
    end do 
   end if   
   
   ! periodic boundaries for dipole charges only
  if (ndipoles /= 0_ik) then 
   do ii = 1, npart 
    if ( ii .gt. nions + nhs ) then 
      chargeid = return_charge_index (ii)   
      do i = 1,3  
        if(charge_positions(i,chargeid(1)).gt. boxlength)charge_positions(i,chargeid(1))=charge_positions(i,chargeid(1)) - boxlength
        if(charge_positions(i,chargeid(1)).lt. 0.0_rk) charge_positions(i,chargeid(1))= charge_positions(i,chargeid(1)) + boxlength 
        if(charge_positions(i,chargeid(2)).gt. boxlength)charge_positions(i,chargeid(2))=charge_positions(i,chargeid(2)) - boxlength 
        if(charge_positions(i,chargeid(2)).lt. 0.0_rk) charge_positions(i,chargeid(2))= charge_positions(i,chargeid(2)) + boxlength  
      end do  
    end if 
   end do 
  end if 
   
           
    
WRITE (6, 99001) itel 
RETURN
99001 FORMAT (' Initialisation On Lattice: ', /, i10,  &
    ' Particles Placed On A Lattice')
END SUBROUTINE lattice_fcc



! !-------------------------------------------------------------------------------------------------
! !   Sorting an array with the Heapsort method  
! !------------------------------------------------------------------------------------------------- 
SUBROUTINE SORT(N,A,B)
integer,intent(in) :: N
real,intent(out)::    A(N), B(N)     !Table to be sorted
real, parameter::    MAX_VALUE  = 1000.0  !Maximum value of table
integer:: Count(12) !statt: Count(1)   !Current value of system clock
integer:: i
real:: x

  !initialize random number generator
!  call System_Clock(Count(12)) 
!  call Random_Seed(Put=Count)
  
  !generate random table of numbers (from 0 to 1000)
  do i=1, N
    call Random_Number(x)  !returns random x >= 0 and <1  
    A(i)=MAX_VALUE*x
  end do

!   print *,' '
!   print *,'Table to be sorted:'
!   call TWRIT(N,A)
  B(:) = A(:)
  !call heapsort subroutine
  call HPSORT(N,B)

!   print *,' '
!   print *,'Sorted table (Heapsort method):'
!   call TWRIT(N,B) 
!   print *,' '
!   stop

END SUBROUTINE sort


!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*          RA	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT(N,RA)
  integer, intent(in):: N
  real:: RA(N), RRA
  integer:: I, J, L, IR
  
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase. 
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END



