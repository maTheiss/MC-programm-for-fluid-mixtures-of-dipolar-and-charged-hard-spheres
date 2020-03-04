Module adjust_displacement_mod
use working_prec_mod, only: rk, ik  
use globals_mod     , only: boxlength
use sim_parameters_mod
  implicit none
  private

  public :: adjust_disp 
  public :: adjust_rot  
  public :: adjust_max_dispRot    

contains

!-------------------------------------------------------------------------------------------------------
!     Adjusts Maximum Displacement Such That 50% Of The
!     Movels Will Be Accepted 
!       -> a high ratio means small displacement moves so that local min. remain longer!!?
!-------------------------------------------------------------------------------------------------------
subroutine adjust_max_dispRot
INTEGER ::  i   
  REAL(rk) :: SimRatio, TarRatio(2)  
    
  
TarRatio(1) = 0.5
TarRatio(2) = 0.5


do i = 1, 3  !Nsp

    if( NATT(1,i) /= 0 ) then

            SimRatio = real( NACC(1,i) ) / real( NATT(1,i) )

            if( SimRatio > TarRatio(1) ) max_disp(i) = max_disp(i) * 1.05    
            if( SimRatio < TarRatio(1) ) max_disp(i) = max_disp(i) * 0.95    

            max_disp(i) = min(max_disp(i),0.5_rk)       

    end if
end do 

if( NATT(2,3) /= 0 ) then
                
        SimRatio = real( NACC(2,3) ) / real( NATT(2,3) )

        if( SimRatio > TarRatio(2) ) max_rot = max_rot * 1.05             
        if( SimRatio < TarRatio(2) ) max_rot = max_rot * 0.95

        max_rot = min( max_rot, 1.0_rk )            

end if
 
NATT = 0_ik
NACC = 0_ik

end subroutine adjust_max_dispRot
 
!-------------------------------------------------------------------------------------------------------
!     Adjusts Maximum Displacement Such That 50% Of The
!     Movels Will Be Accepted 
!       -> a high ratio means small displacement moves so that local min. remain longer!!?
!-------------------------------------------------------------------------------------------------------
SUBROUTINE adjust_disp(attemp, nacc, dr)

  INTEGER, INTENT(IN)                      :: attemp
  INTEGER, INTENT(IN)                      :: nacc
  REAL(rk), INTENT(IN OUT)                 :: dr
  INTEGER :: attempp, naccp
  REAL(rk) :: dro, frac, hbox

  SAVE naccp, attempp
  DATA naccp/0/
  DATA attempp/0/


  !  Attemp (Input)  Number Of Attemps That Have Been Performed
  !                  To Displace A Particle
  !  Nacc   (Input)  Number Of Successful Attemps To
  !                  Displace A Particle
  !  Dr     (Output) New Maximum Displacement


  hbox = 0.5_rk*boxlength 

  IF (attemp == 0.OR.attempp >= attemp) THEN
    naccp = nacc
    attempp = attemp
  ELSE
    frac = DBLE(nacc-naccp)/DBLE(attemp-attempp)
    dro  = dr
    dr   = dr*ABS(frac/0.5D0)
    
  !        ---Limit The Change:
    
    IF (dr/dro > 1.5D0) dr = dro*1.5D0
    IF (dr/dro < 0.5D0) dr = dro*0.5D0
    IF (dr > hbox/2.d0) dr = hbox/2.d0
! !     WRITE (6, 99001) dr, dro, frac, attemp - attempp, nacc - naccp
    
  !        ---Store Nacc And Attemp For Next Use
    
    naccp = nacc
    attempp = attemp
  END IF
  RETURN
! !   99001 FORMAT (' Max. Displ. Set To : ', f6.3, ' (Old : ', f6.3, ')', /,  &
! !       ' Frac. Acc.: ', f5.2, ' Attempts: ', i7, ' Succes: ', i7)

END SUBROUTINE adjust_disp


SUBROUTINE adjust_rot(attemp, nacc, dr)

  INTEGER, INTENT(IN)                      :: attemp
  INTEGER, INTENT(IN)                      :: nacc
  REAL(rk), INTENT(IN OUT)                 :: dr
  INTEGER :: attempp, naccp
  REAL(rk) :: dro, frac, hbox

  SAVE naccp, attempp
  DATA naccp/0/
  DATA attempp/0/

  !     Adjusts Maximum Displacement Such That 50% Of The
  !     Movels Will Be Accepted

  !  Attemp (Input)  Number Of Attemps That Have Been Performed
  !                  To Displace A Particle
  !  Nacc   (Input)  Number Of Successful Attemps To
  !                  Displace A Particle
  !  Dr     (Output) New Maximum Displacement


  hbox = 0.5_rk*boxlength 

  IF (attemp == 0.OR.attempp >= attemp) THEN
    naccp = nacc
    attempp = attemp
  ELSE
    frac = DBLE(nacc-naccp)/DBLE(attemp-attempp)
    dro  = dr
    dr   = dr*ABS(frac/0.5D0)
    
  !        ---Limit The Change:
    
    IF (dr/dro > 1.5D0) dr = dro*1.5D0
    IF (dr/dro < 0.5D0) dr = dro*0.5D0
    IF (dr > hbox/2.d0) dr = hbox/2.d0
! !     WRITE (6, 99001) dr, dro, frac, attemp - attempp, nacc - naccp
    
  !        ---Store Nacc And Attemp For Next Use
    
    naccp = nacc
    attempp = attemp
  END IF
  RETURN
! !   99001 FORMAT (' Max. Displ. Set To : ', f6.3, ' (Old : ', f6.3, ')', /,  &
! !       ' Frac. Acc.: ', f5.2, ' Attempts: ', i7, ' Succes: ', i7)

END SUBROUTINE adjust_rot

END MODULE adjust_displacement_mod