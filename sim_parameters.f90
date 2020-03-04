
module sim_parameters_mod
use working_prec_mod , only: rk, ik 
use basic_parameters_mod
  implicit none
  
! total number of accepted/attempted displacement(trans/rot) moves per SpID(hs,ion,dip.) 
integer(ik)                                             :: total_accepted_swap, total_attempted_swap
integer(ik), dimension(2,3)                             :: total_natt, total_nacc, natt, nacc 

! EXPX contains the value of exp( i*kx*x ) for a given kx and ion.  
complex(rk), allocatable, dimension(:,:)                :: EXPX
complex(rk), allocatable, dimension(:,:)                :: EXPY, EXPZ

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian.
complex(rk), allocatable, dimension(:)                  :: SUMQEXPV 

                                 
! Nkvec is the number of k-vectors used in the Fourier sum.
! KX, KY, KZ contain the vector identity of the Nkvec vectors.
! CONST contains the constant part of the Fourier summation for a given Nkvec.
integer(ik)                                             :: Nkvec    
real(rk), allocatable, dimension(:)                     :: CONST   
integer(ik), allocatable, dimension(:)                  :: KX, KY, KZ   

! number of sampling blocks 
integer(ik), parameter                                  :: nsamp_blocks = 5

! Names of input/output files 
character (len=40)                                      :: in_file 
character (len=40)                                      :: SimName
character(:), allocatable                               :: ResultsFile
character(:), allocatable                               :: FinalConf_PartPos
character(:), allocatable                               :: FinalConf_ChargePos
character(:), allocatable                               :: UavFile
character(:), allocatable                               :: RunnEnFile

! file numbers for.. 
integer, parameter                                      :: pt_file_num  = 10    !particle positions
integer, parameter                                      :: ch_file_num  = 11    !charge positions 
integer, parameter                                      :: res_file_num = 12    !results file
integer, parameter                                      :: plt_file_num = 13    !plot file
integer, parameter                                      :: re_file_num  = 14    !running energy file

!maximum displacement: translation (hs,ions,dipoles) and rotation
real(rk) :: max_disp(3) 
real(rk) :: max_rot(3)

integer(ik)  :: naccp(4), attempp(4)
  
!for random number
integer:: Seed

end module sim_parameters_mod                            
