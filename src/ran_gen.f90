module ran_gen_mod
  use working_prec_mod, only: ik, rk
  implicit none
  private
  public :: init_random_seed

  logical, public, protected :: random_seed_initialized = .false.
contains

  ! Initialisation of the RNG !
  subroutine init_random_seed()
    implicit none
    integer(ik), allocatable :: seed(:)
    integer(ik) :: n, un, istat

    call RANDOM_SEED(size = n)
    allocate(seed(n))
    ! The OS provides a random number generator:
    open(newunit=un, file="/dev/urandom", access="stream", &
      form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      print*, 'ran_gen_mod:init_random_seed(): ERROR while opening File urandom: ',istat
      stop
    end if

    ! uncomment this to generate same set of rnd for every run
!     seed = 1_ik

    call RANDOM_SEED(put=seed)

    ! set logical
    random_seed_initialized = .true.
    !    if (VERBOSE) print*, 'ran_gen_mod:init_random_seed(): Initialization successfull '
  end subroutine init_random_seed
end module
