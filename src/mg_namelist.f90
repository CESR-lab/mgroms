module mg_namelist

  use mg_mpi
  use mg_tictoc

  implicit none

  !- control integer and real precision (global vparameters)
  integer(kind=4), parameter :: rp = 8, ip = 4

  !- Number of point in the halo of each subdomain
  integer(kind=ip) :: nhalo       =   1

  ! smallest dimension ever for the global domain
  integer(kind=ip) :: nsmall      =   8

  integer(kind=ip) :: ns_coarsest = 40
  integer(kind=ip) :: ns_pre      =   2
  integer(kind=ip) :: ns_post     =   2

  character(len=16) :: cmatrix='real'         !- 'simple' or 'real'

  character(len=16) :: mpiexchange='blocking' !- 'blocking' or 'nonblocking'

  logical           :: aggressive = .false.   !- .false. or .true.

  namelist/nhparam/ &
       nhalo      , &
       nsmall     , &
       ns_coarsest, &
       ns_pre     , &
       ns_post    , &
       cmatrix    , &
       mpiexchange, &
       aggressive

 contains

   !--------------------------------------------------------------------
   subroutine read_namelist(filename, verbose)

     character(len=*), optional, intent(in) :: filename
     logical         , optional, intent(in) :: verbose

     character(len=64) :: fn_nml
     logical           :: vb
     integer(kind=ip)  :: lun_nml = 4

     if (present(verbose)) then
        vb = verbose
     else
        vb = .true.
     endif

     if (present(filename)) then
        fn_nml = filename
     else
        fn_nml = 'nh_namelist'
     endif

     open(unit=lun_nml, File=fn_nml, ACTION='READ')

     rewind(unit=lun_nml)
     read(unit=lun_nml, nml=nhparam)

     if (vb) then
        if (myrank == 0) then
           write(*,*)'Namelist non hydrostatic parameters:'
           write(*,*)'  - nhalo      : ', nhalo
           write(*,*)'  - nsmall     : ', nsmall 
           write(*,*)'  - ns_coarsest: ', ns_coarsest
           write(*,*)'  - ns_pre     : ', ns_pre
           write(*,*)'  - ns_post    : ', ns_post
           write(*,*)'  - cmatrix    : ', trim(cmatrix)
           write(*,*)'  - mpiexchange: ', trim(mpiexchange)
           write(*,*)'  - aggressive : ', aggressive
        endif
     endif

   end subroutine read_namelist


 end module mg_namelist
