module mg_namelist

  use mg_mpi

  implicit none

  !- Number of point in the halo of each subdomain
  integer(kind=4) :: nhalo       =   2

  ! smallest dimension ever for the global domain
  integer(kind=4) :: nsmall      =   8

  integer(kind=4) :: ns_coarsest = 100
  integer(kind=4) :: ns_pre      =   2
  integer(kind=4) :: ns_post     =   2

  logical         :: aggressive = .false.

  namelist/nhparam/ &
       ns_coarsest, &
       ns_pre     , &
       ns_post

 contains

   subroutine read_namelist(filename, verbose)

     character(len=*), optional, intent(in) :: filename
     logical         , optional, intent(in) :: verbose

     character(len=64) :: fn_nml
     logical           :: vb
     integer(kind=4)   :: lun_nml = 4

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

     OPEN(unit=lun_nml, File=fn_nml, ACTION='READ')

     REWIND(unit=lun_nml)
     READ(unit=lun_nml, nml=nhparam)

     if (vb) then
        if (myrank == 0) then
           write(*,*)'Namelist non hydrostatic parameters:'
           write(*,*)'  - nhalo      : ', nhalo
           write(*,*)'  - nsmall     : ', nsmall 
           write(*,*)'  - ns_coarsest: ', ns_coarsest
           write(*,*)'  - ns_pre     : ', ns_pre
           write(*,*)'  - ns_post    : ', ns_post
           write(*,*)'  - aggressive : ', aggressive
        endif
     endif

   end subroutine read_namelist


 end module mg_namelist
