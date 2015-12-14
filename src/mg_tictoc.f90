module mg_tictoc
! intrinsec fortran funciton
! cpu_time(time)  time in second

  implicit none

  integer(kind=4), parameter :: st=4, lg=8

  integer(kind=st), parameter :: levmax=32, submax=32

  real(kind = lg)  , dimension(levmax,submax) :: ntic
  real(kind = lg)  , dimension(levmax,submax) :: ntoc
  real(kind = lg)  , dimension(levmax,submax) :: time
  integer(kind=st) , dimension(levmax,submax) :: calls
  character(len=32),dimension(submax)         :: subname
  integer(kind=st) :: nblev = 0, nbsub = 0

contains

 !------------------------------------------------
  subroutine tic(lev, string)
    integer(kind=st), intent(in) :: lev
    character(len=*), intent(in) :: string

    integer(kind=st) :: ns
    logical :: flag 

    if (nbsub > 0) then

       flag = .true.

       !- Search if subroutine is already timed
       !- if yes -> cpu_time(tic)
       do ns=1, nbsub
          if (TRIM(string) == subname(ns)) then
             call cpu_time(ntic(lev,ns))
             flag = .false.
             exit
          endif
       end do

       !- New subroutine to time
       !- Add its name to "subname"
       !- cpu_time(tic)
       if (flag) then
          nbsub = nbsub + 1
          subname(nbsub)=TRIM(string)
          call cpu_time(ntic(lev,nbsub))
          time(lev,nbsub)  = 0._8
          calls(lev,nbsub) = 0
       endif

    else
       !- First subroutine to time
       !- add its name to "subname"
       !- cpu_time(tic)
       nbsub = 1
       subname(nbsub)=TRIM(string)
       call cpu_time(ntic(lev,nbsub))
       time(lev,nbsub)  = 0._8
       calls(lev,nbsub) = 0
    endif

    if (lev > nblev) nblev = lev

  end subroutine tic

 !------------------------------------------------
  subroutine toc(lev, string)
    integer(kind=st), intent(in) :: lev
    character(len=*), intent(in) :: string

    integer(kind=st) :: ns
    logical :: flag 

    if (nbsub > 0) then

       flag = .true.

       do ns=1, nbsub
          if (TRIM(string) == subname(ns)) then
             call cpu_time(ntoc(lev,ns))
             time(lev,ns) = time(lev,ns) + ntoc(lev,ns) - ntic(lev,ns)
             calls(lev,ns) = calls(lev,ns) + 1
             if (lev > nblev) nblev = lev
             flag = .false.
             exit
          endif
       end do

       if (flag) then
          write(*,*)'Error: tictoc: a toc is calling before a tic !'
          write(*,*)'Error: check if a tic exist for:', TRIM(string)
       endif

    else
       write(*,*)'Error: tictoc: a toc is calling before a tic !'
       write(*,*)'Error: check if a tic exist for:', TRIM(string)
    endif

  end subroutine toc

 !------------------------------------------------
  subroutine print_tictoc(myrank)
    integer(kind=st), optional, intent(in)::myrank

    integer(kind=st)  :: lev
    integer(kind=st)  :: ii
    integer(kind=st)  :: lun
    CHARACTER(len=14) :: cmftf, cmfti

    if (present(myrank)) then
       lun = myrank + 10
    else
       lun = 10
    endif

    WRITE(cmftf , 1000) nblev
    WRITE(cmfti , 1001) nblev

1000 FORMAT('(', I3, '(x,E9.3))')
1001 FORMAT('(', I3, '(x,I9))')

    write(lun,'(t22)', ADVANCE="no")
    do lev=1, nblev
       write(lun,'(x,I9)', ADVANCE="no") lev
    enddo

    write(lun,'(x)', ADVANCE="yes")

    do ii=1, nbsub
       write(lun,'(x,A20)' , ADVANCE="no" ) TRIM(subname(ii))
       write(lun,FMT=cmftf , ADVANCE="no" ) time(1:nblev,ii)
       write(lun,'(x)'     , ADVANCE="yes")
       write(lun,'(t22)'   , ADVANCE="no" )
       write(lun,FMT=cmfti , ADVANCE="no" ) calls(1:nblev,ii)
       write(lun,'(x)'     , ADVANCE="yes")
    end do

  end subroutine print_tictoc

end module mg_tictoc
