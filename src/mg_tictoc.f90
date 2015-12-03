! intrinsec fortran funciton
! cpu_time(time)  time in second

module mg_tictoc

  implicit none

  integer(kind=4), parameter :: st=4, lg=8

  integer(kind=st), parameter :: levmax=32, submax=32

  real(kind = lg), dimension(levmax,submax) :: ntic
  real(kind = lg), dimension(levmax,submax) :: ntoc
  real(kind = lg), dimension(levmax,submax) :: time
  character(len=32),dimension(submax)    :: subname
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
       !- add its name to "subname"
       !- cpu_time(tic)
       if (flag) then
          nbsub = nbsub + 1
          subname(nbsub)=TRIM(string)
          call cpu_time(ntic(lev,nbsub))
          time(lev,nbsub) = 0._8
       endif

    else
       !- First subroutine to time
       !- add its name to "subname"
       !- cpu_time(tic)
       nbsub = 1
       subname(nbsub)=TRIM(string)
       call cpu_time(ntic(lev,nbsub))
       time(lev,nbsub) = 0._8
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

    integer(kind=st) :: lev

    integer(kind=st) :: ii

    integer(kind=st) ::lun

    if (present(myrank)) then
       lun = myrank + 10
    else
       lun = 10
    endif

    write(lun,'(A)',ADVANCE="no")'   '
    do ii=1, nbsub
       write(lun,'(A16)',ADVANCE="no") TRIM(subname(ii))
    end do

    do lev=1, nblev
       write(lun,'(I2,x,10f16.2)')lev, time(lev,1:nbsub)
    enddo

  end subroutine print_tictoc

end module mg_tictoc
