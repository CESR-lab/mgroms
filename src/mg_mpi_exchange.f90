module mg_mpi_exchange

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids

  implicit none

  interface fill_halo
     module procedure   &
          fill_halo_2D, &
          fill_halo_3D, &
          fill_halo_4D
  end interface fill_halo

  type halo
     real(kind=rp), dimension(:,:,:,:), pointer :: sendN,recvN,sendS,recvS
     real(kind=rp), dimension(:,:,:,:), pointer :: sendE,recvE,sendW,recvW
     real(kind=rp), dimension(:,:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
     real(kind=rp), dimension(:,:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE
  end type halo

  type(halo), dimension(40) :: halo4D

  logical :: flag4D = .true.

contains

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_2D(lev,a2D,nhalo)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:), pointer, intent(inout)::a2D
    integer(kind=ip), optional, intent(in):: nhalo

    integer(kind=ip) :: nx, ny
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: ih
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    call tic(lev,'fill_halo_2D_nb')

    nx = grid(lev)%nx
    ny = grid(lev)%ny

    if (present(nhalo)) then
       nh = nhalo
    else
       nh = 1
    endif

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)

    if (nh==1) then
       sendS => grid(lev)%sendS2D1
       recvS => grid(lev)%recvS2D1
       sendN => grid(lev)%sendN2D1
       recvN => grid(lev)%recvN2D1

       sendE => grid(lev)%sendE2D1
       recvE => grid(lev)%recvE2D1
       sendW => grid(lev)%sendW2D1
       recvW => grid(lev)%recvW2D1

    elseif(nh==2) then
       southwest = grid(lev)%neighb(5)
       southeast = grid(lev)%neighb(6)
       northeast = grid(lev)%neighb(7)
       northwest = grid(lev)%neighb(8)

       sendS => grid(lev)%sendS2D2
       recvS => grid(lev)%recvS2D2
       sendN => grid(lev)%sendN2D2
       recvN => grid(lev)%recvN2D2

       sendE => grid(lev)%sendE2D2
       recvE => grid(lev)%recvE2D2
       sendW => grid(lev)%sendW2D2
       recvW => grid(lev)%recvW2D2

       sendSW => grid(lev)%sendSW2D2
       sendSE => grid(lev)%sendSE2D2
       sendNW => grid(lev)%sendNW2D2
       sendNE => grid(lev)%sendNE2D2

       recvSW => grid(lev)%recvSW2D2
       recvSE => grid(lev)%recvSE2D2
       recvNW => grid(lev)%recvNW2D2
       recvNE => grid(lev)%recvNE2D2
    else
       stop
    endif

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    if (nh==2) then
       swnetag = 104
       senwtag = 105
       neswtag = 106
       nwsetag = 107
    endif

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvS,nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
       !! a2D(1-nh:0,1:nx) = a2D(nh:1:-1,1:nx)
       do ih = 1, nh
          if (ih == 1) then
             a2D(0,1:nx) = a2D(1,1:nx)
          elseif(ih == 2) then
             a2D(-1,1:nx) = 2._rp * a2D(1,1:nx) - a2D(2,1:nx)
          else
             stop
          endif
       enddo
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvE,ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
       do ih = 1, nh
          if (ih == 1) then
             a2D(1:ny,nx+1) = a2D(1:ny,nx)
          elseif(ih == 2) then
             a2D(1:ny,nx+2) = 2._rp * a2D(1:ny,nx) - a2D(1:ny,nx-1)
          else
             stop
          endif
       enddo
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                              &
            recvN,nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
       !! a2D(ny+1:ny+nh,1:nx) = a2D(ny:ny-nh+1:-1,1:nx)
       do ih = 1, nh
          if (ih == 1) then
             a2D(ny+1,1:nx) = a2D(ny,1:nx)
          elseif(ih == 2) then
             a2D(ny+2,1:nx) = 2._rp * a2D(ny,1:nx) - a2D(ny-1,1:nx)
          else
             stop
          endif
       enddo
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                             &
            recvW,ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else !!Homogenous Neumann
       !! a2D(1:ny,1-nh:0) = a2D(1:ny,nh:1:-1)
       do ih = 1, nh
          if (ih == 1) then
             a2D(1:ny,0) = a2D(1:ny,1)
          elseif(ih == 2) then
             a2D(1:ny,-1) = 2._rp * a2D(1:ny,1) - a2D(1:ny,2)
          else
             stop
          endif
       enddo
    endif

    if (nh==2) then

       flag_sw_s = .false.
       flag_sw_w = .false.
       if (southwest.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                   &
               recvSW,nh*nh,MPI_DOUBLE_PRECISION,southwest, &
               neswtag,MPI_COMM_WORLD,req(5),ierr)
          comm(5)=5
       elseif (south.ne.MPI_PROC_NULL) then
          flag_sw_s = .true.
       elseif (west.ne.MPI_PROC_NULL) then
          flag_sw_w = .true.
       else !!Homogenous Neumann  
          a2D(1-nh:0,1-nh:0) = a2D(nh:1:-1,nh:1:-1)
       endif

       flag_se_s = .false.
       flag_se_e = .false.
       if (southeast.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                   &
               recvSE,nh*nh,MPI_DOUBLE_PRECISION,southeast, &
               nwsetag,MPI_COMM_WORLD,req(6),ierr)
          comm(6)=6
       elseif (south.ne.MPI_PROC_NULL) then
          flag_se_s = .true.
       elseif (east.ne.MPI_PROC_NULL) then
          flag_se_e = .true.
       else !!Homogenous Neumann  
          a2D(1-nh:0,nx+1:nx+nh) = a2D(nh:1:-1,nx:nx-nh+1:-1)
       endif

       flag_ne_n = .false.
       flag_ne_e = .false.
       if (northeast.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                   &
               recvNE,nh*nh,MPI_DOUBLE_PRECISION,northeast, &
               swnetag,MPI_COMM_WORLD,req(7),ierr)
          comm(7)=7
       elseif (north.ne.MPI_PROC_NULL) then
          flag_ne_n = .true.
       elseif (east.ne.MPI_PROC_NULL) then
          flag_ne_e = .true.
       else !!Homogenous Neumann  
          a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny:ny-nh+1:-1,nx:nx-nh+1:-1)
       endif

       flag_nw_n = .false.
       flag_nw_w = .false.
       if (northwest.ne.MPI_PROC_NULL) then
          call MPI_IRecv(                                   &
               recvNW,nh*nh,MPI_DOUBLE_PRECISION,northwest, &
               senwtag,MPI_COMM_WORLD,req(8),ierr)
          comm(8)=8
       elseif (north.ne.MPI_PROC_NULL) then
          flag_nw_n = .true.
       elseif (west.ne.MPI_PROC_NULL) then
          flag_nw_w = .true.
       else !!Homogenous Neumann  
          a2D(ny+1:ny+nh,1-nh:0) = a2D(ny:ny-nh+1:-1,nh:1:-1)
       endif

    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS(:,:) = a2D(1:nh,1:nx)  
       call MPI_ISend(                              &
            sendS,nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE(:,:) = a2D(1:ny,nx-nh+1:nx) 
       call MPI_ISend(                             &
            sendE,ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN(:,:) = a2D(ny-nh+1:ny,1:nx)
       call MPI_ISend(                              &
            sendN,nx*nh,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW(:,:) = a2D(1:ny,1:nh)  
       call MPI_ISend(                             &
            sendW,ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (nh==2) then
       if (southwest.ne.MPI_PROC_NULL) then
          sendSW(:,:) = a2D(1:nh,1:nh)  
          call MPI_ISend(                                   &
               sendSW,nh*nh,MPI_DOUBLE_PRECISION,southwest, &
               swnetag,MPI_COMM_WORLD,req(13),ierr)
          comm(13)=13
       endif

       if (southeast.ne.MPI_PROC_NULL) then
          sendSE(:,:) = a2D(1:nh,nx-nh+1:nx)  
          call MPI_ISend(                                   &
               sendSE,nh*nh,MPI_DOUBLE_PRECISION,southeast, &
               senwtag,MPI_COMM_WORLD,req(14),ierr)
          comm(14)=14
       endif

       if (northeast.ne.MPI_PROC_NULL) then
          sendNE(:,:) = a2D(ny-nh+1:ny,nx-nh+1:nx) 
          call MPI_ISend(                                   &
               sendNE,nh*nh,MPI_DOUBLE_PRECISION,northeast, &
               neswtag,MPI_COMM_WORLD,req(15),ierr)
          comm(15)=15
       endif

       if (northwest.ne.MPI_PROC_NULL) then
          sendNW(:,:) = a2D(ny-nh+1:ny,1:nh)
          call MPI_ISend(                                   &
               sendNW,nh*nh,MPI_DOUBLE_PRECISION,northwest, &
               nwsetag,MPI_COMM_WORLD,req(16),ierr)
          comm(16)=16
       endif
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          a2D(1-nh:0,1:nx)  = recvS(:,:)

       elseif (indx.eq.2) then ! east
          a2D(1:ny,nx+1:nx+nh) = recvE(:,:)

       elseif (indx.eq.3) then ! north
          a2D(ny+1:ny+nh,1:nx)  = recvN (:,:)

       elseif (indx.eq.4) then ! west
          a2D(1:ny,1-nh:0) = recvW(:,:)

       elseif (indx.eq.5) then ! southwest
          a2D(1-nh:0,1-nh:0) = recvSW(:,:)

       elseif (indx.eq.6) then ! southeast
          a2D(1-nh:0,nx+1:nx+nh) = recvSE(:,:)

       elseif (indx.eq.7) then ! northeast
          a2D(ny+1:ny+nh,nx+1:nx+nh) = recvNE(:,:)

       elseif (indx.eq.8) then ! northwest
          a2D(ny+1:ny+nh,1-nh:0) = recvNW(:,:)

       endif

    enddo      !<-- while

    if (nh==2) then
       !- corners on physicical boundarie if a flag is true-!
       if (flag_sw_s) then
          a2D(1-nh:0,1-nh:0) = a2D(1-nh:0,nh:1:-1)
       elseif (flag_sw_w) then
          a2D(1-nh:0,1-nh:0) = a2D(nh:1:-1,1-nh:0)
       endif

       if (flag_se_s) then
          a2D(1-nh:0,nx+1:nx+nh) = a2D(1-nh:0,nx:nx-nh+1:-1)
       elseif (flag_se_e) then
          a2D(1-nh:0,nx+1:nx+nh) = a2D(nh:1:-1,nx+1:nx+nh)
       endif

       if (flag_ne_n) then
          a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny+1:ny+nh,nx:nx-nh+1:-1)
       elseif (flag_ne_e) then
          a2D(ny+1:ny+nh,nx+1:nx+nh) = a2D(ny:ny-nh+1:-1,nx+1:nx+nh)
       endif

       if (flag_nw_n) then
          a2D(ny+1:ny+nh,1-nh:0) = a2D(ny+1:ny+nh,nh:1:-1)
       elseif (flag_nw_w) then
          a2D(ny+1:ny+nh,1-nh:0) = a2D(ny:ny-nh+1:-1,1-nh:0)
       endif
    endif

    call toc(lev,'fill_halo_2D_nb')

  end subroutine fill_halo_2D

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_3D(lev,p)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:,:), pointer, intent(inout)::p

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    logical :: flag_sw_s, flag_sw_w
    logical :: flag_se_s, flag_se_e
    logical :: flag_ne_n, flag_ne_e
    logical :: flag_nw_n, flag_nw_w

    real(kind=rp), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    call tic(lev,'fill_halo_3D_nb')

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = size(p,dim=1)

    nh = 1

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    call verify_array_shape_3D(grid(lev)%sendS,nz,nh,nx)
    sendS => grid(lev)%sendS
    call verify_array_shape_3D(grid(lev)%recvS,nz,nh,nx)
    recvS => grid(lev)%recvS
    call verify_array_shape_3D(grid(lev)%sendN,nz,nh,nx)
    sendN => grid(lev)%sendN
    call verify_array_shape_3D(grid(lev)%recvN,nz,nh,nx)
    recvN => grid(lev)%recvN

    call verify_array_shape_3D(grid(lev)%sendE,nz,ny,nh)
    sendE => grid(lev)%sendE
    call verify_array_shape_3D(grid(lev)%recvE,nz,ny,nh)
    recvE => grid(lev)%recvE
    call verify_array_shape_3D(grid(lev)%sendW,nz,ny,nh)
    sendW => grid(lev)%sendW
    call verify_array_shape_3D(grid(lev)%recvW,nz,ny,nh)
    recvW => grid(lev)%recvW

    call verify_array_shape_3D(grid(lev)%sendSW,nz,nh,nh)
    sendSW => grid(lev)%sendSW
    call verify_array_shape_3D(grid(lev)%sendSE,nz,nh,nh)
    sendSE => grid(lev)%sendSE
    call verify_array_shape_3D(grid(lev)%sendNW,nz,nh,nh)
    sendNW => grid(lev)%sendNW
    call verify_array_shape_3D(grid(lev)%sendNE,nz,nh,nh)
    sendNE => grid(lev)%sendNE

    call verify_array_shape_3D(grid(lev)%recvSW,nz,nh,nh)
    recvSW => grid(lev)%recvSW
    call verify_array_shape_3D(grid(lev)%recvSE,nz,nh,nh)
    recvSE => grid(lev)%recvSE
    call verify_array_shape_3D(grid(lev)%recvNW,nz,nh,nh)
    recvNW => grid(lev)%recvNW
    call verify_array_shape_3D(grid(lev)%recvNE,nz,nh,nh)
    recvNE => grid(lev)%recvNE

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                 &
            recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                 &
            recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                &
            recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else !!Homogenous Neumann
       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
    endif

    flag_sw_s = .false.
    flag_sw_w = .false.
    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                      &
            recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    elseif (south.ne.MPI_PROC_NULL) then
       flag_sw_s = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_sw_w = .true.
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    flag_se_s = .false.
    flag_se_e = .false.
    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                       &
            recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    elseif (south.ne.MPI_PROC_NULL) then
       flag_se_s = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_se_e = .true.
    else !!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    endif

    flag_ne_n = .false.
    flag_ne_e = .false.
    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                      &
            recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    elseif (north.ne.MPI_PROC_NULL) then
       flag_ne_n = .true.
    elseif (east.ne.MPI_PROC_NULL) then
       flag_ne_e = .true.
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    endif

    flag_nw_n = .false.
    flag_nw_w = .false.
    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                      &
            recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    elseif (north.ne.MPI_PROC_NULL) then
       flag_nw_n = .true.
    elseif (west.ne.MPI_PROC_NULL) then
       flag_nw_w = .true.
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
       call MPI_ISend(                                 &
            sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                &
            sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh+1:ny,1:nx)
       call MPI_ISend(                                 &
            sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
       call MPI_ISend(                                &
            sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
       call MPI_ISend(                                      &
            sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh+1:nx)  
       call MPI_ISend(                                      &
            sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                      &
            sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh+1:ny,1:nh)
       call MPI_ISend(                                      &
            sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          p(:,1-nh:0,1:nx)  = recvS

       elseif (indx.eq.2) then ! east
          p(:,1:ny,nx+1:nx+nh) = recvE

       elseif (indx.eq.3) then ! north
          p(:,ny+1:ny+nh,1:nx)  = recvN 

       elseif (indx.eq.4) then ! west
          p(:,1:ny,1-nh:0) = recvW

       elseif (indx.eq.5) then ! southwest
          p(:,1-nh:0,1-nh:0) = recvSW

       elseif (indx.eq.6) then ! southeast
          p(:,1-nh:0,nx+1:nx+nh) = recvSE

       elseif (indx.eq.7) then ! northeast
          p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE

       elseif (indx.eq.8) then ! northwest
          p(:,ny+1:ny+nh,1-nh:0) = recvNW
       endif

    enddo      !<-- while

    !- corners on physicical boundarie if a flag is true-!
    if (flag_sw_s) then
       p(:,1-nh:0,1-nh:0) = p(:,1-nh:0,nh:1:-1)
    elseif (flag_sw_w) then
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,1-nh:0)
    endif

    if (flag_se_s) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,1-nh:0,nx:nx-nh+1:-1)
    elseif (flag_se_e) then
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx+1:nx+nh)
    endif

    if (flag_ne_n) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny+1:ny+nh,nx:nx-nh+1:-1)
    elseif (flag_ne_e) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx+1:nx+nh)
    endif

    if (flag_nw_n) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny+1:ny+nh,nh:1:-1)
    elseif (flag_nw_w) then
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,1-nh:0)
    endif

    call toc(lev,'fill_halo_3D_nb')

  end subroutine fill_halo_3D

  !----------------------------------------
  subroutine verify_array_shape_3D(pt,dim1,dim2, dim3)

    real(kind=rp), dimension(:,:,:), allocatable, pointer , intent(inout):: pt
    integer(kind=ip), intent(in) :: dim1
    integer(kind=ip), intent(in) :: dim2
    integer(kind=ip), intent(in) :: dim3

    if (allocated(pt)) then
       if (all(shape(pt) == [dim1,dim2,dim3])) then
          ! nothing todo
       else
          deallocate(pt)
          allocate(pt(dim1,dim2,dim3))
       endif
    else
       allocate(pt(dim1,dim2,dim3))
    endif

  end subroutine verify_array_shape_3D

  !----------------------------------------
  subroutine fill_halo_4D(lev,cA)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:,:,:), pointer, intent(inout)::cA

    integer(kind=ip) :: nd

    nd=size(cA,1)

    call alloc_halo4D(lev,nd)

    call fill_halo_4D_nb(lev,cA)

  end subroutine fill_halo_4D

  !----------------------------------------------------------------------------
  !- Nonblocking MPI exchanges -!
  !-----------------------------!
  subroutine fill_halo_4D_nb(lev,cA)

    integer(kind=ip), intent(in):: lev
    real(kind=rp), dimension(:,:,:,:), pointer, intent(inout)::cA

    integer(kind=ip) :: nx, ny, nz, nd
    integer(kind=ip) :: nh
    integer(kind=ip) :: south, east, north, west
    integer(kind=ip) :: southwest, southeast, northeast, northwest

    integer(kind=ip) :: sntag, ewtag, nstag, wetag
    integer(kind=ip) :: swnetag, senwtag, nwsetag, neswtag

    integer(kind=ip) :: i, j
    integer(kind=ip) :: icount
    integer(kind=ip) :: indx
    integer(kind=ip),dimension(16) :: req
    integer(kind=ip),dimension(16) :: comm
    integer(kind=ip),dimension(MPI_STATUS_SIZE) :: status
    integer(kind=ip) :: ierr

    real(kind=rp), dimension(:,:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rp), dimension(:,:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rp), dimension(:,:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rp), dimension(:,:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    call tic(lev,'fill_halo_4D_nb')

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = 1

    nd = size(cA(:,:,:,:),dim=1)

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    sendS => halo4D(lev)%sendS
    recvS => halo4D(lev)%recvS
    sendN => halo4D(lev)%sendN
    recvN => halo4D(lev)%recvN

    sendE => halo4D(lev)%sendE
    recvE => halo4D(lev)%recvE
    sendW => halo4D(lev)%sendW
    recvW => halo4D(lev)%recvW

    sendSW => halo4D(lev)%sendSW
    sendSE => halo4D(lev)%sendSE
    sendNW => halo4D(lev)%sendNW
    sendNE => halo4D(lev)%sendNE

    recvSW => halo4D(lev)%recvSW
    recvSE => halo4D(lev)%recvSE
    recvNW => halo4D(lev)%recvNW
    recvNE => halo4D(lev)%recvNE

    comm(:) = 0
    req(:)  = MPI_REQUEST_NULL

    !- Tag coherency is very important between isend and irecv -!
    sntag   = 100
    ewtag   = 101
    nstag   = 102
    wetag   = 103
    swnetag = 104
    senwtag = 105
    neswtag = 106
    nwsetag = 107

    !-----------------------!
    !- Nonblocking RECEIVE -!
    !-----------------------!

    if (south.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                    &
            recvS,nd*nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            nstag,MPI_COMM_WORLD,req(1),ierr)
       comm(1)=1
    else !!Homogenous Neumann  
       cA(:,:,1-nh:0,1:nx) = 0._8 
    endif

    if (east.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvE,nd*nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            wetag,MPI_COMM_WORLD,req(2),ierr)
       comm(2)=2
    else !!Homogenous Neumann
       cA(:,:,1:ny,nx+1:nx+nh) = 0._8
    endif

    if (north.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                 &
            recvN,nd*nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            sntag,MPI_COMM_WORLD,req(3),ierr)
       comm(3)=3
    else !!Homogenous Neumann  
       cA(:,:,ny+1:ny+nh,1:nx) = 0._8
    endif

    if (west.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                   &
            recvW,nd*nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            ewtag,MPI_COMM_WORLD,req(4),ierr)
       comm(4)=4
    else !!Homogenous Neumann
       cA(:,:,1:ny,1-nh:0) = 0._8
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                      &
            recvSW,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            neswtag,MPI_COMM_WORLD,req(5),ierr)
       comm(5)=5
    else !!Homogenous Neumann  
       cA(:,:,1-nh:0,1-nh:0) = 0._8
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                         &
            recvSE,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            nwsetag,MPI_COMM_WORLD,req(6),ierr)
       comm(6)=6
    else !!Homogenous Neumann  
       cA(:,:,1-nh:0,nx+1:nx+nh) = 0._8
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                      &
            recvNE,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            swnetag,MPI_COMM_WORLD,req(7),ierr)
       comm(7)=7
    else !!Homogenous Neumann  
       cA(:,:,ny+1:ny+nh,nx+1:nx+nh) = 0._8
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       call MPI_IRecv(                                         &
            recvNW,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            senwtag,MPI_COMM_WORLD,req(8),ierr)
       comm(8)=8
    else !!Homogenous Neumann  
       cA(:,:,ny+1:ny+nh,1-nh:0) = 0._8
    endif

    !--------------------!
    !- Nonblocking SEND -!
    !--------------------!

    if (south.ne.MPI_PROC_NULL) then
       sendS = cA(:,:,1:nh,1:nx)  
       call MPI_ISend(                                    &
            sendS,nd*nz*nx*nh,MPI_DOUBLE_PRECISION,south, &
            sntag,MPI_COMM_WORLD,req(9),ierr)
       comm(9)=9
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = cA(:,:,1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                &
            sendE,nd*nz*ny*nh,MPI_DOUBLE_PRECISION,east, &
            ewtag,MPI_COMM_WORLD,req(10),ierr)
       comm(10)=10
    endif

    if (north.ne.MPI_PROC_NULL) then
       sendN = cA(:,:,ny-nh+1:ny,1:nx)
       call MPI_ISend(                                    &
            sendN,nd*nz*nx*nh,MPI_DOUBLE_PRECISION,north, &
            nstag,MPI_COMM_WORLD,req(11),ierr)
       comm(11)=11
    endif

    if (west.ne.MPI_PROC_NULL) then
       sendW = cA(:,:,1:ny,1:nh)  
       call MPI_ISend(                                &
            sendW,nd*nz*ny*nh,MPI_DOUBLE_PRECISION,west, &
            wetag,MPI_COMM_WORLD,req(12),ierr)
       comm(12)=12
    endif

    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = cA(:,:,1:nh,1:nh)  
       call MPI_ISend(                                         &
            sendSW,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,southwest, &
            swnetag,MPI_COMM_WORLD,req(13),ierr)
       comm(13)=13
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = cA(:,:,1:nh,nx-nh+1:nx)  
       call MPI_ISend(                                      &
            sendSE,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,southeast, &
            senwtag,MPI_COMM_WORLD,req(14),ierr)
       comm(14)=14
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = cA(:,:,ny-nh+1:ny,nx-nh+1:nx) 
       call MPI_ISend(                                         &
            sendNE,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,northeast, &
            neswtag,MPI_COMM_WORLD,req(15),ierr)
       comm(15)=15
    endif

    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = cA(:,:,ny-nh+1:ny,1:nh)
       call MPI_ISend(                                         &
            sendNW,nd*nz*nh*nh,MPI_DOUBLE_PRECISION,northwest, &
            nwsetag,MPI_COMM_WORLD,req(16),ierr)
       comm(16)=16
    endif

    !- Wait for completion of receive and fill ghost points

    icount=0                       ! Compress arrays "comm" and
    do i=1,16                      ! "req" to disregard directions
       if (comm(i).gt.0) then      ! in which no message was sent
          icount=icount+1          ! or is expected from.  At the
          if (icount.lt.i) then    ! end of this segment icount
             comm(icount)=comm(i)  ! is equal to the actual number
             req(icount)=req(i)    ! of messages sent and received, 
          endif                    ! arrays comm,req(1:icount)
       endif                       ! store directional indices
    enddo

    do while (icount > 0)

       call MPI_Waitany(icount, req, j, status, ierr)

       indx=comm(j)           ! Save directional index for
       icount=icount-1        ! message received and ready to
       do i=j,icount          ! unpack, then erase its "req"
          req(i)=req(i+1)     ! and "comm" and "req" by 
          comm(i)=comm(i+1)   ! by compressing the arrays, so
       enddo                  ! that the same message will 

       ! be unpacked only once.
       if (indx.eq.1) then ! south
          cA(:,:,1-nh:0,1:nx)  = recvS

       elseif (indx.eq.2) then ! east
          cA(:,:,1:ny,nx+1:nx+nh) = recvE

       elseif (indx.eq.3) then ! north
          cA(:,:,ny+1:ny+nh,1:nx)  = recvN 

       elseif (indx.eq.4) then ! west
          cA(:,:,1:ny,1-nh:0) = recvW

       elseif (indx.eq.5) then ! southwest
          cA(:,:,1-nh:0,1-nh:0) = recvSW

       elseif (indx.eq.6) then ! southeast
          cA(:,:,1-nh:0,nx+1:nx+nh) = recvSE

       elseif (indx.eq.7) then ! northeast
          cA(:,:,ny+1:ny+nh,nx+1:nx+nh) = recvNE

       elseif (indx.eq.8) then ! northwest
          cA(:,:,ny+1:ny+nh,1-nh:0) = recvNW
       endif

    enddo      !<-- while  

    call toc(lev,'fill_halo_4D_nb')

  end subroutine fill_halo_4D_nb

  !----------------------------------------
  subroutine alloc_halo4D(lev,nd)

    integer(kind=ip) :: nx, ny, nz
    integer(kind=ip) :: nh, nd
    integer(kind=ip) :: lev


    !    do lev = 1, nlevs

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = 1

    allocate(halo4D(lev)%sendS(nd,nz,nh,nx))
    allocate(halo4D(lev)%recvS(nd,nz,nh,nx))
    allocate(halo4D(lev)%sendN(nd,nz,nh,nx))
    allocate(halo4D(lev)%recvN(nd,nz,nh,nx))

    allocate(halo4D(lev)%sendE(nd,nz,ny,nh))
    allocate(halo4D(lev)%recvE(nd,nz,ny,nh))
    allocate(halo4D(lev)%sendW(nd,nz,ny,nh))
    allocate(halo4D(lev)%recvW(nd,nz,ny,nh))

    allocate(halo4D(lev)%sendSW(nd,nz,nh,nh))
    allocate(halo4D(lev)%sendSE(nd,nz,nh,nh))
    allocate(halo4D(lev)%sendNW(nd,nz,nh,nh))
    allocate(halo4D(lev)%sendNE(nd,nz,nh,nh))

    allocate(halo4D(lev)%recvSW(nd,nz,nh,nh))
    allocate(halo4D(lev)%recvSE(nd,nz,nh,nh))
    allocate(halo4D(lev)%recvNW(nd,nz,nh,nh))
    allocate(halo4D(lev)%recvNE(nd,nz,nh,nh))

    !    enddo

  end subroutine alloc_halo4D

  !----------------------------------------
  subroutine dealloc_halo4D()

    integer(kind=ip) :: lev

    do lev = 1, nlevs
       deallocate(halo4D(lev)%sendS)
       deallocate(halo4D(lev)%recvS)
       deallocate(halo4D(lev)%sendN)
       deallocate(halo4D(lev)%recvN)

       deallocate(halo4D(lev)%sendE)
       deallocate(halo4D(lev)%recvE)
       deallocate(halo4D(lev)%sendW)
       deallocate(halo4D(lev)%recvW)

       deallocate(halo4D(lev)%sendSW)
       deallocate(halo4D(lev)%sendSE)
       deallocate(halo4D(lev)%sendNW)
       deallocate(halo4D(lev)%sendNE)

       deallocate(halo4D(lev)%recvSW)
       deallocate(halo4D(lev)%recvSE)
       deallocate(halo4D(lev)%recvNW)
       deallocate(halo4D(lev)%recvNE)
    enddo

  end subroutine dealloc_halo4D

  !----------------------------------------
  subroutine global_max(maxloc)
    ! return the global max: maxglo
    ! using the local max on each subdomain
    real(kind=rp),intent(inout) :: maxloc

    real(kind=rp)    :: maxglo
    integer(kind=ip) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    ! this is not optimal, but not wrong
    call MPI_ALLREDUCE(maxloc,maxglo,1,MPI_DOUBLE_PRECISION,MPI_max,MPI_COMM_WORLD,ierr)   

    maxloc = maxglo

  end subroutine global_max

  !----------------------------------------
  subroutine global_sum(lev,sumloc,sumglo)
    ! return the global sum: sumglo
    ! using the local sum on each subdomain
    integer(kind=ip),intent(in) :: lev
    real(kind=rp),intent(in) :: sumloc
    real(kind=rp),intent(out) :: sumglo

    integer(kind=ip) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    call MPI_ALLREDUCE(sumloc,sumglo,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)   

    ! therefore we need to rescale the global sum
    sumglo = sumglo * (grid(lev)%npx*grid(lev)%npy)/(grid(1)%npx*grid(1)%npy)

  end subroutine global_sum


end module mg_mpi_exchange
