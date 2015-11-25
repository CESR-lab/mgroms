module mg_grids

  use mg_mpi

  implicit none

  integer(kind=4), parameter:: maxlev=1

  type grid_type
     real(kind=rl),dimension(:,:,:)  ,pointer :: p,b,r
     real(kind=rl),dimension(:,:,:,:),pointer :: cA
     integer(kind=is) :: nx,ny, nz
     integer(kind=is) :: npx, npy
     integer(kind=is) :: pj_offset, pi_offset
     integer(kind=is):: nh                 ! number of points in halo
     integer(kind=is),dimension(8)::neighb
  end type grid_type

  type type_halo
    real(kind=rl), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rl), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rl), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rl), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE
  end type type_halo

  type(type_halo),dimension(maxlev) :: halo

  type(grid_type), dimension(maxlev) :: grid
  integer:: nlevs ! index of the coarsest level (1 is the finest)
  integer:: nhalo

contains

  !----------------------------------------
  subroutine define_grids(nhalo,npxg,npyg)

    integer(kind=is), intent(in):: nhalo         ! number of halo points
    integer(kind=is), intent(in):: npxg,npyg  ! global CPU topology

    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: nh
    integer(kind=is) :: npx, npy

    integer(kind=is) :: pi, pj
    integer(kind=is) :: lev

    ! 1rst loop about the grid dimensions at deifferent levels
    ! at the end of that loop we have the number of levels 

    grid(1)%nx = nxo
    grid(1)%ny = nyo
    grid(1)%nz = nzo
    grid(1)%nh = nhalo

    grid(1)%npx = npxg
    grid(1)%npy = npyg

    grid(1)%pj_offset = 0
    grid(1)%pi_offset = 0

    !!do lev=2, maxlev
    ! figure out additional grid dimensions
    ! TODO: figure out process numbers for grid levels
    !!nlevs=lev
    nlevs=1
    !!enddo

    do lev=1,nlevs
       ! Allocate memory
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       nh = grid(lev)%nh
       allocate(grid(lev)%p(   nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%b(   nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%cA(8,nz,1-nh:ny+nh,1-nh:nx+nh))

       allocate(halo(lev)%sendS(nz,nh,nx))
       allocate(halo(lev)%recvS(nz,nh,nx))
       allocate(halo(lev)%sendN(nz,nh,nx))
       allocate(halo(lev)%recvN(nz,nh,nx))

       allocate(halo(lev)%sendE(nz,ny,nh))
       allocate(halo(lev)%recvE(nz,ny,nh))
       allocate(halo(lev)%sendW(nz,ny,nh))
       allocate(halo(lev)%recvW(nz,ny,nh))

       allocate(halo(lev)%sendSW(nz,nh,nh))
       allocate(halo(lev)%sendSE(nz,nh,nh))
       allocate(halo(lev)%sendNW(nz,nh,nh))
       allocate(halo(lev)%sendNE(nz,nh,nh))

       allocate(halo(lev)%recvSW(nz,nh,nh))
       allocate(halo(lev)%recvSE(nz,nh,nh))
       allocate(halo(lev)%recvNW(nz,nh,nh))
       allocate(halo(lev)%recvNE(nz,nh,nh))
    enddo

    ! Neighbours
    do lev=1,nlevs
       npx = grid(lev)%npx
       npy = grid(lev)%npy

       !! rank/grid(1)%npx = floor(myrank/grid(1)%npx) with integers
       pj = myrank/grid(1)%npx    + grid(1)%pj_offset
       pi = myrank-pj*grid(1)%npx + grid(1)%pi_offset

       if (pj > 0) then ! south
          grid(lev)%neighb(1) = (pj-1)*npx+pi
       else
          grid(lev)%neighb(1) = MPI_PROC_NULL
       endif

       if (pi < npx-1) then ! east
          grid(lev)%neighb(2) = pj*npx+pi+1
       else
          grid(lev)%neighb(2) = MPI_PROC_NULL
       endif

       if (pj < npy-1) then ! north
          grid(lev)%neighb(3) = (pj+1)*npx+pi
       else
          grid(lev)%neighb(3) = MPI_PROC_NULL
       endif

       if (pi >0) then ! west
          grid(lev)%neighb(4) = pj*npx+pi-1
       else
          grid(lev)%neighb(4) = MPI_PROC_NULL
       endif

       if ((pj > 0).and.(pi > 0)) then ! south west
          grid(lev)%neighb(5) = (pj-1)*npx+ pi-1
       else
          grid(lev)%neighb(5) = MPI_PROC_NULL
       endif

       if ((pj > 0).and.(pi < npx-1)) then ! south east
          grid(lev)%neighb(6) = (pj-1)*npx+ pi+1
       else
          grid(lev)%neighb(6) = MPI_PROC_NULL
       endif

       if ((pj < npy-1).and.(pi < npx-1)) then ! north east
          grid(lev)%neighb(7) = (pj+1)*npx + pi+1
       else
          grid(lev)%neighb(7) = MPI_PROC_NULL
       endif
       if ((pj < npy-1).and.(pi >0)) then ! north west
          grid(lev)%neighb(8) = (pj+1)*npx + pi-1
       else
          grid(lev)%neighb(8) = MPI_PROC_NULL
       endif
    enddo

  end subroutine define_grids

  !----------------------------------------
  subroutine fill_halo(lev)

    integer(kind=is), intent(in):: lev

    real(kind=rl),dimension(:,:,:), pointer :: p

    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: nh
    integer(kind=is) :: south, east, north, west
    integer(kind=is) :: southwest, southeast, northeast, northwest

    integer(kind=is) :: etag, wtag, ntag, stag
    integer(kind=is) :: ierr,status1 
    !!integer(kind=is) :: ilev

    real(kind=rl), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rl), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rl), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rl), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    p => grid(lev)%p

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    sendS => halo(lev)%sendS
    recvS => halo(lev)%recvS
    sendN => halo(lev)%sendN
    recvN => halo(lev)%recvN

    sendE => halo(lev)%sendE
    recvE => halo(lev)%recvE
    sendW => halo(lev)%sendW
    recvW => halo(lev)%recvW

    sendSW => halo(lev)%sendSW
    sendSE => halo(lev)%sendSE
    sendNW => halo(lev)%sendNW
    sendNE => halo(lev)%sendNE

    recvSW => halo(lev)%recvSW
    recvSE => halo(lev)%recvSE
    recvNW => halo(lev)%recvNW
    recvNE => halo(lev)%recvNE

    !-----------
    !    Fill West-East ghost cells
    !
    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh:nx)  
    endif

    etag = 3
    call MPI_SendRecv(                                  &
         sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,etag, &
         recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 4
    call MPI_SendRecv(                                  &
         sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,wtag, &
         recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,wtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (west.ne.MPI_PROC_NULL) then
       p(:,1:ny,1-nh:0) = recvW
    else !!Homogenous Neumann  
       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
    endif

    if (east.ne.MPI_PROC_NULL) then
       p(:,1:ny,nx+1:nx+nh) = recvE
    else !!Homogenous Neumann  
       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
    end if
    !

    !-----------
    !    Fill North-South ghost cells
    !
    !     Pack:
    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
    endif
    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh:ny,1:nx)  
    endif
    ntag = 1
    call MPI_SendRecv(&
         sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,ntag, &
         recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,ntag, &
         MPI_COMM_WORLD,status1,ierr)
    stag = 2
    call MPI_SendRecv&
         (sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,stag, &
         recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,stag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack:
    if (south.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1:nx)  = recvS
    else !!Homogenous Neumann  
       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
    end if

    if (north.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1:nx)  = recvN
    else!!Homogenous Neumann  
       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
    end if
    !
    !-----------
    !    Fill Corner ghost cells
    !
    ! SW and SE !
    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh:nx)  
    endif

    etag = 3
    call MPI_SendRecv(&
         sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,etag, &
         recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 4
    call MPI_SendRecv(&
         sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,wtag, &
         recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,wtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (southwest.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1-nh:0) = recvSW
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,nx+1:nx+nh) = recvSE
    else !!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    end if
    !

    ! NW and NE !
    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh:ny,1:nh)  
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh:ny,nx-nh:nx)  
    endif

    etag = 3
    call MPI_SendRecv( &
         sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,etag, &
         recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 4
    call MPI_SendRecv(&
         sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,wtag, &
         recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,wtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (northwest.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1-nh:0) = recvNW
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    end if
    !
  end subroutine fill_halo


end module mg_grids
