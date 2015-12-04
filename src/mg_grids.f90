module mg_grids

  use mg_mpi
  use mg_tictoc

  implicit none

  integer(kind=4), parameter:: maxlev=10

  type grid_type
     real(kind=rl),dimension(:,:,:)  ,pointer :: p,b,r
     real(kind=rl),dimension(:,:,:,:),pointer :: cA
     integer(kind=is) :: nx,ny, nz
     integer(kind=is) :: npx, npy, incx, incy
     integer(kind=is) :: nh                 ! number of points in halo
     integer(kind=is) :: gather
     integer(kind=is),dimension(8)::neighb
     real(kind=rl), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
     real(kind=rl), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
     real(kind=rl), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
     real(kind=rl), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE
  end type grid_type

  type(grid_type), dimension(:), pointer :: grid

  integer(kind=is):: nlevs ! index of the coarsest level (1 is the finest)
  integer(kind=is):: nhalo

  !- put it in namelist file !
  logical :: aggressive = .false.

contains

  !----------------------------------------
  subroutine define_grids(npxg, npyg, nxl, nyl, nzl)

    integer(kind=is), intent(in) :: npxg,npyg  ! global CPU topology
    integer(kind=is), intent(in) :: nxl, nyl, nzl ! local dims

    integer(kind=is) :: nhalo        ! number of halo points
    integer(kind=is) :: nh, nd
    integer(kind=is) :: npx, npy
    integer(kind=is) :: incx, incy

    integer(kind=is) :: nx, ny, nz


    integer(kind=is) :: lev

    call  find_grid_levels(npxg, npyg, nzl, nyl, nxl)

    allocate(grid(nlevs))

    nhalo = 2

    grid(1)%nx = nxl 
    grid(1)%ny = nyl
    grid(1)%nz = nzl
    grid(1)%nh = nhalo

    grid(1)%npx = npxg
    grid(1)%npy = npyg

    grid(1)%incx=1
    grid(1)%incy=1

    ! define grid dimensions at each level
    call define_grid_dims()

    ! Allocate memory
    do lev=1,nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       nh = grid(lev)%nh

       if (nz == 1) then
          nd = 3
       else
          nd = 8
       endif

       allocate(grid(lev)%p(    nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%b(    nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%r(    nz,1-nh:ny+nh,1-nh:nx+nh)) ! Need or not ?
       allocate(grid(lev)%cA(nd,nz,1-nh:ny+nh,1-nh:nx+nh))

       allocate(grid(lev)%sendS(nz,nh,nx))
       allocate(grid(lev)%recvS(nz,nh,nx))
       allocate(grid(lev)%sendN(nz,nh,nx))
       allocate(grid(lev)%recvN(nz,nh,nx))

       allocate(grid(lev)%sendE(nz,ny,nh))
       allocate(grid(lev)%recvE(nz,ny,nh))
       allocate(grid(lev)%sendW(nz,ny,nh))
       allocate(grid(lev)%recvW(nz,ny,nh))

       allocate(grid(lev)%sendSW(nz,nh,nh))
       allocate(grid(lev)%sendSE(nz,nh,nh))
       allocate(grid(lev)%sendNW(nz,nh,nh))
       allocate(grid(lev)%sendNE(nz,nh,nh))

       allocate(grid(lev)%recvSW(nz,nh,nh))
       allocate(grid(lev)%recvSE(nz,nh,nh))
       allocate(grid(lev)%recvNW(nz,nh,nh))
       allocate(grid(lev)%recvNE(nz,nh,nh))
    enddo

  end subroutine define_grids

 !----------------------------------------
  subroutine find_grid_levels(npxg, npyg, nz,ny,nx)

    integer(kind=4) , intent(in) :: npxg, npyg
    integer(kind=is), intent(in) :: nx, ny, nz

    integer(kind=is) :: nxg, nyg, nzg
    integer(kind=is) :: npx, npy
    integer(kind=is) :: nsmall


    nxg = npxg * nx
    nyg = npyg * ny
    nzg = nz

    nlevs = 1
    nsmall=8 ! smallest dimension ever for the global domain 

    do 
       ! stop criterion
       if (min(nxg,nyg).lt.nsmall) exit

       nlevs = nlevs+1

       if(nz.eq.1)then 
          ! 2D coarsening
          if ((mod(nxg,2) == 0).or.(mod(nyg,2) == 0)) then
             nxg = nxg/2
             nyg = nyg/2
          else
             write(*,*)'Error:grid dimension not a power of 2!'
             stop -1
          end if
       else
          ! regular 3D coarsening
          if ((mod(nxg,2) == 0).or.(mod(nyg,2) == 0).or.(mod(nzg,2) == 0)) then
             nxg = nxg/2
             nyg = nyg/2
             nzg = nzg/2
          else
             write(*,*)'Error:grid dimension not a power of 2!'
             stop -1
          end if
       endif

    enddo

  end subroutine find_grid_levels

  !----------------------------------------
  subroutine define_grid_dims()

    integer(kind=is) :: nx, ny, nz, nd, nh
    integer(kind=is) :: npx, npy
    integer(kind=is) :: lev, n2d, incx, incy, nsmall
    integer(kind=is) :: coarsen, smooth

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh
    npx = grid(1)%npx
    npy = grid(1)%npy

    incx = 1
    incy = 1
    nsmall=8 ! smallest dimension ever for the global domain 

    do lev = 2, nlevs

       if (aggressive.and.(lev==2)) then
          if (mod(nz,8) == 0) then
             nz = nz/8
          else
             write(*,*)'Error: aggressive coarsening not possible'
             stop -1
          endif
       else
          if (nz.eq.1) then ! 2D coarsening
             nx = nx/2
             ny = ny/2
          else              ! regular 3D coarsening
             nx = nx/2
             ny = ny/2
             nz = nz/2
          endif
       endif

       ! determine if gathering is needed
       !- assumes squarish nxg nyg dimensions !

       grid(lev)%gather=0

       if((nx < nsmall).and.(npx > 1))then
          npx  = npx/2
          nx   = nx*2
          incx = incx*2
          grid(lev)%gather = grid(lev)%gather + 1
       endif

       if((ny < nsmall).and.(npy > 1))then
          npy  = npy/2
          ny   = ny*2
          incy = incy*2
          grid(lev)%gather = grid(lev)%gather + 2
       endif

       grid(lev)%nx   = nx
       grid(lev)%ny   = ny
       grid(lev)%nz   = nz
       grid(lev)%npx  = npx
       grid(lev)%npy  = npy
       grid(lev)%incx = incx
       grid(lev)%incy = incy
       grid(lev)%nh   = nh

    enddo

  end subroutine define_grid_dims

  !----------------------------------------
  subroutine define_neighbours(neighb)
    integer(kind=4), dimension(4), optional, intent(in) :: neighb ! S, E, N, W

    integer(kind=4) :: lev
    integer(kind=is) :: npx, npy
    integer(kind=is) :: incx, incy
    integer(kind=is) :: pi, pj

    ! Neighbours
    do lev=1,nlevs
       npx = grid(1)%npx
       npy = grid(1)%npy
       incx = grid(lev)%incx
       incy = grid(lev)%incy

       pj = myrank/npx
       pi = mod(myrank,npx)

       if (pj >= incy) then ! south
          grid(lev)%neighb(1) = (pj-incy)*npx+pi
       else
          grid(lev)%neighb(1) = MPI_PROC_NULL
       endif

       if (pi < npx-incx) then ! east
          grid(lev)%neighb(2) = pj*npx+pi+incx
       else
          grid(lev)%neighb(2) = MPI_PROC_NULL
       endif

       if (pj < npy-incy) then ! north
          grid(lev)%neighb(3) = (pj+incy)*npx+pi
       else
          grid(lev)%neighb(3) = MPI_PROC_NULL
       endif

       if (pi >= incx) then ! west
          grid(lev)%neighb(4) = pj*npx+pi-incx
       else
          grid(lev)%neighb(4) = MPI_PROC_NULL
       endif

       if ((pj >= incy).and.(pi >= incx)) then ! south west
          grid(lev)%neighb(5) = (pj-incy)*npx+ pi-incx
       else
          grid(lev)%neighb(5) = MPI_PROC_NULL
       endif

       if ((pj >= incy).and.(pi < npx-incx)) then ! south east
          grid(lev)%neighb(6) = (pj-incy)*npx+ pi+incx
       else
          grid(lev)%neighb(6) = MPI_PROC_NULL
       endif

       if ((pj < npy-incy).and.(pi < npx-incx)) then ! north east
          grid(lev)%neighb(7) = (pj+incy)*npx + pi+incx
       else
          grid(lev)%neighb(7) = MPI_PROC_NULL
       endif
       if ((pj < npy-incy).and.(pi >= incx)) then ! north west
          grid(lev)%neighb(8) = (pj+incy)*npx + pi-incx
       else
          grid(lev)%neighb(8) = MPI_PROC_NULL
       endif
    enddo

    if (present(neighb)) then
       ! Test for level 1 the coherency with the ocean model
       If  ((grid(1)%neighb(1) /= neighb(1)).or. &
            (grid(1)%neighb(2) /= neighb(2)).or. &
            (grid(1)%neighb(3) /= neighb(3)).or. &
            (grid(1)%neighb(4) /= neighb(4))) then
          write(*,*)'Error: neighbour definition problem !'
          stop -1
       endif
    endif

  end subroutine define_neighbours

  !----------------------------------------
  subroutine global_max(maxloc,maxglo)
    ! return the global max: maxglo
    ! using the local max on each subdomain
    real(kind=rl),intent(in) :: maxloc
    real(kind=rl),intent(out) :: maxglo

    integer(kind=is) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    ! this is not optimal, but not wrong
    call MPI_ALLREDUCE(maxloc,maxglo,1,MPI_DOUBLE_PRECISION,MPI_max,MPI_COMM_WORLD,ierr)   

  end subroutine global_max

  !----------------------------------------
  subroutine global_sum(lev,sumloc,sumglo)
    ! return the global sum: sumglo
    ! using the local sum on each subdomain
    integer(kind=is),intent(in) :: lev
    real(kind=rl),intent(in) :: sumloc
    real(kind=rl),intent(out) :: sumglo

    integer(kind=is) :: ierr

    ! note: the global comm using MPI_COMM_WORLD is over-kill for levels 
    ! where subdomains are gathered
    call MPI_ALLREDUCE(sumloc,sumglo,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)   
    ! therefore we need to rescale the global sum
    sumglo = sumglo * (grid(lev)%npx*grid(lev)%npy)/(grid(1)%npx*grid(1)%npy)

  end subroutine global_sum


end module mg_grids
