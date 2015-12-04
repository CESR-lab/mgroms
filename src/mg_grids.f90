module mg_grids

  use mg_mpi
  use mg_tictoc

  implicit none

  integer(kind=4), parameter:: maxlev=10

  type grid_type
     ! GR: why are those arrays declared as 'pointer'???
     ! why not regular arrays?
     real(kind=rl),dimension(:,:,:)  ,pointer :: p,b,r,dummy3
     real(kind=rl),dimension(:,:,:,:),pointer :: cA
     !GR real(kind=rl),dimension(:,:,:,:,:),pointer :: gatherbuffer ! a 5D(!) buffer for the gathering
     integer(kind=is) :: nx,ny, nz
     integer(kind=is) :: npx, npy, incx, incy
     integer(kind=is) :: nh                 ! number of points in halo
     integer(kind=is) :: gather
!!$=======
!!$     integer(kind=is) :: Ng, ngx, ngy
!!$     integer:: localcomm ! should be integer (output of MPI_SPLIT)
!!$     integer(kind=is) :: coarsening_method, smoothing_method, gather
!!$     integer(kind=is) :: color,family,key
!!$>>>>>>> 5d76062d572541a52c0574dba607c2e2d63cb883
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
!!$=======
!!$    ! for the gathering
!!$    integer(kind=is) :: ngx, ngy
!!$    integer::     N, ff, family, prevfamily, nextfamily, color, key, localcomm, ierr
!!$    
!!$    ! 1rst loop about the grid dimensions at deifferent levels
!!$    ! at the end of that loop we have the number of levels 
!!$>>>>>>> 5d76062d572541a52c0574dba607c2e2d63cb883

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
!!$=======
!!$    ! Watch out, I continue to use the global indexing
!!$    ! to locate each core
!!$    ! a core that has coordinates (2,3) on the finest decomposition
!!$    ! will remain at this location (2,3) after gathering
!!$    npx = grid(1)%npx ! grid(1) is not a bug!
!!$    npy = grid(1)%npy
!!$    
!!$    pj = myrank/npx
!!$    pi = mod(myrank,npx)
!!$>>>>>>> 5d76062d572541a52c0574dba607c2e2d63cb883

    ! Neighbours
    do lev=1,nlevs       
	   npx = grid(1)%npx
       npy = grid(1)%npy
       ! incx is the distance to my neighbours in x (1, 2, 4, ...)
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
!GR =======
!GR!    return ! not yet ready to go through
!GR    
!GR    ! prepare the informations for the gathering 
 !GR   do lev=1,nlevs-1
!GR       if(grid(lev)%gather.eq.1)then
!GR          
!GR          nx = grid(lev)%nx
!GR          ny = grid(lev)%ny
 !GR         nz = grid(lev)%nz
!GR          nh = grid(lev)%nh
!GR          incx=grid(lev)%incx / 2
!GR          incy=grid(lev)%incy / 2          
!GR          ngx=grid(lev)%ngx
!GR          ngy=grid(lev)%ngy
 

 !GR         !gather cores by quadruplets (and marginally by pair, for the coarsest grid)
!GR
!GR         ! cores having the same family index share the same subdomain
!GR          family=(pi/incx)*incx*incy + (npx)*incy*(pj/incy)
!GR
 !GR         nextfamily = (pi/(2*incx))*incx*incy*4 + (npx)*2*incy*(pj/(incy*2))
!GR
!GR          ! - assign a color to each core: make a cycling ramp index
!GR          ! through 2 or 4 close families 
!GR          ! - cores having the same color should be a pair or a quadruplet 
!GR          ! - colors are all distinct *within* a family
!GR          color=nextfamily + mod(pi,incx)+mod(pj,incy)*incx
!GR
!GR          
!          prevfamily = (pi/(incx/2))*incx*incy/4 + (npx/2)*(incy/2)*(pj/(incy/2))
!GR          N=incx*npx;
!GR          key = mod(mod(family,N)/(incx*incy),2)+2*mod( (family/N),2)
!GR         
!GR
!GR          grid(lev)%color=color
!GR          grid(lev)%family=nextfamily
!GR          grid(lev)%key=key
!GR
!GR          call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, localcomm, ierr)
!GR          grid(lev)%localcomm = localcomm
!GR
!GR
!GR          ! this dummy 3D array is to store the restriction from lev-1, before the gathering
!GR          ! its size can be deduced from the size after the gathering
!GR          
!GR          nx = nx/ngx ! ngx is 1 or 2 (and generally 2)
!GR          ny = ny/ngy ! ngy is 1 or 2 (and generally 2)
!GR          allocate(grid(lev)%dummy3(nz,1-nh:ny+nh,1-nh:nx+nh)) 
!GR          allocate(grid(lev)%gatherbuffer(nz,1-nh:ny+nh,1-nh:nx+nh,ngx,ngy))
!GR          ! number of elements of dummy3
!GR          grid(lev)%Ng=(nx+2*nh)*(ny+2*nh)*nz
!GR
!          if(myrank.eq.0)then
!             write(*,*)"incx=",incx,"Ng=",grid(lev)%Ng,"size(dummy3)=",&
!   size(grid(lev)%dummy3),"size(gatherbuffer)=",size(grid(lev)%gatherbuffer)
!      endif
!GR
!GR      endif
!GR   enddo
!GR
!GR  end subroutine define_grids

  end subroutine define_neighbours

end module mg_grids
