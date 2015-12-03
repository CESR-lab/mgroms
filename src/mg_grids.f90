module mg_grids

  use mg_mpi

  implicit none

  integer(kind=4), parameter:: maxlev=10

  type grid_type
     ! GR: why are those arrays declared as 'pointer'???
     ! why not regular arrays?
     real(kind=rl),dimension(:,:,:)  ,pointer :: p,b,r,dummy3
     real(kind=rl),dimension(:,:,:,:),pointer :: cA
     real(kind=rl),dimension(:,:,:,:,:),pointer :: gatherbuffer ! a 5D(!) buffer for the gathering
     integer(kind=is) :: nx,ny, nz
     integer(kind=is) :: nd ! number of diagonals in cA
     integer(kind=is) :: npx, npy, incx, incy
     integer(kind=is) :: pj_offset, pi_offset
     integer(kind=is) :: nh                 ! number of points in halo
     integer(kind=is) :: Ng, ngx, ngy
     integer:: localcomm ! should be integer (output of MPI_SPLIT)
     integer(kind=is) :: coarsening_method, smoothing_method, gather
     integer(kind=is) :: color,family,key
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
    integer(kind=is) :: nh, nd
    integer(kind=is) :: npx, npy
    integer(kind=is) :: incx, incy

    integer(kind=is) :: pi, pj
    integer(kind=is) :: lev

    ! for the gathering
    integer(kind=is) :: ngx, ngy
    integer::     N, ff, family, prevfamily, nextfamily, color, key, localcomm, ierr
    
    ! 1rst loop about the grid dimensions at deifferent levels
    ! at the end of that loop we have the number of levels 

    grid(1)%nx = nxo
    grid(1)%ny = nyo
    grid(1)%nz = nzo
    grid(1)%nh = nhalo
    grid(1)%nd = 8

    grid(1)%npx = npxg
    grid(1)%npy = npyg

    grid(1)%incx=1
    grid(1)%incy=1

    grid(1)%pj_offset = 0
    grid(1)%pi_offset = 0

    call define_grid_hierarchy()

!!$    nlevs=4
!!$    do lev=2, nlevs
!!$    ! figure out additional grid dimensions
!!$    ! TODO: figure out process numbers for grid levels
!!$       grid(lev)%nx = grid(lev-1)%nx/2
!!$       grid(lev)%ny = grid(lev-1)%ny/2
!!$       grid(lev)%nz = grid(lev-1)%nz/2
!!$       grid(lev)%nh = nhalo
!!$       grid(lev)%npx = npxg
!!$       grid(lev)%npy = npyg
!!$
!!$       grid(lev)%pj_offset = 0
!!$       grid(lev)%pi_offset = 0
!!$    enddo

    do lev=1,nlevs
       ! Allocate memory
       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       nh = grid(lev)%nh
       nd = grid(lev)%nd
       allocate(grid(lev)%p(   nz,1-nh:ny+nh,1-nh:nx+nh))
       grid(lev)%p(:,:,:)=0._8
       allocate(grid(lev)%b(   nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%r(   nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%cA(nd,nz,1-nh:ny+nh,1-nh:nx+nh))

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

    ! Watch out, I continue to use the global indexing
    ! to locate each core
    ! a core that has coordinates (2,3) on the finest decomposition
    ! will remain at this location (2,3) after gathering
    npx = grid(1)%npx ! grid(1) is not a bug!
    npy = grid(1)%npy
    
    pj = myrank/npx
    pi = mod(myrank,npx)

    ! Neighbours
    do lev=1,nlevs       

       ! incx is the distance to my neighbours in x (1, 2, 4, ...)
       incx = grid(lev)%incx
       incy = grid(lev)%incy

       !! rank/grid(1)%npx = floor(myrank/grid(1)%npx) with integers
!       pj = myrank/grid(lev)%npx    + grid(lev)%pj_offset
!       pi = myrank-pj*grid(lev)%npx + grid(lev)%pi_offset


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

!    return ! not yet ready to go through
    
    ! prepare the informations for the gathering 
    do lev=1,nlevs-1
       if(grid(lev)%gather.eq.1)then
          
          nx = grid(lev)%nx
          ny = grid(lev)%ny
          nz = grid(lev)%nz
          nh = grid(lev)%nh
          incx=grid(lev)%incx / 2
          incy=grid(lev)%incy / 2
          ngx=grid(lev)%ngx
          ngy=grid(lev)%ngy
 

          !gather cores by quadruplets (and marginally by pair, for the coarsest grid)

         ! cores having the same family index share the same subdomain
          family=(pi/incx)*incx*incy + (npx)*incy*(pj/incy)

          nextfamily = (pi/(2*incx))*incx*incy*4 + (npx)*2*incy*(pj/(incy*2))

          ! - assign a color to each core: make a cycling ramp index
          ! through 2 or 4 close families 
          ! - cores having the same color should be a pair or a quadruplet 
          ! - colors are all distinct *within* a family
          color=nextfamily + mod(pi,incx)+mod(pj,incy)*incx

          
!          prevfamily = (pi/(incx/2))*incx*incy/4 + (npx/2)*(incy/2)*(pj/(incy/2))
          N=incx*npx;
          key = mod(mod(family,N)/(incx*incy),2)+2*mod( (family/N),2)
         

          grid(lev)%color=color
          grid(lev)%family=nextfamily
          grid(lev)%key=key

          call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, localcomm, ierr)
          grid(lev)%localcomm = localcomm


          ! this dummy 3D array is to store the restriction from lev-1, before the gathering
          ! its size can be deduced from the size after the gathering
          
          nx = nx/ngx ! ngx is 1 or 2 (and generally 2)
          ny = ny/ngy ! ngy is 1 or 2 (and generally 2)
          allocate(grid(lev)%dummy3(nz,1-nh:ny+nh,1-nh:nx+nh)) 
          allocate(grid(lev)%gatherbuffer(nz,1-nh:ny+nh,1-nh:nx+nh,ngx,ngy))
          ! number of elements of dummy3
          grid(lev)%Ng=(nx+2*nh)*(ny+2*nh)*nz

!          if(myrank.eq.0)then
!             write(*,*)"incx=",incx,"Ng=",grid(lev)%Ng,"size(dummy3)=",&
!   size(grid(lev)%dummy3),"size(gatherbuffer)=",size(grid(lev)%gatherbuffer)
!      endif

       endif
    enddo

  end subroutine define_grids

  !----------------------------------------
  subroutine define_grid_hierarchy

    integer(kind=is) :: nx, ny, nz, nd, nh
    integer(kind=is) :: npx, npy
    integer(kind=is) :: ngx, ngy
    integer(kind=is) :: lev, n2d, incx, incy, nsmall
    integer(kind=is) :: coarsen, smooth
    logical:: ok

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    npx = grid(1)%npx
    npy = grid(1)%npy

    ok = .false.
    lev=2
    n2d = 0 ! number of 2D coarsening
    incx = 1
    incy = 1
    nsmall=8 ! smallest dimension ever for the global domain 
    do while (.not.(ok))
       if(lev>1) then! coarsen
          if(nz.eq.1)then 
             ! 2D coarsening
             n2d=n2d+1
             nx=nx/2
             ny=ny/2
             coarsen=2
          else
!!$             if(aspect_ratio.le.0.25)then !pure vertical coarsening
!!$                aspect_ratio=aspect_ratio*8
!!$                nz=nz/8
!!$                coarsen=5
!!$             else 
             ! regular 3D coarsening
             nx=nx/2
             ny=ny/2
             nz=nz/2
             coarsen=1
!!$             endif
          endif
          grid(lev-1)%coarsening_method=coarsen
       endif

       ! determine if gathering is needed
       if(((nx.lt.nsmall).or.(ny.lt.nsmall)).and.(npx*npy.ge.2))then
          if(npx.ge.2)then
             npx=npx/2
             nx=nx*2
             incx=incx*2
             ngx=2
          else
             ngx=1
          endif
          if(npy.ge.2)then
             npy=npy/2
             ny=ny*2
             incy=incy*2
             ngy=2
          else
             ngy=1
          endif
          ! NOTE THIS IMPORTANT POINT
          ! if level=lev has the gather flag activated
          ! it means that the gathering occurs between lev-1 and lev
          ! (not the other way around)
          grid(lev)%gather=1
       else
          grid(lev)%gather=0
       endif

       !now determine what is the appropriate smoothing on this grid
       if(nz.eq.1)then
          smooth=2
          nd=3
       else
!!$          if(aspect_ratio.le.0.25)then
!!$             smooth=1
!!$             nd=8
!!$          else
          smooth=0
          nd=8
!!$          endif
       endif
       
       grid(lev)%nx=nx
       grid(lev)%ny=ny
       grid(lev)%nz=nz
       grid(lev)%npx=npx
       grid(lev)%npy=npy
       grid(lev)%ngx=ngx
       grid(lev)%ngy=ngy
       grid(lev)%incx=incx
       grid(lev)%incy=incy
       grid(lev)%nd=nd ! number of coefficients for cA
       grid(lev)%nh=nh
       grid(lev)%smoothing_method=smooth

       ! stop criterion
       ok = (nz.eq.1).or.((n2d.eq.2).or.(min(nx,ny).lt.nsmall))
!       ok = (nx.le.8)
       ! stop after 2 pure horizontal coarsenings
       lev = lev+1
    enddo
    nlevs = lev-1

  end subroutine define_grid_hierarchy

  !----------------------------------------
  subroutine fill_halo(lev,p)

    integer(kind=is), intent(in):: lev

    real(kind=rl),dimension( &
         grid(lev)%nz,       &
         1-grid(lev)%nh:grid(lev)%ny+grid(lev)%nh, &
         1-grid(lev)%nh:grid(lev)%nx+grid(lev)%nh), intent(inout) :: p
 
!    real(kind=rl),dimension(:,:,:),intent(inout)::p

    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: nh
    integer(kind=is) :: south, east, north, west
    integer(kind=is) :: southwest, southeast, northeast, northwest

    integer(kind=is) :: etag, wtag, ntag, stag
    integer(kind=is) :: swtag, setag, nwtag, netag
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

!    p => grid(lev)%p

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
       sendE = p(:,1:ny,nx-nh+1:nx)  
    endif

    etag = 1
    call MPI_SendRecv(                                  &
         sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,etag, &
         recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 1
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
       sendN = p(:,ny-nh+1:ny,1:nx)  
    endif
    ntag = 1
    call MPI_SendRecv(&
         sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,ntag, &
         recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,ntag, &
         MPI_COMM_WORLD,status1,ierr)
    stag = 1
    call MPI_SendRecv(&
         sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,stag, &
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
       sendSE = p(:,1:nh,nx-nh+1:nx)  
    endif

    setag = 1
    call MPI_SendRecv(&
         sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,setag, &
         recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,setag, &
         MPI_COMM_WORLD,status1,ierr)

    swtag = 1
    call MPI_SendRecv(&
         sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,swtag, &
         recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,swtag, &
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
       sendNW = p(:,ny-nh+1:ny,1:nh)  
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx)  
    endif

    netag = 1
    call MPI_SendRecv( &
         sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,netag, &
         recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,netag, &
         MPI_COMM_WORLD,status1,ierr)

    nwtag = 1
    call MPI_SendRecv(&
         sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,nwtag, &
         recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,nwtag, &
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

    !write(*,*)'rank- p(nz/2,0:1,0)       :', myrank, p(nz/2,0:1,0), p(nz/2,1,1)
    !write(*,*)'rank- p(nz/2,0:1,nx/2)    :', myrank, p(nz/2,0:1,nx/2)
    !write(*,*)'rank- p(nz/2,0:1,nx+1)    :', myrank, p(nz/2,0:1,nx+1), p(nz/2,1,nx)
    !write(*,*)'rank- p(nz/2,ny:ny+1,0)   :', myrank, p(nz/2,ny:ny+1,0), p(nz/2,ny,1)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx/2):', myrank, p(nz/2,ny:ny+1,nx/2)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx+1):', myrank, p(nz/2,ny:ny+1,nx+1), p(nz/2,ny,nx)

  end subroutine fill_halo

  !----------------------------------------
  subroutine global_max(lev,maxloc,maxglo)
    ! return the global max: maxglo
    ! using the local max on each subdomain
    integer(kind=is),intent(in) :: lev
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
