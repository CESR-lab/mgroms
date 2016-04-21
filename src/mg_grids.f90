module mg_grids

  use mg_mpi
  use mg_tictoc
  use mg_namelist

  implicit none

  type grid_type
     real(kind=rp),dimension(:,:,:)  ,pointer :: p,b,r,dummy3
     integer(kind=1),dimension(:,:),pointer :: rmask
     real(kind=rp),dimension(:,:,:,:),pointer :: cA
     real(kind=rp),dimension(:,:,:),pointer :: zr
     real(kind=rp),dimension(:,:,:),pointer :: zw
     integer(kind=ip) :: nx,ny, nz
     integer(kind=ip) :: npx, npy, incx, incy
     integer(kind=ip) :: nh                    ! number of points in halo
     integer(kind=ip) :: gather
     integer(kind=ip) :: Ng, ngx, ngy
     integer(kind=4) :: localcomm ! should be integer (output of MPI_SPLIT)
     integer(kind=ip) :: coarsening_method, smoothing_method
     integer(kind=ip) :: color,family,key
     integer(kind=ip),dimension(8)::neighb
     real(kind=rp), dimension(:,:,:,:,:),pointer :: gatherbuffer
     real(kind=rp), dimension(:,:,:,:,:),pointer :: gatherbufferp
     real(kind=rp), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
     real(kind=rp), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
     real(kind=rp), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
     real(kind=rp), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE
  end type grid_type

  type(grid_type), dimension(:), pointer :: grid

  integer(kind=ip):: nlevs ! index of the coarsest level (1 is the finest)

  integer(kind=ip),dimension(27,3) :: loc

contains

  !----------------------------------------
  subroutine define_grids(npxg, npyg, nxl, nyl, nzl)

    integer(kind=ip), intent(in) :: npxg,npyg  ! global CPU topology
    integer(kind=ip), intent(in) :: nxl, nyl, nzl ! local dims

    integer(kind=ip) :: nh, nd

    integer(kind=ip) :: nx, ny, nz

    integer(kind=ip) :: lev

    if (myrank==0) write(*,*)'- define grids:'
    if (myrank==0) write(*,*)'  - define grid levels'
    call find_grid_levels(npxg, npyg, nxl, nyl, nzl)

    allocate(grid(nlevs))

    grid(1)%nx = nxl 
    grid(1)%ny = nyl
    grid(1)%nz = nzl
    grid(1)%nh = nhalo

    grid(1)%npx = npxg
    grid(1)%npy = npyg

    grid(1)%incx=1
    grid(1)%incy=1

    if (myrank==0) write(*,*)'  - define grid dims'
    ! define grid dimensions at each level
    call define_grid_dims()

    ! Allocate memory
    do lev=1,nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz
       nh = grid(lev)%nh

       if ((trim(interp_type)=='nearest') .and. (trim(restrict_type)=='avg')) then

          if (nz == 1) then
             nd = 5
          else
             nd = 8
          endif

       elseif (( trim(interp_type)=='linear') .and. (trim(restrict_type)=='avg')) then

          if (nz == 1) then
             nd = 9
          else
             if (lev == 1) then
                nd = 8
             else
                nd =27
             endif
          endif

       elseif (( trim(interp_type)=='nearest') .and. (trim(restrict_type)=='linear')) then
          ! todo 
       endif

       allocate(grid(lev)%p(    nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%b(    nz,1-nh:ny+nh,1-nh:nx+nh))
       allocate(grid(lev)%r(    nz,1-nh:ny+nh,1-nh:nx+nh)) ! Need or not ?
       allocate(grid(lev)%cA(nd,nz,1-nh:ny+nh,1-nh:nx+nh))

       allocate(grid(lev)%zr(nz,-1:ny+2,-1:nx+2))
       allocate(grid(lev)%zw(nz+1,-1:ny+2,-1:nx+2))

       allocate(grid(lev)%rmask(1-nh:ny+nh,1-nh:nx+nh))

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

    grid(1)%p(:,:,:) = 0._8

    if (myrank==0) write(*,*)'  - define gather informations'
    call define_gather_informations()

  end subroutine define_grids

 !----------------------------------------
  subroutine find_grid_levels(npxg, npyg, nx,ny,nz)

    integer(kind=ip) , intent(in) :: npxg, npyg
    integer(kind=ip), intent(in) :: nx, ny, nz

    integer(kind=ip) :: nxg, nyg, nzg

    integer(kind=ip) :: ncoarsest,nhoriz,nzmin, nl1,nl2
    
    nxg = npxg * nx
    nyg = npyg * ny
    nzg = nz

    ! smallest horizontal dimension of the coarsest grid
    ncoarsest = 4 ! TODO: put it into the namelist

    nzmin = 2

    ! smallest horizontal dimension of the finest grid
    nhoriz = min(nxg,nyg)

    ! we have 
    ! nhoriz = ncoarsest * 2^(nlevs-1)
    ! thus nlevs = ...
    nl1 = 1+floor( log( nhoriz*1._8 / ncoarsest*1._8) / log(2._8) )

    nl2 = 1+floor( log( nzg*1._8 / nzmin*1._8) / log(2._8) )

    nlevs=min(nl1,nl2)

    return
    ! the code below can be removed
    ! nlevs = 1   

    ! do 
    !    ! stop criterion
    !    if (min(nxg,nyg).le.nsmall*2) exit

    !    nlevs = nlevs+1

    !    if(nz.eq.1)then 
    !       ! 2D coarsening
    !       if ((mod(nxg,2) == 0).or.(mod(nyg,2) == 0)) then
    !          nxg = nxg/2
    !          nyg = nyg/2
    !       else
    !          write(*,*)'Error:grid dimension not a power of 2!'
    !          stop -1
    !       end if
    !    else
    !       ! regular 3D coarsening
    !       if ((mod(nxg,2) == 0).or.(mod(nyg,2) == 0).or.(mod(nzg,2) == 0)) then
    !          nxg = nxg/2
    !          nyg = nyg/2
    !          nzg = nzg/2
    !       else
    !          write(*,*)'Error:grid dimension not a power of 2!'
    !          stop -1
    !       end if
    !    endif

    ! enddo

  end subroutine find_grid_levels

  !----------------------------------------
  subroutine define_grid_dims()

    integer(kind=ip) :: nx, ny, nz, nh
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: lev, incx, incy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh
    npx = grid(1)%npx
    npy = grid(1)%npy

    incx = 1
    incy = 1

    lev=1
    grid(lev)%gather=0
    grid(lev)%ngx = 1
    grid(lev)%ngy = 1

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
       grid(lev)%ngx = 1
       grid(lev)%ngy = 1

       if((min(nx,ny)<nsmall).and.(npx*npy>1))then
          grid(lev)%gather = 1
          if (npx > 1)then
             npx  = npx/2
             nx   = nx*2             
             grid(lev)%ngx = 2
          endif
          if (npy > 1)then
             npy  = npy/2
             ny   = ny*2             
             grid(lev)%ngy = 2
          endif
          incx=incx*2
          incy=incy*2

       endif

!!$       if((nx < nsmall).and.(npx > 1))then
!!$          npx  = npx/2
!!$          nx   = nx*2
!!$!          incx = incx*2
!!$          grid(lev)%gather = 1
!!$          grid(lev)%ngx = 2
!!$       endif
!!$
!!$       if((ny < nsmall).and.(npy > 1))then
!!$          npy  = npy/2
!!$          ny   = ny*2
!!$!          incy = incy*2
!!$          grid(lev)%gather = 1
!!$          grid(lev)%ngy = 2
!!$       endif
!!$
!!$       if(grid(lev)%gather == 1)then
!!$          incx=incx*2
!!$          incy=incy*2
!!$       endif

       grid(lev)%nx   = nx
       grid(lev)%ny   = ny
       grid(lev)%nz   = nz
       grid(lev)%npx  = npx
       grid(lev)%npy  = npy
       grid(lev)%incx = incx
       grid(lev)%incy = incy
       grid(lev)%nh   = nh


!       if(myrank==0)write(*,*)'11/12/15: lev,npx,npy=',lev,npx,npy
    enddo

  end subroutine define_grid_dims

  !----------------------------------------
  subroutine define_neighbours(neighb)
    integer(kind=ip), dimension(4), optional, intent(in) :: neighb ! S, E, N, W

    integer(kind=ip) :: lev
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: incx, incy
    integer(kind=ip) :: pi, pj
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


    if (myrank==0) write(*,*)'  - define neighbours'

    npx = grid(1)%npx
    npy = grid(1)%npy

    pj = myrank/npx
    pi = mod(myrank,npx)

    ! Neighbours
    do lev=1,nlevs       
       ! incx is the distance to my neighbours in x (1, 2, 4, ...)
       incx = grid(lev)%incx
       incy = grid(lev)%incy

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
    !GR    do lev=1,nlevs-1
    !GR       if(grid(lev)%gather.eq.1)then
    !GR          
    !GR          nx = grid(lev)%nx
    !GR          ny = grid(lev)%ny
    !GR          nz = grid(lev)%nz
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
    !GR          nextfamily = (pi/(2*incx))*incx*incy*4 + (npx)*2*incy*(pj/(incy*2))
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

  !---------------------------------------------------------------------
  subroutine define_gather_informations()

    integer(kind=ip) :: nx, ny, nz, nd
    integer(kind=ip) :: nh
    integer(kind=ip) :: npx, npy
    integer(kind=ip) :: incx, incy

    integer(kind=ip) :: pi, pj 
    integer(kind=ip) :: lev

    ! for the gathering
    integer(kind=ip) :: ngx, ngy
    integer(kind=ip) :: N, family, nextfamily, color, key, localcomm, ierr

    ! Watch out, I continue to use the global indexing
    ! to locate each core
    ! a core that has coordinates (2,3) on the finest decomposition
    ! will remain at this location (2,3) after gathering
    npx = grid(1)%npx ! grid(1) is not a bug!
    npy = grid(1)%npy
    pj = myrank/npx
    pi = mod(myrank,npx)

    do lev=2,nlevs
       if(grid(lev)%gather.eq.1)then

          nx = grid(lev)%nx
          ny = grid(lev)%ny
          nz = grid(lev)%nz
          nh = grid(lev)%nh
          nd = size(grid(lev)%cA,1)
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

          allocate(grid(lev)%gatherbuffer(nz,1-nh:ny+nh,1-nh:nx+nh,0:ngx-1,0:ngy-1))
          allocate(grid(lev)%gatherbufferp(nz+1,1-nh:ny+nh,1-nh:nx+nh,0:ngx-1,0:ngy-1))

          ! number of elements of dummy3
          grid(lev)%Ng=(nx+2*nh)*(ny+2*nh)*nz

       endif
    enddo

  end subroutine define_gather_informations

  !---------------------------------------------------------------------
  subroutine print_grids()
    integer(kind=ip) :: lev,ierr

    if (myrank==0) write(*,*)'- print grid information:'

    call MPI_Barrier( MPI_COMM_WORLD ,ierr)
    if (myrank.eq.0)then
       do lev=1,nlevs
          if (grid(lev)%gather.eq.0)then
             write(*,100)"  lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
          else
             write(*,100)"  lev=",lev,": ", &
                  grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                  " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather"
          endif
       enddo
    endif
100 format (A6,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)

    end subroutine print_grids

  !---------------------------------------------------------------------
  subroutine grids_dealloc()

    integer(kind=ip) :: lev

    do lev=1,nlevs
       deallocate(grid(lev)%p)
       deallocate(grid(lev)%b)
       deallocate(grid(lev)%r)
       deallocate(grid(lev)%cA)

       deallocate(grid(lev)%sendS)
       deallocate(grid(lev)%recvS)
       deallocate(grid(lev)%sendN)
       deallocate(grid(lev)%recvN)

       deallocate(grid(lev)%sendE)
       deallocate(grid(lev)%recvE)
       deallocate(grid(lev)%sendW)
       deallocate(grid(lev)%recvW)

       deallocate(grid(lev)%sendSW)
       deallocate(grid(lev)%sendSE)
       deallocate(grid(lev)%sendNW)
       deallocate(grid(lev)%sendNE)

       deallocate(grid(lev)%recvSW)
       deallocate(grid(lev)%recvSE)
       deallocate(grid(lev)%recvNW)
       deallocate(grid(lev)%recvNE) 
    end do

  end subroutine grids_dealloc

end module mg_grids
