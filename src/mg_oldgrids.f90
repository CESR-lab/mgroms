module mg_grids

  implicit none

  type grid_type
     real*8,dimension(:,:,:)  ,allocatable :: x,b,r,dummy
     real*8,dimension(:,:,:,:),allocatable :: A
     integer:: nx,ny,nz,np,mp,nd
     integer::coarsening_method,smoothing_method,gather
     integer::ix,iy
     integer,dimension(8)::neighb
  end type grid_type


  type(grid_type),dimension(40) :: grid
  integer:: nlevs ! index of the coarsest level (1 is the finest)
  integer:: nhalo

contains

  !----------------------------------------
  subroutine define_grids(nx,ny,nz,np,mp,aspect_ratio,nsmall)
    ! define the size of each grid, this routine makes decisions on
    !
    ! * how should be done the coarsening: 
    !  - R_z^3 (if hz/hx << 1)
    !  - R_x*R_y*R_z (regular case)
    !  - R_x*R_y (if there is only one level in z)
    !
    ! * which smoother should be used on each grid
    !
    integer:: nx,ny,nz,np,mp,nsmall
    real*8:: aspect_ratio

    logical:: ok
    integer:: lev,n2d,coarsen,smooth,nd

    integer :: ix,iy

    nhalo=2

    ok = .false.
    lev=1
    n2d = 0 ! number of 2D coarsening
    ix = 1
    iy = 1
    do while (.not.(ok))

       if(lev>1) then! coarsen
          if(nz.eq.1)then ! 2D coarsening
             n2d=n2d+1
             nx=nx/2
             ny=ny/2
             coarsen=2
          else
             if(aspect_ratio.le.0.25)then !pure vertical coarsening
                aspect_ratio=aspect_ratio*8
                nz=nz/8
                coarsen=5
             else ! regular 3D coarsening
                nx=nx/2
                ny=ny/2
                nz=nz/2
                coarsen=1
             endif
          endif
          grid(lev-1)%coarsening_method=coarsen
       endif

       ! determine if gathering is needed
       if(((nx.lt.nsmall).or.(ny.lt.nsmall)).and.(np*mp.ge.2))then
          if(np.ge.2)then
             np=np/2
             nx=nx*2
          endif
          if(mp.ge.2)then
             mp=mp/2
             ny=ny*2
          endif
          grid(lev)%gather=1
       else
          grid(lev)%gather=0
       endif

       !now determine what is the appropriate smoothing on this grid
       if(nz.eq.1)then
          smooth=2
          nd=3
       else
          if(aspect_ratio.le.0.25)then
             smooth=1
             nd=8
          else
             smooth=0
             nd=8
          endif
       endif

       grid(lev)%nx=nx
       grid(lev)%ny=ny
       grid(lev)%nz=nz
       grid(lev)%np=np
       grid(lev)%mp=mp
       grid(lev)%nd=nd
       grid(lev)%smoothing_method=smooth

       allocate( grid(lev)%x(1-nh:nx+nh,ny,nz) )
       allocate( grid(lev)%b(nx,ny,nz) )
       allocate( grid(lev)%r(nx,ny,nz) )
       allocate( grid(lev)%A(nd,nx,ny,nz) ) ! matrix is symmetric
       ! one the finest grid it has 15 coefficients=>8 independant
       ! it can have 27 coefficients on coarser grids=>13 independant

       if (grid(lev)%gather.eq.1)then
          ! dummy is the array resulting from the restriction on a
          ! gather level (i.e. just prior gathering it with 4 others
          !  cores)

          ! todo: think hard on the dimensions of this dummy array
          ! this is subtle
          !!NG: 16 nov 2015 comment this part of program
          !!$ allocate( grid(lev)%dummy(.,.,.) )
       endif

       ! stop criterion
       ok = (nz.eq.1).and.((n2d.eq.2).or.(min(nx,ny).lt.nsmall))
       ! stop after 2 pure horizontal coarsenings
       lev = lev+1
    enddo
    nlevs = lev-1

  end subroutine define_grids

  !----------------------------------------
  subroutine define_neighbours(lev,myrank,np,mp,ix,iy)

    integer:: lev,myrank,np,mp,ix,iy

    integer:: i,j,im,ip,jm,jp

    i=mod(myrank,np)
    j=myrank/np
    im=mod(i-ix,np)
    jm=mod(j-iy,mp)
    ip=mod(i+ix,np)
    jp=mod(j+iy,mp)

    grid(lev)%neighb(1)=im+jm*np ! sw
    grid(lev)%neighb(2)=i +jm*np ! s
    grid(lev)%neighb(3)=ip+jm*np ! se

    grid(lev)%neighb(4)=im+j*np  ! w
    grid(lev)%neighb(5)=ip+j*np  ! e

    grid(lev)%neighb(6)=im+jp*np ! nw
    grid(lev)%neighb(7)=i +jp*np ! n
    grid(lev)%neighb(8)=ip+jp*np ! ne

  end subroutine define_neighbours

  !----------------------------------------
  subroutine print_grids
    integer:: lev

    do lev=1,nlevs
       if (grid(lev)%gather.eq.0)then
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%np,' x',grid(lev)%mp," procs"
       else
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%np,' x',grid(lev)%mp," procs / gather"
       endif
    enddo
100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)
  end subroutine print_grids

end module mg_grids
