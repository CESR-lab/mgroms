module mg_gather

  use mg_mpi 

  use mg_grids

  implicit none

  contains

    !----------------------------------------
    subroutine gather(lev,x,y)
      
      ! lev is the level of y, after the gathering

      ! the grid information for x is not available because x does not belong to a grid
      ! x is a temporary array (see mg_intergrids.f90)

      integer(kind=is),intent(in) :: lev
      real(kind=rl),dimension(:,:,:),intent(in) :: x
      real(kind=rl),dimension(:,:,:),intent(out) :: y

      integer(kind=is):: nx,ny,nz,nh,N
      integer(kind=is):: ngx,ngy,Ng
      integer(kind=is):: i,j,k,l,m,ii,jj
      real(kind=rl),dimension(:,:,:,:,:),pointer :: buffer
      integer(kind=is) ::ierr
      
      real(kind=rl):: z
      real(kind=rl),dimension(2,2)::b

      buffer => grid(lev)%gatherbuffer

      ! number of cores per direction involved in this gathering (1 or 2)
      ngx = grid(lev)%ngx
      ngy = grid(lev)%ngy

      nx = grid(lev)%nx
      ny = grid(lev)%ny
      nz = grid(lev)%nz
      nh = grid(lev)%nh

      ! numel(x)
      Ng = grid(lev)%Ng

       call MPI_ALLGATHER( x, Ng, MPI_DOUBLE_PRECISION, buffer, Ng, MPI_DOUBLE_PRECISION, grid(lev)%localcomm,ierr)

      ! I can see two possibilities to copy the 4 buffers into y
      !
      ! 1/ either sweep across each buffer and copy it to y
      !
      ! 2/ sweep across y and pick the value from the proper buffer
      !
      ! what's the fastest???
      !
      nx = nx / ngx
      ny = ny / ngy
      do m=0,ngy-1
         do l=0,ngx-1
            ii = 1-nh+l*nx
            do i=1-nh,nx+nh
               jj = 1-nh+m*ny
               do j=1-nh,ny+nh
                  do k=1,nz
                     y(k,jj,ii) = buffer(k,j,i,l,m)
                  enddo
                  jj=jj+1
               enddo
               ii=ii+1
            enddo
         enddo
      enddo

    end subroutine gather

    !----------------------------------------
    subroutine split(lev,x,y)
      
      ! lev is the level of x, where it has to be split
      ! y is the dummy 3D intermediate array, before interpolation

      integer(kind=is),intent(in) :: lev
      real(kind=rl),dimension(:,:,:),intent(in) :: x
!    real(kind=rl),dimension( &
!         grid(lev)%nz,       &
!         1-grid(lev)%nh:grid(lev)%ny+grid(lev)%nh, &
!         1-grid(lev)%nh:grid(lev)%nx+grid(lev)%nh), intent(in) :: x

      real(kind=rl),dimension(:,:,:),intent(out) :: y
!    real(kind=rl),dimension( &
!         grid(lev)%nz,       &
!         1-grid(lev)%nh:grid(lev)%ny/grid(lev)%ngy+grid(lev)%nh, &
!         1-grid(lev)%nh:grid(lev)%nx/grid(lev)%ngx+grid(lev)%nh), intent(out) :: y


      integer(kind=is):: nx,ny,nz,nh
      integer(kind=is):: ngx,ngy
      integer(kind=is):: i,j,k,l,m,ii,jj,key,ierr
      real(kind=rl):: z

      ! number of cores per direction involved in this gathering (1 or 2)
      ngx = grid(lev)%ngx
      ngy = grid(lev)%ngy

      nx = grid(lev)%nx / ngx
      ny = grid(lev)%ny / ngy
      nz = grid(lev)%nz
      nh = grid(lev)%nh

      key = grid(lev)%key

      l = mod(key,2)
      m = key/2

!      write(*,*)myrank,nx,ny,nz,l,m,nh!size(x),size(y),l,m
 !     l=0
 !     m=0
!      call MPI_Barrier( MPI_COMM_WORLD)

!      y = 0.
!      z = 0.
!      return
!      return
      ii = 1-nh+l*nx
      do i=1-nh,nx+nh
         jj = 1-nh+m*ny
         do j=1-nh,ny+nh
            do k=1,nz
               z=z+1.
               y(1,1,1) = x(k,jj,ii)
            enddo
            jj=jj+1
         enddo
         ii=ii+1
      enddo
!      call MPI_Barrier( MPI_COMM_WORLD)
!      write(*,*)"done",z
!      call MPI_Barrier( MPI_COMM_WORLD,ierr)
      
    end subroutine split

  end module mg_gather
