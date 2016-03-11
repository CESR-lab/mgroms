module mg_gather

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids

  implicit none

  contains

    !----------------------------------------
    subroutine gather(lev,x,y)
      
      ! lev is the level of y, after the gathering

      ! the grid information for x is not available because x does not belong to a grid
      ! x is a temporary array (see mg_intergrids.f90)

      integer(kind=ip),intent(in) :: lev
      real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
      real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y

      integer(kind=ip):: nx,ny,nz,nh
      integer(kind=ip):: ngx,ngy,Ng
      integer(kind=ip):: i,j,k,l,m,ii,jj
      integer(kind=ip):: i0,i1
      integer(kind=ip):: ierr
      real(kind=rp),dimension(:,:,:,:,:),pointer :: buffer

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
!       if(myrank==0)write(*,*)'gather lev, Ng=',lev,Ng,ngx,ngy,nx,ny
!       if(myrank==0)write(*,*)buffer(1,1,1),buffer(1,1,1),buffer(1,1,1),buffer(1,1,1,1,1)

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
         ! copy only the inner points of x into y because 
         ! the halo of x is corrupted
         ! Indeed, x comes from the coarsening
         ! after which we didn't update the halo
         ! because fill_halo is not available for the intermediate grid...
         !
         do l=0,ngx-1
            !
            if(l==0)then
               i0=1-nh
            else
               i0=1
            endif
            if(l==ngx-1)then
               i1=nx+nh
            else
               i1=nx
            endif         
            !
            ii = 1+l*nx
            do i=1,nx
               jj = 1+m*ny
               do j=1,ny
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

      integer(kind=ip),intent(in) :: lev
      real(kind=rp),dimension(:,:,:),pointer,intent(in) :: x
!    real(kind=rp),dimension( &
!         grid(lev)%nz,       &
!         1-grid(lev)%nh:grid(lev)%ny+grid(lev)%nh, &
!         1-grid(lev)%nh:grid(lev)%nx+grid(lev)%nh), intent(in) :: x

      real(kind=rp),dimension(:,:,:),pointer,intent(out) :: y
!    real(kind=rp),dimension( &
!         grid(lev)%nz,       &
!         1-grid(lev)%nh:grid(lev)%ny/grid(lev)%ngy+grid(lev)%nh, &
!         1-grid(lev)%nh:grid(lev)%nx/grid(lev)%ngx+grid(lev)%nh), intent(out) :: y


      integer(kind=ip):: nx,ny,nz,nh
      integer(kind=ip):: ngx,ngy
      integer(kind=ip):: i,j,k,l,m,ii,jj,key

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

!      write(*,'(I3,I4,I4,I4,I3,I3,I3,I3)')myrank,nx,ny,nz,l,m,key,grid(lev)%color!size(x),size(y),l,m
 !     l=0
 !     m=0
!      call MPI_Barrier( MPI_COMM_WORLD,ierr)

!      y = 0.
!      z = 0.
!      return
!      return
      ii = 1-nh+l*nx
      do i=1-nh,nx+nh
         jj = 1-nh+m*ny
         do j=1-nh,ny+nh
            do k=1,nz
               y(k,j,i) = x(k,jj,ii)
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
