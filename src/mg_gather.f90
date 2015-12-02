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
      
      buffer => grid(lev)%gatherbuffer

      ! number of cores per direction involved in this gathering (1 or 2)
      ngx = grid(lev)%ngx
      ngy = grid(lev)%ngy

      nx = grid(lev)%nx
      ny = grid(lev)%ny
      nz = grid(lev)%nz
      nh = grid(lev)%nh

      ! numel(x)
!      Ng = (nx*ngx+2*nh)*(ny*ngy+2*nh)*nz
      Ng = grid(lev)%Ng


      call MPI_ALLGATHER( x, Ng, MPI_DOUBLE_PRECISION, buffer, Ng, MPI_DOUBLE_PRECISION, grid(lev)%localcomm)

      do m=1,ngy
         ii = 1-nh+(l-1)*nx
         do i=1-nh,nx+nh
            do l=1,ngx
               jj = 1-nh+(k-1)*ny
               do j=1-nh,ny+nh
                  do k=1,nz
                     y(k,jj,ii) = buffer(k,j,i,l,m)
                     jj=jj+1
                  enddo
               enddo
            enddo
            ii=ii+1
         enddo
      enddo

    end subroutine gather

  end module mg_gather
