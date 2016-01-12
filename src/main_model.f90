program main_model

  use mpi
  use mg_namelist
  use nhydro
  use mg_netcdf_out

  implicit none

  integer(kind=4):: nxg    ! global x dimension
  integer(kind=4):: nyg    ! global y dimension
  integer(kind=4):: nzg    ! z dimension

  integer(kind=4):: npxg   ! number of processes in x
  integer(kind=4):: npyg   ! number of processes in y

  integer(kind=4), dimension(4) :: neighb ! S, E, N, W
  integer(kind=4) :: nx, ny, nz  ! local dimensions

  real(kind=8), dimension(:,:), allocatable   :: dx, dy
  real(kind=8), dimension(:,:,:), allocatable :: zr, zw

  real(kind=8), dimension(:,:,:), allocatable :: u,v,w
  real(kind=8), dimension(:,:), allocatable :: h

  integer(kind=4) :: k
  integer(kind=4) :: np, ierr, rank
  integer(kind=4) :: nh
  integer(kind=4) :: pi, pj
  integer(kind=4) :: i, j 
  real(kind=8)    :: x, y, x0, y0

  ! global domain dimensions
  nxg   = 128*2
  nyg   = 128*2
  nzg   = 64!128/8

  npxg  = 2
  npyg  = 2

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)

  if (np /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !"
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  !- read the NonHydro namelist file if it is present 
  !- else default values and print them (or not).
  call read_nhnamelist(vbrank=rank)

  ! fill the neighbours array
  call mm_define_heighbours(npxg, npyg, neighb)

  allocate(u(nz,ny,1:nx+1))
  allocate(v(nz,1:ny+1,nx))
  allocate(w(1:nz+1,ny,nx))

  u(:,:,:)    =  0._8
  v(:,:,:)    =  0._8
  w(2:nz,:,:) = -1._8
  w(1,:,:)    =  0._8
  w(nz+1,:,:) =  0._8

  if (cmatrix == 'simple') then
     ! grid definition
     allocate(dx(ny,nx))
     allocate(dy(ny,nx))
     allocate(zr(  nz,ny,nx))
     allocate(zw(nz+1,ny,nx))

     dx(:,:) = 1._8
     dy(:,:) = 1._8

     do k=1, nz
        zr(k,:,:) = k - 0.5_8
        zw(k,:,:) = k - 1._8
     end do

     zw(nz+1,:,:) = nz

  elseif (cmatrix == 'real') then

     ! grid definition
     allocate(h(0:ny+1,0:nx+1))
     allocate(dx(0:ny+1,0:nx+1))
     allocate(dy(0:ny+1,0:nx+1))
     allocate(zr(nz,0:ny+1,0:nx+1))
     allocate(zw(nz+1,0:ny+1,0:nx+1))

     dx(:,:) = 1._8
     dy(:,:) = 1._8

     nh = nhalo

     pj = rank/npxg   
     pi = rank-pj*npxg

     x0 = real(nx,kind=rp)
     y0 = real(ny,kind=rp)

!!!  I need to know my global index range
     do i = 0,nx+1
        do j = 0,ny+1
           x = ( real(i+(pi*nx),kind=rp)- 0.5_rp) * dx(i,j)
           y = ( real(j+(pj*ny),kind=rp)- 0.5_rp) * dy(i,j)

           h(j,i) = 4000._rp - 2000._rp * &
                exp(-0.001_rp*((x-x0)**2._rp+(y-y0)**2._rp))
        enddo
     enddo

     call write_netcdf(h,vname='h',netcdf_file_name='h.nc',rank=rank)

     do i = 0,nx+1
        do j = 0,ny+1
           do k = 1,nz
              zr(k,j,i) = (k-0.5_rp)*h(j,i)/nz - h(j,i)
              zw(k,j,i) = (k-1.0_rp)*h(j,i)/nz - h(j,i)
           enddo
           zw(nz+1,j,i) = 0.0_rp
        enddo
     enddo

  else
     stop

  endif

  if (rank.eq.0)  write(*,*)'Start main model!'
  ! Everything above this point mimics the calling ocean model 
  !-----------------------------------------------------------

  call nhydro_init(nx, ny, nz, npxg, npyg, neighb, dx, dy, zr, zw)

  ! define rhs
  grid(1)%b(1:nz,1:ny,1:nx) = 0._8
  grid(1)%b=0.
  call random_number(grid(1)%p(1:nz,1:ny,1:nx))
  call fill_halo(1,grid(1)%p)
  call relax(1,4)
  grid(1)%b=grid(1)%p

  call nhydro_solve(u,v,w)

  if(myrank == 0) call print_tictoc(myrank)

  !  call nhydro_clean()

  call mpi_finalize(ierr)

contains

  !--------------------------------------------------------------
  subroutine mm_define_heighbours(npxg, npyg, neighb)

    use mpi
    implicit none

    integer(kind=4), intent(in) :: npxg, npyg
    integer(kind=4),dimension(:), intent(out)::neighb

    integer(kind=4):: myrank
    integer(kind=4) :: ierr
    integer(kind=4) :: pi, pj

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    ! Neighbour

    pj = myrank/npxg
    pi = mod(myrank,npxg)

    if (pj > 0) then ! south
       neighb(1) = (pj-1)*npxg+pi
    else
       neighb(1) = MPI_PROC_NULL
    endif

    if (pi < npxg-1) then ! east
       neighb(2) = pj*npxg+pi+1
    else
       neighb(2) = MPI_PROC_NULL
    endif

    if (pj < npyg-1) then ! north
       neighb(3) = (pj+1)*npxg+pi
    else
       neighb(3) = MPI_PROC_NULL
    endif

    if (pi > 0) then ! west
       neighb(4) = pj*npxg+pi-1
    else
       neighb(4) = MPI_PROC_NULL
    endif

  end subroutine mm_define_heighbours

end program main_model

