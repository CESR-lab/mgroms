program main_model

  use nhydro
  use mpi

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

  integer(kind=4) :: k
  integer(kind=4) :: ierr

  ! global domain dimensions
  nxg   = 128
  nyg   = 128
  nzg   = 128

  npxg  = 2
  npyg  = 2

  call mpi_init(ierr)

  call mpi_comm_size(mpi_comm_world, nprocs, ierr)

  if (nprocs /= (npxg*npyg)) then
     write(*,*) "Error: in number of processes !"
     stop -1
  endif

  nx = nxg / npxg
  ny = nyg / npyg
  nz = nzg

  ! fill the neighbours array
  call define_heighbours(npxg, npyg, neighb)

  ! grid definition
  allocate(dx(ny,nx))
  allocate(dy(ny,nx))
  allocate(zr(  nz,ny,nx))
  allocate(zw(nz+1,ny,nx))

  allocate(u(nz,ny,1:nx+1))
  allocate(v(nz,1:ny+1,nx))
  allocate(w(1:nz+1,ny,nx))

  u(:,:,:)    =  0._8
  v(:,:,:)    =  0._8
  w(2:nz,:,:) = -1._8
  w(1,:,:)    =  0._8
  w(nz+1,:,:) =  0._8

  dx(:,:) = 1._8
  dy(:,:) = 1._8

  do k=1, nz
     zr(k,:,:) = k - 0.5_8
     zw(k,:,:) = k - 1._8
  end do
  
  zw(nz+1,:,:) = nz
  
  write(*,*)'Start main model!'
  ! Everything above this point mimics the calling ocean model 
  !-----------------------------------------------------------
  call nhydro_init(npxg, npyg, neighb, dx, dy, zr, zw)

  !p(k,j,i) and rhs(k,j,i)
  call nhydro_solve(u,v,w)

  call nhydro_clean()
  write(*,*)'Start main model!'

end program main_model

!--------------------------------------------------------------
subroutine define_heighbours(npxg, npyg, neighb)

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

end subroutine define_heighbours
