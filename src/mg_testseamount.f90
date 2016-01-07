program mg_testseamount

  use mpi
  use nhydro

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

  real(kind=8)    :: Lx, Ly, Hc
  real(kind=8)    :: x, x0

  ! global domain dimensions
  nxg   = 128
  nyg   = 128
  nzg   = 128

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

  ! grid definition
  allocate(h(0:ny+1,0:nx+1))
  allocate(dx(0:ny+1,0:nx+1))
  allocate(dy(0:ny+1,0:nx+1))
  allocate(zr(nz,0:ny+1,0:nx+1))
  allocate(zw(nz+1,0:ny+1,0:nx+1))

  Lx = 1.e4
  Ly = 1.e4
  Hc = 4.e3

  dx(:,:) = Lx/nxg
  dy(:,:) = Ly/nyg

  nh = nhalo

  pj = rank/npxg   
  pi = rank-pj*npxg

  x0 = Lx/2._rp

  !!!  I need to know my global index range
  do i = 0,nx+1
     do j = 0,ny+1
        x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(i,j)
        h(j,i) = Hc * (  1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5)**2._rp) )
     enddo
  enddo

  do i = 0,nx+1
     do j = 0,ny+1
        do k = 1,nz
           zr(k,j,i) = (k-0.5_rp)*h(j,i)/nz - h(j,i)
           zw(k,j,i) = (k-1.0_rp)*h(j,i)/nz - h(j,i)
        enddo
        zw(nz+1,j,i) = 0.0_rp
     enddo
  enddo

  if (rank.eq.0)  write(*,*)'Start main model!'

  ! Everything above this point mimics the calling ocean model 
  !-----------------------------------------------------------

  call nhydro_init(nx, ny, nz, npxg, npyg, neighb, dx, dy, zr, zw)

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

end program mg_testseamount

