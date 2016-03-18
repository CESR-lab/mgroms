program mg_testintergrids

  use mg_mpi
  use mg_grids
  use mg_intergrids
  use mg_gather
  use mg_netcdf_out

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension
  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: nit    ! number of iterations
  integer(kind=ip):: ngx,ngy
  integer(kind=ip):: nsweeps
  integer(kind=ip):: nx,ny,nz
  integer(kind=ip):: lev, ierr, np, rank
  character(len = 16) :: filen
  real(kind=rp) :: x
  integer(kind=ip):: i,px

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg   = 128*2
  nyg   = 128*2
  nzg   = 128/4

  npxg  = 4
  npyg  = 2

  nit     = 10
  nsweeps = 1
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

  !-------------------!
  !- Enter in nhydro -!
  !-------------------!
  call mg_mpi_init()
  call define_grids(npxg, npyg, nx, ny, nz)
  call define_neighbours()

  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)write(*,*)"---------- grid hierarchy ----------"
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)then
     do lev=1,nlevs
        if (grid(lev)%gather.eq.0)then
           write(*,100)"lev=",lev,": ", &
                grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
        else
           ngx=grid(lev)%ngx
           ngy=grid(lev)%ngy
           write(*,110)"lev=",lev,": ", &
                grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather ",ngx,"x",ngy
        endif
     enddo
  endif
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)


  if (myrank.eq.0)write(*,*)"---------- check fine2coarse ----------"

  grid(1)%b = 3.14_8
  call fill_halo(1,grid(1)%b)

  do lev=1,nlevs-1
     grid(lev)%rmask=1
  enddo
  do lev=1,nlevs-1
     grid(lev)%r=myrank!grid(lev)%b
     call fine2coarse(lev)
     write(filen,'("b_",i1,".nc")') lev+1
     call write_netcdf(grid(lev+1)%b,vname='b',netcdf_file_name=filen,rank=myrank)

     if(myrank==0)write(*,*)"lev =",lev+1," is ok / b(1,1,1)=",grid(lev+1)%b(1,1,1)
!     if((myrank==0).and.(lev==nlevs-2))write(*,*)grid(lev+1)%b
  enddo

  if (myrank.eq.0)write(*,*)"---------- check coarse2fine ----------"
  grid(nlevs)%p=1.2345_8
  grid(5)%p=myrank
  do lev=4,1,-1!nlevs-1,1,-1
!     grid(lev+1)%p=1.2345_8
     grid(lev)%p=0.
     nx = grid(lev+1)%nx
     px= mod(myrank,grid(lev+1)%npx)
!     px= myrank/grid(lev+1)%npx
     do i=0,nx+1
        x = real(i+nx*px)-0.5_8
        grid(lev+1)%p(:,:,i)=x
     enddo

     call coarse2fine(lev)
    write(filen,'("p_",i1,".nc")') lev
    call write_netcdf(grid(lev)%p,vname='p',netcdf_file_name=filen,rank=myrank)

     nx=grid(lev)%nx
     ny=grid(lev)%ny
     if(myrank==0)write(*,*)"lev =",lev," is ok / b(1,1,1)=",grid(lev)%r(1,ny/2,nx/2)
!     if((myrank==0).and.(lev==nlevs-2))write(*,*)grid(lev)%r(1,:,:)
  enddo
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)

100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)
110 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A,I1,A,I1)

  call mpi_finalize(ierr)

end program mg_testintergrids
