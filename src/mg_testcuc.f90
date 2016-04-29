program mg_testcuc

  use mg_mpi 
  use mg_tictoc
  use mg_grids
  !  use mg_define_rhs
  use mg_define_matrix
  !  use mg_relax
  !  use mg_intergrids
  use mg_cuc
  !  use mg_solvers
  use nhydro

  implicit none

  integer(kind=4):: nxg    ! global x dimension
  integer(kind=4):: nyg    ! global y dimension
  integer(kind=4):: nzg    ! z dimension
  integer(kind=4):: npxg   ! number of processes in x
  integer(kind=4):: npyg   ! number of processes in y
  integer(kind=4):: it     ! iteration loop number
  integer(kind=4):: nit    ! number of iterations
  integer(kind=4) :: nx, ny, nz ! local dimensions

  integer(kind=4):: nsweeps

  integer(kind=4):: lev,ierr, np, nh, rank,k,inc
  real(kind=8)    :: res,conv,res0,bnorm

  real(kind = lg) :: tstart,tend,perf
  real(kind=8) :: rnxg,rnyg,rnzg
  real(kind=8) :: rnpxg,rnpyg

  character(len = 16) :: filen

  call tic(1,'mg_testcuc')

  !---------------!
  !- Ocean model -!
  !---------------!
  inc  = 1
  nxg  = 1024/inc
  nyg  = 1024/inc
  nzg  =   64/inc

  Lx   =  200d3
  Ly   =  200d3
  Htot = 4d3

  ! global variables define in mg_grids 
  hlim    = 250._8
  theta_b =   6._8
  theta_s =   6._8

  npxg  = 2
  npyg  = 2

  nit     = 8
  nsweeps = 100

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
  call print_grids()
  !!call define_rhs(nxg, nyg, npxg)

  call setup_cuc(inc)

  nh = grid(1)%nh

  do lev=1,1!nlevs
     if (myrank.eq.0)write(*,*)'----------------------------------------'

     !--------------------!
     !- P initialisation -!
     !--------------------!
     grid(lev)%p = 0._8

     !--------------------------!
     !- RHS initialisation (b) -!
     !--------------------------!
     ! call setup_rhs(nx,ny,nz,nh,grid(lev)%b)
     ! call setup_smooth_rhs(4)
     call setup_random_patches()

     if (netcdf_output) then
        call write_netcdf(grid(lev)%b,vname='rhs',netcdf_file_name='rhs.nc',rank=myrank)
     endif

     call norm(lev,grid(lev)%b,grid(lev)%b,nx,ny,nz,bnorm)
     if (myrank.eq.0)write(*,*)"bnorm=",bnorm

     call cpu_time(tstart)

     do it=0, nit

        if (it > 0) call Fcycle()

        call compute_residual(lev,res)

        res=sqrt(res/bnorm)

        if (netcdf_output) then
           write(filen,'("r_",i0.2,".nc")') it
           call write_netcdf(grid(lev)%r,vname='r',netcdf_file_name=filen,rank=myrank)
           write(filen,'("p_",i0.2,".nc")') it
           call write_netcdf(grid(lev)%p,vname='p',netcdf_file_name=filen,rank=myrank)
        endif

        if(it>0)then
           conv=res0/res
           !           conv=log(res0/res)/log(10._8)
        else
           conv=0.
        endif
        res0=res
        if (myrank.eq.0)write(*,1000)"lev=",lev," - ite=",it," - res=",res," - conv=",conv
     enddo
     call cpu_time(tend)

  enddo
1000 format(A,I2,A,I5,A,G16.3,A,F8.3)

  if (myrank == 0) then
     rnpxg=real(grid(1)%npx,kind=8)
     rnpyg=real(grid(1)%npy,kind=8)
     rnxg=real(grid(1)%nx,kind=8)*rnpxg
     rnyg=real(grid(1)%ny,kind=8)*rnpyg
     rnzg=real(grid(1)%nz,kind=8)
     ! the rescaled time should be expressed in terms of error reduction,
     ! therefore the ratio rnorm/rnorm0 [the rnorm0 was missing prior Dec 11th]
     perf = (tend-tstart)*(rnpxg*rnpyg)/(-log(res)/log(10._8))/(rnxg*rnyg*rnzg)
     write(*,*)'--- summary ---'
     write(*,'(A,F8.3,A)')"time spent to solve :",tend-tstart," s"
     write(*,'(A,E10.3)')"rescaled performance:",perf
     write(*,*)'---------------'
  end if

  call mpi_finalize(ierr)

  call toc(1,'mg_testcuc')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testcuc
