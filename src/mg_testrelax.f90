program mg_testrelax

  use mg_mpi 
  use mg_tictoc
  !  use mg_grids
  !  use mg_define_rhs
  !  use mg_define_matrix
  !  use mg_relax
  !  use mg_intergrids
  use mg_seamount
  !  use mg_solvers
  use nhydro

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension
  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: it     ! iteration loop number
  integer(kind=ip):: nit    ! number of iterations
  integer(kind=ip) :: nx, ny, nz ! local dimensions

  integer(kind=ip):: nsweeps

  integer(kind=ip):: lev,ierr, np, nh, rank,i,j,k
  real(kind=rp)    :: res,conv,res0,bnorm

  real(kind = lg) :: tstart,tend,perf
  real(kind=rp) :: rnxg,rnyg,rnzg
  real(kind=rp) :: rnpxg,rnpyg

  character(len = 16) :: filen

  call tic(1,'mg_testrelax')

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg   = 512
  nyg   = 512
  nzg   = 32

  Lx =  200d3
  Ly =  200d3
  Htot = 4d3

!!$  nxg   = 128
!!$  nyg   = 128
!!$  nzg   = 128
!!$
!!$  Lx =  20d3
!!$  Ly =  20d3
!!$  Htot = 4d3


  hc      = 100.
  theta_b = 0.
  theta_s = 3.

  npxg  = 2
  npyg  = 2

  nit     = 4
  nsweeps = 2

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


  call setup_realistic_matrix()!2d4,2d4,4d3)
  nh = grid(1)%nh

  do lev=1,1!nlevs
     if (myrank.eq.0)write(*,*)'----------------------------------------'


     grid(lev)%p = 0._8
     !     grid(lev)%b = 0._8
     !     call random_number(grid(lev)%p)
     !     call random_number(grid(lev)%b)
     !     grid(lev)%b(grid(lev)%nz,:,:)= 0.
     do k=1,grid(lev)%nz
        grid(lev)%p(k,:,:)= grid(lev)%p(k,:,:)*grid(lev)%rmask
        !        grid(lev)%b(k,:,:)= (2*grid(lev)%b(k,:,:)-1.)*grid(lev)%rmask
     enddo
     call fill_halo(lev,grid(lev)%p)
     !     call fill_halo(lev,grid(lev)%b)
     call setup_rhs(nx,ny,nz,nh,grid(lev)%b)

     if (netcdf_output) then
        call write_netcdf(grid(lev)%b,vname='rhs',netcdf_file_name='rhs.nc',rank=myrank)
     endif

     call norm(lev,grid(lev)%b,grid(lev)%b,nx,ny,nz,bnorm)


     call cpu_time(tstart)
     do it=0, nit
        !if(it>0) call Vcycle(lev)
        !if(it>0) call Vcycle2(lev,lev+5)
        if(it>0) call Fcycle()
        !if(it>0)call relax(lev,nsweeps)
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
     rnpxg=real(grid(1)%npx,kind=rp)
     rnpyg=real(grid(1)%npy,kind=rp)
     rnxg=real(grid(1)%nx,kind=rp)*rnpxg
     rnyg=real(grid(1)%ny,kind=rp)*rnpyg
     rnzg=real(grid(1)%nz,kind=rp)
     ! the rescaled time should be expressed in terms of error reduction,
     ! therefore the ratio rnorm/rnorm0 [the rnorm0 was missing prior Dec 11th]
     perf = (tend-tstart)*(rnpxg*rnpyg)/(-log(res)/log(10._8))/(rnxg*rnyg*rnzg)
     write(*,*)'--- summary ---'
     write(*,'(A,F8.3,A)')"time spent to solve :",tend-tstart," s"
     write(*,'(A,E10.3)')"rescaled performance:",perf
     write(*,*)'---------------'
  end if



  call mpi_finalize(ierr)

  call toc(1,'mg_testrelax')
  if(myrank == 0) call print_tictoc(myrank)

end program mg_testrelax
