program mg_testgather

  use mg_mpi
  use mg_grids
  use mg_gather

  implicit none

  integer(kind=ip):: nxg    ! global x dimension
  integer(kind=ip):: nyg    ! global y dimension
  integer(kind=ip):: nzg    ! z dimension
  integer(kind=ip):: npxg   ! number of processes in x
  integer(kind=ip):: npyg   ! number of processes in y
  integer(kind=ip):: nit    ! number of iterations
  integer(kind=ip):: ngx,ngy
  integer(kind=ip):: nsweeps
  integer(kind=ip):: nx,ny,nz,nh
  integer(kind=ip):: lev, ierr, np, rank

  !---------------!
  !- Ocean model -!
  !---------------!
  nxg   = 128
  nyg   = 128
  nzg   = 128

  npxg  = 4
  npyg  = 4

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


  if (myrank.eq.0)write(*,*)"---------- check families ----------"
  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  do lev=1,nlevs-1
     if (grid(lev)%gather.eq.1)then
        call MPI_Barrier( MPI_COMM_WORLD ,ierr)
        if (myrank.eq.0)then
           write(*,'(A,I2,A)')"    ---  Level =",lev,' ---'
        endif
        call MPI_Barrier( MPI_COMM_WORLD ,ierr)
        write(*,'(A,I3,A,I3,A,I2,A,I1)')"rank=",myrank," /family=",grid(lev)%family," /color=",grid(lev)%color," /key=",grid(lev)%key
     endif
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  enddo

  if (myrank.eq.0)write(*,*)"---------- check gathering ----------"
  do lev=1,nlevs-1
     nx = grid(lev)%nx
     ny = grid(lev)%ny
     nz = grid(lev)%nz
     nh = grid(lev)%nh
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)

     if (grid(lev)%gather.eq.1)then

        if (myrank.eq.0)then
           write(*,'(A,I2,A)')"    ---  Level =",lev,' ---'
        endif
        call MPI_Barrier( MPI_COMM_WORLD ,ierr)

        grid(lev)%dummy3(:,:,:)=myrank*1._8
        !       write(*,*)size(grid(lev)%dummy3),grid(lev)%Ng,size(grid(lev)%p)
        call gather(lev,grid(lev)%dummy3,grid(lev)%p)
        !       call fill_halo(lev,grid(lev)%p)
        !       grid(lev)%p=myrank*1._8
        call split(lev,grid(lev)%p,grid(lev)%dummy3)


        !       write(*,'(A,I2,A,F3.0,A,F3.0)')'rank=',myrank,' / p(1,1,1)=',grid(lev)%p(1,1,1), &
        !            ' / p(1,ny/2,nx/2)=',grid(lev)%p(1,ny/2,nx/2)

        !        if(myrank==0)write(*,*)(grid(lev)%p(1,:,:))
        write(*,'(A,I2,A,F3.0)')'rank=',myrank,' / dummy3=',grid(lev)%dummy3(1,1,1)

     endif
  enddo

100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)
110 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A,I1,A,I1)

  call mpi_finalize(ierr)

end program mg_testgather
