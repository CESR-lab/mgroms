program mg_testgather

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids
  use mg_gather

  implicit none

  integer(kind=is):: nxg    ! global x dimension
  integer(kind=is):: nyg    ! global y dimension
  integer(kind=is):: nzg    ! z dimension
  !!integer(kind=is):: nhalo  ! number of halo points
  integer(kind=is):: npxg   ! number of processes in x
  integer(kind=is):: npyg   ! number of processes in y
  integer(kind=is):: it     ! iteration loop number
  integer(kind=is):: nit    ! number of iterations

  integer(kind=is):: ngx,ngy
  integer(kind=is):: nsweeps
  integer(kind=is):: nx,ny,nz,nh

  integer(kind=is):: lev, ierr
  real(kind=8)    :: res,res0,conv

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1

  npxg  = 8
  npyg  = 8

  nit     = 10
  nsweeps = 1

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)


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

!       write(*,'(A,I2,A,F3.0)')'rank=',myrank,' / dummy(1,ny,nx)=',grid(lev)%dummy3(1,ny/2,nx/2)
       write(*,'(A,I2,A,F3.0,A,F3.0)')'rank=',myrank,' / p(1,1,1)=',grid(lev)%p(1,1,1), &
            ' / p(1,ny/2,nx/2)=',grid(lev)%p(1,ny/2,nx/2)
!       if (myrank.eq.3)write(*,*)grid(lev)%p(1,:,:)
    endif
 enddo
 

100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)
110 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A,I1,A,I1)




  call mg_mpi_finalize()

end program mg_testgather
