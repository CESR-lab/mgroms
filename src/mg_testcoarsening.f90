program mg_testcoarsening

  use mg_mpi ! everything will come from the outside !!!

  use mg_grids
  use mg_define_rhs
  use mg_define_matrix
  use mg_relax
  use mg_restrict

  implicit none

  integer(kind=is):: nxg    ! global x dimension
  integer(kind=is):: nyg    ! global y dimension
  integer(kind=is):: nzg    ! z dimension
  !!integer(kind=is):: nhalo  ! number of halo points
  integer(kind=is):: npxg   ! number of processes in x
  integer(kind=is):: npyg   ! number of processes in y
  integer(kind=is):: it     ! iteration loop number
  integer(kind=is):: nit    ! number of iterations

  integer(kind=is):: nsweeps

  integer(kind=is):: lev, ierr
  real(kind=8)    :: res,res0,conv

  nxg   = 128
  nyg   = 128
  nzg   = 128
  nhalo = 1

  npxg  = 2
  npyg  = 2

  nit     = 10
  nsweeps = 1

  call init_mpi(nxg, nyg, nzg, npxg, npyg)

  call define_grids(nhalo,npxg,npyg)


  call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  if (myrank.eq.0)then
    do lev=1,nlevs
       if (grid(lev)%gather.eq.0)then
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%npx,' x',grid(lev)%npy," procs"
       else
          write(*,100)"lev=",lev,": ", &
                     grid(lev)%nx,' x',grid(lev)%ny,' x',grid(lev)%nz, &
                     " on ",grid(lev)%npx,' x',grid(lev)%npy," procs / gather"
       endif
    enddo
 endif
100 format (A4,I2,A,I3,A,I3,A,I3,A,I3,A,I3,A)



  call define_rhs(nxg, nyg, npxg)

  write(*,*)"rhs  done"

  call define_matrix_simple()

  write(*,*)"matrix simple done"

  ! coarsen matrix on all grids
  do lev=1,nlevs-1
     call coarsen_matrix(lev)
  enddo
  write(*,*)"coarsening done"
  do lev=1,nlevs
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
     if (myrank.eq.0)then
        write(*,'(A,I2,A,8F6.3)')"stencil on lev =",lev," / ",&
             grid(lev)%cA(:,grid(lev).nz/2,grid(lev).ny/2,grid(lev).nx/2)
!        write(*,1010)grid(lev)%cA(:,grid(lev).nz/2,grid(lev).ny/2,grid(lev).nx/2)
     endif
     call MPI_Barrier( MPI_COMM_WORLD ,ierr)
  end do
1010 format(8F6.3)

  ! coarsen RHS on all grids
  do lev=1,nlevs-1
     call restrict_xyz(lev,lev+1,grid(lev)%b,grid(lev+1)%b)
  enddo

  ! check smoothing on all grids
  do lev=1,nlevs
     if (myrank.eq.0)then
        write(*,'(A,I2,A,I3,A,I3,A,I3,A)')"check relaxation on lev = ",lev,&
             " / (",grid(lev)%nx,",",grid(lev)%ny,",",grid(lev)%nz,")"
     endif
     res0=0.
     do it=1, nit
        call relax_line(lev,nsweeps)
        call compute_residual(lev,res)
        conv = log(res0/res)/log(10.)
        if (myrank.eq.0)then
           write(*,'(A,I2,A,F6.3,A,F6.3)')"  ite=",it," - res=",res," - conv rate =",conv
        endif
        res0=res
     enddo
  enddo

!  call check_solution(lev)


  call mg_mpi_finalize()

end program mg_testcoarsening
