program mg_testgalerkin

  use mpi
  use mg_namelist
  use mg_grids
  use mg_relax
  use mg_solvers
  use mg_gather
  use mg_intergrids

  implicit none
  
  ! test that <xc,Ac*xc>=<I*xc,Af*(I*xc)>

  real(kind=8) :: norm_c,norm_f,dummy
  integer(kind=4) :: lev,nx,ny,nz,nh


  do lev=nlevs,2,-1

     nx = grid(lev)%nx
     ny = grid(lev)%ny
     nz = grid(lev)%nz
     nh = grid(lev)%nh

     grid(lev)%p = 1._8
     grid(lev)%b = 0._8
     call compute_residual(lev,dummy)    
     call norm(lev-1,grid(lev)%p,grid(lev)%r,nx,ny,nz,norm_c)

     grid(lev-1)%p = 0._8 
     call coarse2fine(lev-1) ! interpolate p to r and add r to p
     grid(lev-1)%b = 0._8
     call compute_residual(lev-1,dummy)

     nx = grid(lev-1)%nx
     ny = grid(lev-1)%ny
     nz = grid(lev-1)%nz
     nh = grid(lev-1)%nh
     call norm(lev-1,grid(lev-1)%p,grid(lev-1)%r,nx,ny,nz,norm_f)


     if (grid(lev)%gather == 1) norm_f = norm_f * 4._8

     if (myrank==0) then
        write(*,*)"======== lev ",lev,"==========="
        write(*,*)"norm coarse = ",norm_c
        write(*,*)"norm fine   = ",norm_f/16._8
     endif

  enddo

end program mg_testgalerkin


subroutine norm(lev,x,y,nx,ny,nz,res)
  use mg_mpi_exchange
  integer(kind=4) :: lev,i,j,k
  integer(kind=4) :: nx,ny,nz
  real(kind=8) :: r,res
  real(kind=8),dimension(:,:,:)  , pointer :: x,y

  r=0._8
  do i=1,nx
     do j=1,ny
        do k=1,nz
           r=r+x(k,j,i)*y(k,j,i)
        enddo
     enddo
  enddo
  call global_sum(lev,r,res)
end subroutine norm

