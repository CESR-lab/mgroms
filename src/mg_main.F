      module mg_main

      implicit none

      use mg_grids
      use mg_mpi
      use mg_smoother
      use mg_interpolation      

      contains

      !----------------------------------------
      subroutine mg_init()

      ! the model metrics should be transfered or made known to the mg
      ! module

      call define_grids()

      ! define here the operator on the finest grid
      ! I see no other way than to do it manually

      cxx=dy*dz/dx
      cyy=dx*dz/dy
      czz=dx*dy/dz
      do k=2,nz-1               ! treat separatately the lower and upper level
         do j=1,ny
            do i=1,nx
               grid(1)%A(1,i,j,k)=czz
               grid(1)%A(2,i,j,k)=cyy
               grid(1)%A(3,i,j,k)=cxx
               grid(1)%A(4,i,j,k)=-2.*(cxx+cyy+czz)
            enddo
         enddo
      enddo

      call define_mpi() ! on all grids
      call define_smoothers() ! on all grids, assuming grid(1)%A is set

      end subroutine

      !----------------------------------------
      subroutine mg_solve(x,b,tol,maxite,res,nite)
!     input: x=first guess / b=RHS
!     output: x(overwritten) = solution

      integer:: maxite,nite
      real*8:: tol,res     
      real*8,dimension(grid(1)%nx,grid(1)%ny,grid(1)%nz) :: x,b
      
      ! local
      real*8:: rnorm,bnorm,res0,conv

      call norm(1,b,bnorm) ! norm of b on level=1

      call residual(1,x,b,grid(1)%b,rnorm) ! residual returns both 'r' and its norm
      res0 = rnorm/bnorm
      
      nite=0

      do while ((nite.lt.maxite).and.(res0.gt.tol))
         call mg_Fcycle(1)
         call add_to(1,x,grid(1)%x)
         call residual(1,x,b,grid(1)%b,rnorm)
         res = rnorm/bnorm
         conv=res0/res ! error reduction after this iteration
         res0=res
         nite=nite+1
         write(*,10) nite,res,conv
      enddo

 10   format("ite = ",I," / res = ",G," / conv = ",G)

      end subroutine
      
      !----------------------------------------
      subroutine mg_Vcycle(lev1)
      
      integer:: lev1,lev
      real*8:: rnorm

      do lev=lev1,nlevs-1
         call smooth(lev,grid(lev)%x,grid(lev)%b,npre)
         call residual(lev,grid(lev)%x,grid(lev)%b,grid(lev)%r,rnorm)
         call finetocoarse(lev,lev+1,grid(lev)%r,grid(lev+1)%b)
         call set_to_zero(lev+1,grid(lev+1)%x)
      enddo

      lev=nlevs
      call smooth(lev,grid(lev)%x,grid(lev)%b,ndeepest)
      
      do lev=nlevs-1,lev1,-1
         call coarsetofine(lev+1,lev,grid(lev+1)%x,grid(lev)%r)
         call add_to(lev,grid(lev)%x,grid(lev)%r) ! add x to r
         call smooth(lev,grid(lev)%x,grid(lev)%b,npost)
      enddo
      
      end subroutine
      
      !----------------------------------------
      subroutine mg_Fcycle(lev1)

      integer:: lev1,lev

      do lev=lev1,nlevs-1
         call finetocoarse(lev,lev+1,grid(lev)%b,grid(lev+1)%b)
      enddo

      lev=nlevs
      call set_to_zero(lev,grid(lev)%x)
      call smooth(lev,grid(lev)%x,grid(lev)%b,ndeepest)

      do lev=nlevs-1,lev1,-1
         call coarsetofine(lev+1,lev,grid(lev+1)%x,grid(lev)%x)
         call mg_Vcycle(lev)
      enddo
      
      end subroutine
      
      

      end module
