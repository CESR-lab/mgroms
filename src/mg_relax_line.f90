********************************************************************************
      subroutine relax(p,b,nx,ny,nz,px,py,pz,bcx,bcy,j0,j1,nu)
      implicit none
*     include 'mpi.com'
*
*     IMPORT/EXPORT
      integer                                     ,intent(in)   :: nx,ny,nz
      real(kind=8),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout):: p
      real(kind=8),dimension(nx,ny,nz)            ,intent(in)   :: b
      real(kind=8),dimension(nx+1,ny,nz)          ,intent(in)   :: px
      real(kind=8),dimension(nx,ny+1,nz)          ,intent(in)   :: py
      real(kind=8),dimension(nx,ny,nz+1)          ,intent(in)   :: pz
      integer                                     ,intent(in)   :: bcx,bcy
      integer                                     ,intent(in)   :: j0,j1
      integer                                     ,intent(in)   :: nu
*     LOCAL 
      integer                    :: i,j,red_black,iters,k
      real(kind=8),dimension(nz) :: rhs,d,ud,p1d
*
      do iters = 1, nu
        do red_black = 1,1
          do j = j0,j1
*           do i = 1 + mod(j+red_black,2),nx, 2
            do i = 1,nx
              do k = 1,nz
                rhs(k) = b(i,j,k) 
     &            - px(i,j,k)*p(i-1,j,k) - px(i+1,j,k)*p(i+1,j,k)
     &            - py(i,j,k)*p(i,j-1,k) - py(i,j+1,k)*p(i,j+1,k) 

                d(k) = -px(i,j,k)-px(i+1,j,k)-py(i,j,k)-py(i,j+1,k)-pz(i,j,k)-pz(i,j,k+1)
                ud(k) = pz(i,j,k+1)
              enddo
              d( 1) = d( 1) + pz(i,j,1 )
              d(nz) = d(nz) + pz(i,j,nz+1)  !! neumann bc at surface
                                            !! ( -sign for dirichlet)
                                            !! We can do this when
                                            !! defining pz !!
*             call tridiag(nz,d,ud,rhs,p(i,j,1:nz))
              p1d = p(i,j,1:nz)
              call tridiag(nz,d,ud,rhs,p1d)
              p(i,j,1:nz) = p1d 

            enddo
          enddo
        enddo
*         call update(p,nx,ny,nz,bcx,bcy)
!       if (j0.eq.1) then
!         p(:,0,:) = p(:,1,:)
!       endif
!       if (j1.eq.ny) then
!          p(:,ny+1,:) = p(:,ny,:)
!       endif
        do k = 1,nz
         do j = j0,j1
           p(   0,j,k) = p( 1,j,k)
           p(nx+1,j,k) = p(nx,j,k)
         enddo
        enddo
        do k = 1,nz
          do i = 0,nx+1
           p(i,   0,k) = p(i, 1,k)
           p(i,ny+1,k) = p(i,ny,k)
          enddo
         enddo
*       do j = j0,j1
        do j = 0,ny+1
          do i = 0,nx+1
           p(i,j,   0) = p( i,j,1 )
           p(i,j,nz+1) = p( i,j,nz)
         enddo
        enddo
!$OMP BARRIER
      enddo
*       
      END
******************************************************************************
      subroutine tridiag(l,d,dd,b,xc)
*     Axc = b
*     Solve tridiagonal system
      implicit none
*     IMPORT/EXPORT
      integer                  ,intent(in)  :: l
      real(kind=8),dimension(l),intent(in)  :: d,b
      real(kind=8),dimension(l),intent(in)  :: dd
      real(kind=8),dimension(l),intent(out) :: xc
*     LOCAL
      integer                  :: i,k
      real(kind=8),dimension(l):: gam
      real(kind=8)             :: bet
*
*     print *, 'hoi'
      bet   = 1._8/d(1)
      xc(1) = b(1)*bet
      do k=2,l
        gam(k)= dd(k-1)*bet
        bet     = 1._8/(d(k)-dd(k-1)*gam(k))
        xc(k) = (b(k)-dd(k-1)*xc(k-1))*bet
      enddo
      do k=l-1,1,-1
        xc(k) = xc(k)-gam(k+1)*xc(k+1)
      enddo
*     print *, 'klaar'
*
      end
