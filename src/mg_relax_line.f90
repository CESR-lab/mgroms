********************************************************************************
      subroutine relax_line(p,b,nsweeps)
      implicit none

      use poisson_coeff
      ! Coefficients are stored in order of diagonals
      ! cA(1,:,:,:)      -> p(k,j,i)
      ! cA(2,:,:,:)      -> p(k-1,j,i)
      ! cA(3,:,:,:)      -> p(k+1,j-1,i)
      ! cA(4,:,:,:)      -> p(k,j-1,i)
      ! cA(5,:,:,:)      -> p(k-1,j-1,i)
      ! cA(6,:,:,:)      -> p(k+1,j,i-1)
      ! cA(7,:,:,:)      -> p(k,j,i-1)
      ! cA(8,:,:,:)      -> p(k-1,j,i-1)
!
      real(kind=8),dimension(0:nz+1,0:ny+1,0:nx+1),intent(inout):: p
      real(kind=8),dimension(nx,ny,nz)            ,intent(in)   :: b
!     LOCAL 
      integer                    :: i,j,k,red_black
      real(kind=8),dimension(nz) :: rhs,d,ud,p1d
!
      ! add a loop on smoothing
      do red_black = 1,1
         do i = 1,nx
!           do i = 1 + mod(j+red_black,2),nx, 2
            do j = 1,ny
               k=1!lower level
               do k = 2,nz-1
                  rhs(k) = b(k,j,i) 
     &              - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i)
     &              - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)
     &              - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i)
     &              - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1)
     &              - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)
     &              - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1)
                  d(k)   = cA(1,k,j,i)
                  ud(k)  = cA(2,k+1,j,i)
                  p1d(k) = p(k,j,i)
               enddo
               k=nz
               ! ...
               call tridiag(nz,d,ud,rhs,p1d)
               do k = 1,nz
                  p(i,j,1:nz) = p1d (k)
               enddo
            enddo
         enddo
      enddo
      ! don't call mpi at every pass if nh>1
      call mpi_update(p)
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
