
      subroutine define_matrix_simple
      implicit none

      ! Define matrix coefficients cA
      ! Coefficients are stored in order of diagonals
      ! cA(1,:,:,:)      -> p(k,j,i)
      ! cA(2,:,:,:)      -> p(k-1,j,i)
      ! cA(3,:,:,:)      -> p(k+1,j-1,i)
      ! cA(4,:,:,:)      -> p(k,j-1,i)
      ! cA(5,:,:,:)      -> p(k-1,j-1,i)
      ! cA(6,:,:,:)      -> p(k+1,j,i-1)
      ! cA(7,:,:,:)      -> p(k,j,i-1)
      ! cA(8,:,:,:)      -> p(k-1,j,i-1)

      use mg_variables
!     provides the array of structure 'grid'
!     each entry of grid contains all the variables of a given level
      ! local
      integer:: k,j,i

      pointer cA <= grid(lev)%cA ! check the syntax / lighten the writing

      
      dxi=1.
      dyi=1.
      dzi=1.

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               cA(1,k,j,i) = -dxi*dxi-dyi*dyi-dzi*dzi
               cA(2,k,j,i) = dzi*dzi
               cA(3,k,j,i) = 0.0
               cA(4,k,j,i) = dyi*dyi
               cA(5,k,j,i) = 0.0
               cA(6,k,j,i) = 0.0
               cA(7,k,j,i) = dxi*dxi
               cA(8,k,j,i) = 0.0
            enddo
         enddo
       enddo

      end subroutine
!-------------------------------------------------------------------------     
      subroutine define_matrix
      implicit none

      ! Define matrix coefficients cA
      ! Coefficients are stored in order of diagonals
      ! cA(1,:,:,:)      -> p(k,j,i)
      ! cA(2,:,:,:)      -> p(k-1,j,i)
      ! cA(3,:,:,:)      -> p(k+1,j-1,i)
      ! cA(4,:,:,:)      -> p(k,j-1,i)
      ! cA(5,:,:,:)      -> p(k-1,j-1,i)
      ! cA(6,:,:,:)      -> p(k+1,j,i-1)
      ! cA(7,:,:,:)      -> p(k,j,i-1)
      ! cA(8,:,:,:)      -> p(k-1,j,i-1)

      use mg_variables
!     provides the array of structure 'grid'
!     each entry of grid contains all the variables of a given level
      ! local
      integer:: k,j,i

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz
               cA(1,k,j,i) = -dxi(j,i)*dxi(j,i)-dyi(j,i)*dyi(j,i)-dzi(k,j,i)*dzi(k,j,i)
               cA(2,k,j,i) = dzi(k,j,i)*dzi(k,j,i)
               cA(3,k,j,i) = 0.0
               cA(4,k,j,i) = dyi(j,i)*dyi(j,i)
               cA(5,k,j,i) = 0.0
               cA(6,k,j,i) = 0.0
               cA(7,k,j,i) = dxi(j,i)*dxi(j,i)
               cA(8,k,j,i) = 0.0
            enddo
         enddo
       enddo

      end subroutine
      
