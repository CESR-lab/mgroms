      program mg_test


      implicit none

      integer:: nx,ny,nz,np,mp,nsmall
      real*8:: aspect_ratio

      nx=128
      ny=128
      nz=128
      aspect_ratio = 0.2
      np = 16
      mp = 16
      nsmall = 32

      call test_grids(nx,ny,nz,np,mp,aspect_ratio,nsmall)

      call test_restrict()

      end

      !----------------------------------------
      subroutine test_restrict()
      use mg_grids
      use mg_optimized

      integer::i,j,k,lev
      
      lev=6

      do k=1,grid(lev)%nz
         do j=1,grid(lev)%ny
            do i=1,grid(lev)%nx
               grid(lev)%x(i,j,k)=(i-nhalo-1)/2+(j-nhalo-1)/2+(k-1)/2
            enddo
         enddo
      enddo
      
      if (grid(lev)%coarsening_method.eq.1)then
         write(*,*)'regular 3D reduction'
         call restrict_xyz(lev,lev+1,grid(lev)%x,grid(lev+1)%x)
      endif
      if (grid(lev)%coarsening_method.eq.2)then
         write(*,*)'horizontal 2D reduction'
         call restrict_xy(lev,lev+1,grid(lev)%x,grid(lev+1)%x)
      endif
      if (grid(lev)%coarsening_method.eq.5)then
         write(*,*)'triple vertical reduction'
         call restrict_zzz(lev,lev+1,grid(lev)%x,grid(lev+1)%x)
      endif
      
      lev=lev
      write(*,*)"fine level",lev
      do k=1,min(2,grid(lev)%nz)
         write(*,"(A,I2,A)")"------------k=",k,"-------------"
         do j=nhalo+1,nhalo+8
            write(*,"(8F5.1)")grid(lev)%x(nhalo+1:nhalo+8,j,k)
         enddo
      enddo

      write(*,*)"coarse level"
      lev=lev+1
      do k=1,1
         write(*,"(A,I2,A)")"------------k=",k,"-------------"
         do j=nhalo+1,nhalo+4
            write(*,"(4F5.1)")grid(lev)%x(nhalo+1:nhalo+4,j,k)
         enddo
      enddo

      end subroutine

      !----------------------------------------
      subroutine test_grids(nx,ny,nz,np,mp,aspect_ratio,nsmall)
      !
      use mg_grids

      integer:: nx,ny,nz,np,mp,nsmall
      real*8:: aspect_ratio
      !
      call define_grids(nx,ny,nz,np,mp,aspect_ratio,nsmall)
      call print_grids()
      !
      end subroutine
