module mg_define_matrix

  use mg_grids
  use mg_define_rhs

  implicit none

contains
  !----------------------------------------
  subroutine define_matrix_simple()

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

    real(kind=rl), dimension(:,:,:,:), pointer :: cA
    integer(kind=is):: k, j, i
    real(kind=rl):: dxi, dyi, dzi
    integer(kind=is):: nx, ny, nz
    integer(kind=is):: nh

    integer(kind=is):: lev=1

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    cA => grid(1)%cA ! check the syntax / lighten the writing

    dxi=1._8/dx
    dyi=1._8/dy
    dzi=1._8/dz

    !extended loops will be a pain for the real matrix
    do i = 1-nh,nx+nh
       do j = 1-nh,ny+nh
          do k = 1,nz
             cA(1,k,j,i) = 2._8*(-dxi*dxi-dyi*dyi-dzi*dzi)
             cA(2,k,j,i) = dzi*dzi
             cA(3,k,j,i) = 0.0_8
             cA(4,k,j,i) = dyi*dyi
             cA(5,k,j,i) = 0.0_8
             cA(6,k,j,i) = 0.0_8
             cA(7,k,j,i) = dxi*dxi
             cA(8,k,j,i) = 0.0_8
          enddo
          cA(1,nz,j,i) = cA(1,nz,j,i) - dzi*dzi 
          cA(1,1,j,i)  = cA(1,1,j,i)  + dzi*dzi 
       enddo
    enddo

  end subroutine define_matrix_simple
  !-------------------------------------------------------------------------     
  subroutine define_matrix

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

  end subroutine define_matrix

end module mg_define_matrix
