module mg_define_matrix

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine define_matrices()

    integer(kind=ip)::  lev

    lev = 1
    call define_matrix_simple(lev)

    do lev=1, nlevs-1
       call coarsen_matrix(lev)
    enddo

  end subroutine define_matrices

  !----------------------------------------
  subroutine define_matrix_simple(lev)

    integer(kind=ip),intent(in):: lev

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)xd
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    integer(kind=ip):: k, j, i
    real(kind=rp):: dxi, dyi, dzi
    integer(kind=ip):: nx, ny, nz
    integer(kind=ip):: nh

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    cA => grid(1)%cA ! check the syntax / lighten the writing

    dxi=1._8!/dx
    dyi=1._8!/dy
    dzi=1._8!/dz

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
  subroutine coarsen_matrix(lev)
    integer(kind=ip),intent(in):: lev

    if ((aggressive).and.(lev==1)) then
       call coarsen_matrix_aggressive(lev)

    elseif (grid(lev)%nz == 1) then
       call coarsen_matrix_2D(lev)

    else
       call coarsen_matrix_3D(lev)

    end if

  end subroutine coarsen_matrix

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_aggressive(lev)
    integer(kind=ip),intent(in):: lev

    integer(kind=ip) :: idum
    idum = lev

    write(*,*)'Error: coarsen matrix aggressive not available yet !'
    stop -1
  end subroutine coarsen_matrix_aggressive


  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_2D(lev)
    integer(kind=ip),intent(in):: lev

    integer(kind=ip) :: idum
    idum = lev
    
    write(*,*)'Error: coarsen matrix_2D  not available yet !'
    stop -1

  end subroutine coarsen_matrix_2D

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D(lev) ! from lev to lev+1

    integer(kind=ip),intent(in):: lev

    real(kind=rp), dimension(:,:,:,:), pointer :: cA, cA2
    real(kind=rp), dimension(:,:,:) , pointer :: dummy3D

    integer(kind=ip):: k, j, i
    integer(kind=ip):: kf, jf, if
    integer(kind=ip):: k2, j2, i2
    integer(kind=ip):: d
    integer(kind=ip):: nx2, ny2, nz2, nh
    real(kind=rp)   :: diag,cff

    call tic(lev,'coarsen_matrix_3D')

    cA  => grid(lev)%cA
    cA2 => grid(lev+1)%cA
    nx2 = grid(lev+1)%nx
    ny2 = grid(lev+1)%ny
    nz2 = grid(lev+1)%nz
    nh  = grid(lev+1)%nh

    do i2 = 1,nx2
       i = 2*i2-1
       if = i+1
       do j2 = 1,ny2
          j = 2*j2-1
          jf = j+1     !TODO take into account GATHER case here (else=bug)
          do k2 = 1,nz2
             k = 2*k2-1
             kf = k+1
             ! cA(2,:,:,:)      -> p(k-1,j,i)
             cA2(2,k2,j2,i2) = cA(2,k,j,i)+cA(2,k,jf,i)+cA(2,k,j,if)+cA(2,k,jf,if)
             ! cA(3,:,:,:)      -> p(k+1,j-1,i)
             cA2(3,k2,j2,i2) = cA(3,k,j,i)+cA(3,k,j,if)
             ! cA(4,:,:,:)      -> p(k,j-1,i)
             cA2(4,k2,j2,i2) = cA(4,k,j,i)+cA(4,kf,j,i)+cA(4,k,j,if)+cA(4,kf,j,if)
             ! cA(5,:,:,:)      -> p(k-1,j-1,i)
             cA2(5,k2,j2,i2) = cA(5,k,j,i)+cA(5,k,j,if)
             ! cA(6,:,:,:)      -> p(k+1,j,i-1)
             cA2(6,k2,j2,i2) = cA(6,k,j,i)+cA(6,k,jf,i)
             ! cA(7,:,:,:)      -> p(k,j,i-1)
             cA2(7,k2,j2,i2) = cA(7,k,j,i)+cA(7,kf,j,i)+cA(7,k,jf,i)+cA(7,kf,jf,i)
             ! cA(8,:,:,:)      -> p(k-1,j,i-1)
             cA2(8,k2,j2,i2) = cA(8,k,j,i)+cA(8,k,jf,i)

             ! the diagonal term is the sum of 48 terms ...             
             ! why?
             ! easy to see: the coarse cell (call it the box) is made of 8 fine cells
             ! take one fine cell, lay the 15 points on this cell
             ! count how many fine cells within this box are connected to it
             ! you should find 6
             ! multifly that by the number of fine cells

             ! here is the first 20
             diag = cA(2,kf,j,i)+cA(2,kf,jf,i)+cA(2,kf,j,if)+cA(2,kf,jf,if)
             diag = cA(3,k,jf,i)+cA(3,k,jf,if)                              + diag
             diag = cA(4,k,jf,i)+cA(4,kf,jf,i)+cA(4,k,jf,if)+cA(4,kf,jf,if) + diag
             diag = cA(5,kf,j,if)+cA(5,kf,jf,if)                            + diag
             diag = cA(6,k,j,if)+cA(6,k,jf,if)                              + diag
             diag = cA(7,k,j,if)+cA(7,kf,j,if)+cA(7,k,jf,if)+cA(7,kf,jf,if) + diag
             diag = cA(8,kf,j,if)+cA(8,kf,jf,if)                            + diag

             ! double that to account for symmetry of connections, we've now 40 terms
             diag = diag+diag

             ! add the 8 self-interacting terms
             diag = cA(1,k,j,i) +cA(1,kf,j,i) +cA(1,k,jf,i) +cA(1,kf,jf,i) &
                  +cA(1,k,j,if)+cA(1,kf,j,if)+cA(1,k,jf,if)+cA(1,kf,jf,if) + diag

             ! here we go!
             cA2(1,k2,j2,i2) = diag
          enddo
       enddo
    enddo

    if (myrank.eq.0)write(*,*)"coefficients computed"

    ! fill the halo
    ! the data should be contiguous in memory to use fill_halo 
    ! no need to allocate an extra buffer
    ! use the residual as a dummy variable
    dummy3D => grid(lev+1)%r 

    !- we should consider a specific fill_halo(4D) -!

    ! the coefficients should be rescaled with 1/4
    cff = 1._8/4._8

    do d = 1,8       
       !!if (myrank.eq.0)write(*,*)"updating halo of coef(",d,",:,:,:)"
       do i2 = 1,nx2
          do j2 = 1,ny2
             do k2 = 1,nz2
                dummy3D(k2,j2,i2) = cA2(d,k2,j2,i2)
             enddo
          enddo
       enddo

       call fill_halo(lev+1,dummy3D)

       do i2 = 1-nh,nx2+nh
          do j2 = 1-nh,ny2+nh
             do k2 = 1,nz2
                cA2(d,k2,j2,i2) = dummy3D(k2,j2,i2) * cff
             enddo
          enddo
       enddo
       ! a way of improvement (only if it impacts perfs):
       ! copy from cA2 to dummy3 only the interior ring used to fill the halo
       ! copy from dummy to cA2 only the halo
    enddo

!!$    do d = 1,8       
!!$       if (myrank.eq.0)write(*,*)"updating halo of coef(",d,",:,:,:)"
!!$
!!$       dummy3D(1:nx2,1:ny2,1:nz2) = cA2(d,1:nx2,1:ny2,1:nz2)
!!$
!!$       call fill_halo(lev+1, dummy3D)
!!$
!!$       cA2(d,:,:,:) =  dummy3D(:,:,:) * cff
!!$
!!$       ! a way of improvement (only if it impacts perfs):
!!$       ! copy from cA2 to dummy3 only the interior ring used to fill the halo
!!$       ! copy from dummy to cA2 only the halo
!!$    enddo

    if (myrank.eq.0) write(*,*)"coarsening done"

    call toc(lev,'coarsen_matrix_3D')

  end subroutine coarsen_matrix_3D

end module mg_define_matrix
