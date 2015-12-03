module mg_define_matrix

  use mg_tictoc
  use mg_grids

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine define_matrices()

    integer(kind=is)::  lev

    lev = 1
    call define_matrix_simple(lev)

    do lev=1, nlevs-1
       call coarsen_matrix(lev)
    enddo

  end subroutine define_matrices


  !----------------------------------------
  subroutine define_matrix_simple(lev)

    integer(kind=is),intent(in):: lev

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

    real(kind=rl), dimension(:,:,:,:), pointer :: cA
    integer(kind=is):: k, j, i
    real(kind=rl):: dxi, dyi, dzi
    integer(kind=is):: nx, ny, nz
    integer(kind=is):: nh

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
    integer(kind=is),intent(in):: lev

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
    integer(kind=is),intent(in):: lev
    write(*,*)'Error: coarsen matrix aggressive not available yet !'
    stop -1
  end subroutine coarsen_matrix_aggressive


  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_2D(lev)
    integer(kind=is),intent(in):: lev

  end subroutine coarsen_matrix_2D

  !-------------------------------------------------------------------------     
  subroutine coarsen_matrix_3D(lev) ! from lev to lev+1

    integer(kind=is),intent(in):: lev

    real(kind=rl), dimension(:,:,:,:), pointer :: cA, cA2
    integer(kind=is):: l, k, j, i, kp, jp, ip, k2, j2, i2
    integer(kind=is):: d
    integer(kind=is):: nx2, ny2, nz2, nh
    real(kind=rl):: diag,cff

    call tic(lev,'coarsen_matrix_3D')

    cA  => grid(lev)%cA
    cA2 => grid(lev+1)%cA
    nx2 = grid(lev+1)%nx
    ny2 = grid(lev+1)%ny
    nz2 = grid(lev+1)%nz
    nh  = grid(lev+1)%nh

    do i2 = 1,nx2
       i = 2*i2-1
       ip = i+1
       do j2 = 1,ny2
          j = 2*j2-1
          jp = j+1
          do k2 = 1,nz2
             k = 2*k2-1
             kp = k+1
             ! cA(2,:,:,:)      -> p(k-1,j,i)
             cA2(2,k2,j2,i2) = cA(2,k,j,i)+cA(2,k,jp,i)+cA(2,k,j,ip)+cA(2,k,jp,ip)
             ! cA(3,:,:,:)      -> p(k+1,j-1,i)
             cA2(3,k2,j2,i2) = cA(3,k,j,i)+cA(3,k,j,ip)
             ! cA(4,:,:,:)      -> p(k,j-1,i)
             cA2(4,k2,j2,i2) = cA(4,k,j,i)+cA(4,kp,j,i)+cA(4,k,j,ip)+cA(4,kp,j,ip)
             ! cA(5,:,:,:)      -> p(k-1,j-1,i)
             cA2(5,k2,j2,i2) = cA(5,k,j,i)+cA(5,k,j,ip)
             ! cA(6,:,:,:)      -> p(k+1,j,i-1)
             cA2(6,k2,j2,i2) = cA(6,k,j,i)+cA(6,k,jp,i)
             ! cA(7,:,:,:)      -> p(k,j,i-1)
             cA2(7,k2,j2,i2) = cA(7,k,j,i)+cA(7,kp,j,i)+cA(7,k,jp,i)+cA(7,kp,jp,i)
             ! cA(8,:,:,:)      -> p(k-1,j,i-1)
             cA2(8,k2,j2,i2) = cA(8,k,j,i)+cA(8,k,jp,i)

             ! the diagonal term is the sum of 48 terms ...             
             ! why?
             ! easy to see: the coarse cell (call it the box) is made of 8 fine cells
             ! take one fine cell, lay the 15 points on this cell
             ! count how many fine cells within this box are connected to it
             ! you should find 6
             ! multiply that by the number of fine cells

             ! here is the first 20
             diag = cA(2,kp,j,i)+cA(2,kp,jp,i)+cA(2,kp,j,ip)+cA(2,kp,jp,ip)
             diag = cA(3,k,jp,i)+cA(3,k,jp,ip)                              + diag
             diag = cA(4,k,jp,i)+cA(4,kp,jp,i)+cA(4,k,jp,ip)+cA(4,kp,jp,ip) + diag
             diag = cA(5,kp,j,ip)+cA(5,kp,jp,ip)                            + diag
             diag = cA(6,k,j,ip)+cA(6,k,jp,ip)                              + diag
             diag = cA(7,k,j,ip)+cA(7,kp,j,ip)+cA(7,k,jp,ip)+cA(7,kp,jp,ip) + diag
             diag = cA(8,kp,j,ip)+cA(8,kp,jp,ip)                            + diag

             ! double that to account for symmetry of connections, we've now 40 terms
             diag = diag+diag

             ! add the 8 self-interacting terms
             diag = cA(1,k,j,i) +cA(1,kp,j,i) +cA(1,k,jp,i) +cA(1,kp,jp,i) &
                  +cA(1,k,j,ip)+cA(1,kp,j,ip)+cA(1,k,jp,ip)+cA(1,kp,jp,ip) + diag

             ! here we go!
             cA2(1,k2,j2,i2) = diag
          enddo
       enddo
    enddo

    if (myrank.eq.0)write(*,*)"coefficients computed"

    ! the coefficients should be rescaled with 1/4

    cff = 1._8/4._8

    do d = 1,8       
       if (myrank.eq.0)write(*,*)"updating halo of coef(",l,",:,:,:)"

       !NG !? call fill_halo(lev+1, cA2(d,:,:,:))

       cA2(d,:,:,:) =  cA2(d,:,:,:) * cff

       ! a way of improvement (only if it impacts perfs):
       ! copy from cA2 to dummy3 only the interior ring used to fill the halo
       ! copy from dummy to cA2 only the halo
    enddo
    if (myrank.eq.0)write(*,*)"coarsening done"

    call toc(lev,'coarsen_matrix_3D')

  end subroutine coarsen_matrix_3D

end module mg_define_matrix
