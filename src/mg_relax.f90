module mg_relax

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_define_matrix

  implicit none

contains

  !----------------------------------------
  subroutine relax(lev,nsweeps)
    integer(kind=ip), intent(in):: lev
    integer(kind=ip), intent(in):: nsweeps

    real(kind=rp),dimension(:,:,:), pointer:: p
    real(kind=rp),dimension(:,:,:), pointer:: b
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nh, nd

    p  => grid(lev)%p
    b  => grid(lev)%b
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call relax_2D_5(lev,p,b,cA,nsweeps,nx,ny,nh)

    else

       call relax_line_3D_8(lev,p,b,cA,nsweeps,nx,ny,nz,nh)

    end if

  end subroutine relax

  !----------------------------------------
  subroutine relax_2D_5(lev,p,b,cA,nsweeps,nx,ny,nh)

    integer(kind=ip)                         , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nh

    integer(kind=ip)           :: i,j,k, it,rb
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp) :: z,gamma,g1,g2

    gamma = 1._8
    g1 = gamma
    g2 = 1._8 - gamma

    k=1

    do it = 1,nsweeps
       if (mod(it,nh) == 0) then
          ib = 1 
          ie = nx
          jb = 1
          je = ny
       else
          ib = 0
          ie = nx+1
          jb = 0
          je = ny+1
       endif

       if (red_black) then
          rbb = 1
          rbe = 2
          rbi = 2
       else
          rbb = 0
          rbe = 0
          rbi = 1
       endif
       rbb = 0
       rbe = 0
       rbi = 1

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

                z =    b(k,j,i)                                           &
                     - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
                     - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
                     - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
                     - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

                p(k,j,i) = z / cA(1,k,j,i)

             enddo
          enddo
       enddo

       call fill_halo(lev,p)

    enddo

  end subroutine relax_2D_5


  !----------------------------------------
  subroutine relax_line_3D_8(lev,p,b,cA,nsweeps,nx,ny,nz,nh)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz, nh

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

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
    !     LOCAL 
    integer(kind=ip)            :: i,j,k,it,rb
    integer(kind=ip)            :: ib,ie,jb,je
    integer(kind=ip)            :: rbb,rbe,rbi
!!$    integer(kind=ip)            :: fc1, fc2, fcb, fce, fci
    real(kind=rp),dimension(nz) :: rhs,d,ud

    real(kind=rp)    :: g1,g2,gamma

    call tic(lev,'relax_line')

    ! monitor convergence
    !    res0 = sum(b(1:nz,1:ny,1:nx)**2)
    !    call global_sum(lev,res0,bnorm)
    !    call compute_residual(lev,rnorm)
    !    res0 = rnorm/bnorm

    gamma = 1.8_8
    g1 = gamma
    g2 = 1._8 - gamma

    ! add a loop on smoothing
    do it = 1,nsweeps

       if (mod(it,nh) == 0) then
          ib = 1 
          ie = nx
          jb = 1
          je = ny
       else
          ib = 0
          ie = nx+1
          jb = 0
          je = ny+1
       endif

       if (red_black) then
          rbb = 1
          rbe = 2
          rbi = 2
       else
          rbb = 0
          rbe = 0
          rbi = 1
       endif

       do rb = rbb,rbe ! Red black loop

          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

!!$       ib = 1
!!$       ie = nx
!!$       jb = 1
!!$       je = ny
!!$
!!$       fcb = 1
!!$       fce = 2
!!$       fci = 2
!!$
!!$       do fc1 = fcb, fce ! 
!!$
!!$          do fc2 = fcb, fce ! 
!!$
!!$             do i = ib+mod(fc1-1,fci),ie,fci
!!$
!!$                do j = jb+mod(fc2-1,fci),je,fci

                   k=1 !lower level
                   rhs(k) = b(k,j,i)                                              &
                        - cA(3,k,j,i)*p(k+1,j-1,i)                                &
                        - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
                        - cA(5,k+1,j+1,i)*p(k+1,j+1,i) &
                        - cA(6,k,j,i)*p(k+1,j,i-1)                                &
                        - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
                        - cA(8,k+1,j,i+1)*p(k+1,j,i+1) 

                   if (cmatrix == 'real') then
                      !- Exception for the redefinition of the coef for the bottom level
                      rhs(k) = rhs(k) &
                           - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
                           - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)
                   endif

                   d(k)   = cA(1,k,j,i)
                   ud(k)  = cA(2,k+1,j,i)

                   do k = 2,nz-1 !interior levels
                      rhs(k) = b(k,j,i) &
                           - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
                           - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
                           - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i) &
                           - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
                           - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
                           - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1) 
                      d(k)   = cA(1,k,j,i)
                      ud(k)  = cA(2,k+1,j,i)
                   enddo

                   k=nz !upper level
                   rhs(k) = b(k,j,i)                                              &
                        - cA(3,k-1,j+1,i)*p(k-1,j+1,i) &
                        - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i) &
                        - cA(5,k,j,i)*p(k-1,j-1,i)                                &
                        - cA(6,k-1,j,i+1)*p(k-1,j,i+1) &
                        - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1) &
                        - cA(8,k,j,i)*p(k-1,j,i-1) 
                   d(k)   = cA(1,k,j,i)

                   call tridiag(nz,d,ud,rhs,p(:,j,i)) !solve for vertical_coeff_matrix.p1d=rhs

                   ! December 10th, dev below is to try to by-pass the computation of the residual
                   ! in the dedicated subroutine and to use instead the rhs(k) to compute the
                   ! residual on the fly.
                   ! it seems that what we gain on the one hand, less computation, we lose it on 
                   ! the other hand, worse convergence rate. Overall the rescaled convergence time
                   ! seems to be the same
                   ! obviously, this residual is only an estimate, it is not the correct one

                enddo ! j
             enddo    ! i

             ! don't call mpi at every pass if nh>1
             !if ((mod(it,nh) == 0).or.(it==nsweeps)) then
             call fill_halo(lev,p)
             !endif

!!$          enddo  ! fc2
!!$
!!$       enddo  ! fc1

       enddo ! Red Black

    enddo  !it

    call toc(lev,'relax_line')

  end subroutine relax_line_3D_8

  !----------------------------------------
  subroutine tridiag(l,d,dd,b,xc)
    !     Axc = b
    !     Solve tridiagonal system
    implicit none
    !     IMPORT/EXPORT
    integer                   ,intent(in)  :: l
    real(kind=rp),dimension(l),intent(in)  :: d,b
    real(kind=rp),dimension(l),intent(in)  :: dd
    real(kind=rp),dimension(l),intent(out) :: xc
    !     LOCAL
    integer                  :: k
    real(kind=rp),dimension(l):: gam
    real(kind=rp)             :: bet

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
    !    endif
  end subroutine tridiag

  !----------------------------------------
  subroutine compute_residual(lev,res)
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)   , intent(out):: res

    real(kind=rp),dimension(:,:,:)  , pointer:: p
    real(kind=rp),dimension(:,:,:)  , pointer:: b
    real(kind=rp),dimension(:,:,:)  , pointer:: r
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nh, nd
    real(kind=rp) ::resloc

    p  => grid(lev)%p
    b  => grid(lev)%b
    r  => grid(lev)%r
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    else

       call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    end if

    call fill_halo(lev,r)

    if (lev >-1) then
       resloc=res
       call global_sum(lev,resloc,res)
       res = sqrt(res)
    else
       res = -999._8
    endif

  end subroutine compute_residual

  !----------------------------------------
  subroutine compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k
    real(kind=rp)  :: z

    res = 0._8

    k=1

    do i = 1,nx
       do j = 1,ny

          z = b(k,j,i) - cA(1,k,j,i)*p(k,j,i)                           &
               - cA(2,k,j,i)*p(k,j-1,i  ) - cA(2,k,j+1,i  )*p(k,j+1,  i)&
               - cA(3,k,j,i)*p(k,j  ,i-1) - cA(3,k,j  ,i+1)*p(k,j  ,i+1)&
               - cA(4,k,j,i)*p(k,j-1,i-1) - cA(4,k,j+1,i+1)*p(k,j+1,i+1)&
               - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1)

          r(k,j,i) = z
          !          res = max(res,abs(r(k,j,i)))
          res = res+z*z

       enddo
    enddo

  end subroutine compute_residual_2D_5

  !----------------------------------------
  subroutine compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip)           :: i,j,k

    res = 0._8

    do i = 1,nx
       do j = 1,ny

          k=1 !lower level
          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k+1,j,i)*p(k+1,j,i)                               &
               - cA(3,k,j,i)*p(k+1,j-1,i)                               &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                           &
               - cA(6,k,j,i)*p(k+1,j,i-1)                               &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k+1,j,i+1)*p(k+1,j,i+1)

          if (cmatrix == 'real') then
             !- Exception for the redefinition of the coef for the bottom level
             r(k,j,i) = r(k,j,i) &
                  - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
                  - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)
          endif

          res = res+r(k,j,i)*r(k,j,i)

          do k = 2,nz-1 !interior levels
             r(k,j,i) = b(k,j,i)                                           &
                  - cA(1,k,j,i)*p(k,j,i)                                   &
                  - cA(2,k,j,i)*p(k-1,j,i)   - cA(2,k+1,j,i)*p(k+1,j,i)    &
                  - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
                  - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i)&
                  - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
                  - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                  - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1)

             res = res+r(k,j,i)*r(k,j,i)
          enddo

          k=nz !upper level
          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k,j,i)*p(k-1,j,i)                                 &
               - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k,j,i)*p(k-1,j-1,i)                               &
               - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k,j,i)*p(k-1,j,i-1)

          res = res+r(k,j,i)*r(k,j,i)

       enddo
    enddo

  end subroutine compute_residual_3D_8

end module mg_relax
