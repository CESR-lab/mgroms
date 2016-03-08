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
       if (nd==3) call relax_2D_3(lev,p,b,cA,nsweeps,nx,ny,nh)
       if (nd==9) call relax_2D_9(lev,p,b,cA,nsweeps,nx,ny,nh)
    else
       !- We can add additional 3D relax routines
       !- creteria: based on grid aspect ratio
       if (nd==8)  call relax_line_3D_8(lev,p,b,cA,nsweeps,nx,ny,nz,nh)
       if (nd==27) call relax_line_3D_27(lev,p,b,cA,nsweeps,nx,ny,nz,nh)
!       call relax_3D_alternate(lev,p,b,cA,nsweeps,nx,ny,nz,nh)
    end if

  end subroutine relax

  !----------------------------------------
  subroutine relax_2D_3(lev,p,b,cA,nsweeps,nx,ny,nh)
    integer(kind=ip)                        , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nh

    integer(kind=ip)           :: i,j,k, it

    k=1

    do it = 1,nsweeps
       do i = 1,nx
          do j = 1,ny
             p(k,j,i) = b(k,j,i)                                           &
                  - cA(2,k,j,i)*p(k  ,j-1,i) - cA(2,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(3,k,j,i)*p(k  ,j,i-1) - cA(3,k  ,j,i+1)*p(k  ,j,i+1)
             p(k,j,i) = p(k,j,i) / cA(1,k,j,i)
          enddo
       enddo

        ! don't call mpi at every pass if nh>1
       if (mod(it,nh) == 0) then
          call fill_halo(lev,p)
       endif

    enddo
    if (mod(it-1,nh) /= 0) then
       call fill_halo(lev,p)
    endif

  end subroutine relax_2D_3
  
  !----------------------------------------
  subroutine relax_2D_9(lev,p,b,cA,nsweeps,nx,ny,nh)
    integer(kind=ip)                        , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nh

    integer(kind=ip)           :: i,j,k,l,di,dj, it

    k=1

!    cA(2)(k,j,i) -> p(k+0,j+1,i+1)
!    cA(3)(k,j,i) -> p(k+0,j+1,i+0)
!    cA(4)(k,j,i) -> p(k+0,j+1,i-1)
!    cA(5)(k,j,i) -> p(k+0,j+0,i+1)
!    cA(6)(k,j,i) -> p(k+0,j+0,i-1)
!    cA(7)(k,j,i) -> p(k+0,j-1,i+1)
!    cA(8)(k,j,i) -> p(k+0,j-1,i+0)
!    cA(9)(k,j,i) -> p(k+0,j-1,i-1)
    do it = 1,nsweeps
       do i = 1,nx
          do j = 1,ny
             p(k,j,i) = b(k,j,i) 
             do l=4,11
                dj=loc(l,2) ! because loc is the location index for the 27 case
                di=loc(l,3) ! but the l=4:11 are those we need for the 9 case
                p(k,j,i) = p(k,j,i) - cA(l-2,k,j,i) * p(k,j+dj,i+di)                
             enddo
             p(k,j,i) = p(k,j,i) / cA(1,k,j,i)
          enddo
       enddo

        ! don't call mpi at every pass if nh>1
       if (mod(it,nh) == 0) then
          call fill_halo(lev,p)
       endif

    enddo
    if (mod(it-1,nh) /= 0) then
       call fill_halo(lev,p)
    endif

  end subroutine relax_2D_9

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
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp),dimension(nz) :: rhs,d,ud

    real(kind=rp)    :: rnorm,bnorm,res0,conv

    call tic(lev,'relax_line')

    ! monitor convergence
    res0 = sum(b(1:nz,1:ny,1:nx)**2)
    call global_sum(lev,res0,bnorm)
    call compute_residual(lev,rnorm)
    res0 = rnorm/bnorm

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

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

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

                !             call tic(lev,'tridiag')
                call tridiag(nz,d,ud,rhs,p(:,j,i)) !solve for vertical_coeff_matrix.p1d=rhs
                !             call toc(lev,'tridiag')

                ! December 10th, dev below is to try to by-pass the computation of the residual
                ! in the dedicated subroutine and to use instead the rhs(k) to compute the
                ! residual on the fly.
                ! it seems that what we gain on the one hand, less computation, we lose it on 
                ! the other hand, worse convergence rate. Overall the rescaled convergence time
                ! seems to be the same
                ! obviously, this residual is only an estimate, it is not the correct one
!!$             k=1
!!$             grid(lev)%r(k,j,i) = rhs(k) -cA(1,k,j,i)*p(k,j,i) &
!!$                   - cA(2,k+1,j,i)*p(k+1,j,i)
!!$             do k = 2,nz-1!interior levels
!!$                grid(lev)%r(k,j,i) = rhs(k) -cA(1,k,j,i)*p(k,j,i) &
!!$                     - cA(2,k,j,i)*p(k-1,j,i) - cA(2,k+1,j,i)*p(k+1,j,i)
!!$             enddo
!!$             k=nz
!!$             grid(lev)%r(k,j,i) = rhs(k) -cA(1,k,j,i)*p(k,j,i) &
!!$                  - cA(2,k,j,i)*p(k-1,j,i) 
             enddo
          enddo

          ! don't call mpi at every pass if nh>1
          if ((mod(it,nh) == 0).or.(it==nsweeps)) then
             call fill_halo(lev,p)
          endif

       enddo

       ! monitor convergence
       call compute_residual(lev,rnorm)
       rnorm = rnorm/bnorm
       conv = res0/rnorm
       res0 = rnorm
       if (myrank == 0) write(200,*) lev,it,rnorm,conv

    enddo

    call toc(lev,'relax_line')

  end subroutine relax_line_3D_8

  !----------------------------------------
  subroutine relax_line_3D_27(lev,p,b,cA,nsweeps,nx,ny,nz,nh)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz, nh

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    !     LOCAL 
    integer(kind=ip)            :: i,j,k,l,it,rb
    integer(kind=ip)            :: di,dj,dk
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp),dimension(nz) :: rhs,d,ud

    real(kind=rp)    :: rnorm,bnorm,res0,conv

    call tic(lev,'relax_line')

    ! monitor convergence
    res0 = sum(b(1:nz,1:ny,1:nx)**2)
    call global_sum(lev,res0,bnorm)
    call compute_residual(lev,rnorm)
    res0 = rnorm/bnorm

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

!    cA(01)(k,j,i) -> p(k+0,j+0,i+0)  ***
!    cA(02)(k,j,i) -> p(k+1,j+0,i+0)  ***
!    cA(03)(k,j,i) -> p(k-1,j+0,i+0)  *** ! cA3(k+1,j,i) == cA2(k,j,i) ?

!    cA(04)(k,j,i) -> p(k+1,j+1,i+1)
!    cA(05)(k,j,i) -> p(k+1,j+1,i+0)
!    cA(06)(k,j,i) -> p(k+1,j+1,i-1)
!    cA(07)(k,j,i) -> p(k+1,j+0,i+1)
!    cA(08)(k,j,i) -> p(k+1,j+0,i-1)
!    cA(09)(k,j,i) -> p(k+1,j-1,i+1)
!    cA(10)(k,j,i) -> p(k+1,j-1,i+0)
!    cA(11)(k,j,i) -> p(k+1,j-1,i-1)

!    cA(12)(k,j,i) -> p(k+0,j+1,i+1)
!    cA(13)(k,j,i) -> p(k+0,j+1,i+0)
!    cA(14)(k,j,i) -> p(k+0,j+1,i-1)
!    cA(15)(k,j,i) -> p(k+0,j+0,i+1)
!    cA(16)(k,j,i) -> p(k+0,j+0,i-1)
!    cA(17)(k,j,i) -> p(k+0,j-1,i+1)
!    cA(18)(k,j,i) -> p(k+0,j-1,i+0)
!    cA(19)(k,j,i) -> p(k+0,j-1,i-1)

!    cA(20)(k,j,i) -> p(k-1,j+1,i+1)
!    cA(21)(k,j,i) -> p(k-1,j+1,i+0)
!    cA(22)(k,j,i) -> p(k-1,j+1,i-1)
!    cA(23)(k,j,i) -> p(k-1,j+0,i+1)
!    cA(24)(k,j,i) -> p(k-1,j+0,i-1)
!    cA(25)(k,j,i) -> p(k-1,j-1,i+1)
!    cA(26)(k,j,i) -> p(k-1,j-1,i+0)
!    cA(27)(k,j,i) -> p(k-1,j-1,i-1)

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

                k=1 !lower level
                rhs(k) = b(k,j,i)
                do l=4,19
                   dk=loc(l,1)
                   dj=loc(l,2)
                   di=loc(l,3)
                   rhs(k) = rhs(k) - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
                enddo
                d(k)   = cA(1,k,j,i)
                ud(k)  = cA(2,k,j,i)
                
                do k = 2,nz-1 !interior levels
                   rhs(k) = b(k,j,i)
                   do l=4,27
                      dk=loc(l,1)
                      dj=loc(l,2)
                      di=loc(l,3)
                      rhs(k) = rhs(k) - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
                   enddo
                   d(k)  = cA(1,k,j,i)
                   ud(k) = cA(2,k,j,i)
                enddo

                k=nz !upper level
                rhs(k) = b(k,j,i) 
                do l=12,27
                   dk=loc(l,1)
                   dj=loc(l,2)
                   di=loc(l,3)
                   rhs(k) = rhs(k) - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
                enddo
                d(k)   = cA(1,k,j,i)

                call tridiag(nz,d,ud,rhs,p(:,j,i)) !solve for vertical_coeff_matrix.p1d=rhs
             enddo
          enddo

          ! don't call mpi at every pass if nh>1
          if ((mod(it,nh) == 0).or.(it==nsweeps)) then
             call fill_halo(lev,p)
          endif

       enddo

       ! monitor convergence
       call compute_residual(lev,rnorm)
       rnorm = rnorm/bnorm
       conv = res0/rnorm
       res0 = rnorm
       if (myrank == 0) write(200,*) lev,it,rnorm,conv

    enddo

    call toc(lev,'relax_line')


  end subroutine relax_line_3D_27

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
    
  end subroutine tridiag

  !----------------------------------------
  subroutine relax_3D_alternate(lev,p,b,cA,nsweeps,nx,ny,nz,nh)

    integer(kind=ip)                        , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer              :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz, nh

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
    integer(kind=ip)           :: i,j,k,it,nh2
    real(kind=rp),dimension(nz) :: rhs,d,ud

    real(kind=rp) :: c1,c2,c3,z

    call tic(lev,'relax_alternate')

    !--DEBUG !!!!
!    call fill_halo(lev,p)
    !--DEBUG !!!!

    ! under/over relaxation parameter (c1==1 is the regular Jacobi)
    c1 = 0.8_8
    c2 = 1. - c1

!    r=>grid(lev)%p

    !
    ! add a loop on smoothing

    nh2 = nh-1
    do it = 1,nsweeps

       do i = 1-nh2,nx+nh2
          !           do i = 1 + mod(j+red_black,2),nx, 2
          do j = 1-nh2,ny+nh2

          k=1!lower level
  !        if(cA(1,k,j,i)/=0.)then
             c3 = c1/cA(1,k,j,i)
  !        else
  !           c3 = 0.
  !        endif
                

          p(k,j,i) = p(k,j,i)*c2 + c3*( b(k,j,i)                        &
               - cA(2,k+1,j,i)*p(k+1,j,i)                               &
               - cA(3,k,j,i)*p(k+1,j-1,i)                               &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k+1,j+1,i)*p(k+1,j+1,i)                           &
               - cA(6,k,j,i)*p(k+1,j,i-1)                               &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k+1,j,i+1)*p(k+1,j,i+1) )


          do k = 2,nz-1!interior levels
 !            if(cA(1,k,j,i)/=0.)then
                c3 = c1/cA(1,k,j,i)
 !            else
 !               c3 = 0.
 !            endif
             
             p(k,j,i) = p(k,j,i)*c2 + c3*( b(k,j,i)                        &
                  - cA(2,k,j,i)*p(k-1,j,i)   - cA(2,k+1,j,i)*p(k+1,j,i)    &
                  - cA(3,k,j,i)*p(k+1,j-1,i) - cA(3,k-1,j+1,i)*p(k-1,j+1,i)&
                  - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                  - cA(5,k,j,i)*p(k-1,j-1,i) - cA(5,k+1,j+1,i)*p(k+1,j+1,i)&
                  - cA(6,k,j,i)*p(k+1,j,i-1) - cA(6,k-1,j,i+1)*p(k-1,j,i+1)&
                  - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                  - cA(8,k,j,i)*p(k-1,j,i-1) - cA(8,k+1,j,i+1)*p(k+1,j,i+1))

          enddo

          k=nz!upper level
!          if(cA(1,k,j,i)/=0.)then
             c3 = c1/cA(1,k,j,i)
!          else
!             c3 = 0.
!          endif

          p(k,j,i) = p(k,j,i)*c2 + c3*( b(k,j,i)                        &
               - cA(2,k,j,i)*p(k-1,j,i)                                 &
               - cA(3,k-1,j+1,i)*p(k-1,j+1,i)                           &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(5,k,j,i)*p(k-1,j-1,i)                               &
               - cA(6,k-1,j,i+1)*p(k-1,j,i+1)                           &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
               - cA(8,k,j,i)*p(k-1,j,i-1) )
          
          enddo
       enddo
       
       if((mod(it,nh)==0).or.(it==nsweeps)) call fill_halo(lev,r)

!!$       do i = 1-nh,nx+nh
!!$          do j = 1-nh,ny+nh
!!$             do k=1,nz
!!$!                c3 = cA(1,k,j,i)/c1
!!$                z  = r(k,j,i)
!!$!                r(k,j,i)=(r(k,j,i)-p(k,j,i)*c2)*c3-cA(1,k,j,i)*p(k,j,i)
!!$                p(k,j,i)=z
!!$             enddo
!!$          enddo
!!$       enddo
    enddo


    call toc(lev,'relax_alternate')

  end subroutine relax_3D_alternate

  !----------------------------------------
  subroutine compute_residual(lev,res)
    integer(kind=ip), intent(in):: lev
    real(kind=rp), intent(out):: res

    real(kind=rp),dimension(:,:,:), pointer:: p
    real(kind=rp),dimension(:,:,:), pointer:: b
    real(kind=rp),dimension(:,:,:), pointer:: r
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
       if (nd==3) call compute_residual_2D_3(res,p,b,r,cA,nx,ny)
       if (nd==9) call compute_residual_2D_9(res,p,b,r,cA,nx,ny)
    else
       call tic(lev,'compute_residual_3D')

       if (nd==8)  call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)
       if (nd==27) call compute_residual_3D_27(res,p,b,r,cA,nx,ny,nz)

       call toc(lev,'compute_residual_3D')
    end if
    
    call fill_halo(lev,r)
!    b(:,:,:) = r(:,:,:)

    if (lev >-1) then
!       call global_max(res)
       resloc=res
       call global_sum(lev,resloc,res)
    else
       res = -999._8
    endif

  end subroutine compute_residual

  !----------------------------------------
  subroutine compute_residual_2D_3(res,p,b,r,cA,nx,ny)
    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k

    res = 0._8

    k=1

    do i = 1,nx
       do j = 1,ny

          r(k,j,i) = b(k,j,i)                                           &
               - cA(1,k,j,i)*p(k,j,i)                                   &
               - cA(2,k,j,i)*p(k  ,j-1,i) - cA(2,k  ,j+1,i)*p(k  ,j+1,i)&
               - cA(3,k,j,i)*p(k  ,j,i-1) - cA(3,k  ,j,i+1)*p(k  ,j,i+1)

!          res = max(res,abs(r(k,j,i)))
          res = res+r(k,j,i)*r(k,j,i)
          
       enddo
    enddo

  end subroutine compute_residual_2D_3

  !----------------------------------------
  subroutine compute_residual_2D_9(res,p,b,r,cA,nx,ny)
    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k,l
    integer(kind=ip) :: di,dj
    real(kind=rp)              :: rr

    res = 0._8

    k=1


    do i = 1,nx
       do j = 1,ny

          rr = b(k,j,i) - cA(1,k,j,i)*p(k,j,i)    
          do l=4,11
             dj=loc(l,2) ! we use loc(4:11)
             di=loc(l,3) ! 
             rr = rr - cA(l-2,k,j,i) * p(k,j+dj,i+di)                
          enddo
          r(k,j,i) = rr
!          res = max(res,abs(r(k,j,i)))
          res = res+rr*rr
          
       enddo
    enddo

  end subroutine compute_residual_2D_9

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
                                          - cA(2,k+1,j,i)*p(k+1,j,i)    &
               - cA(3,k,j,i)*p(k+1,j-1,i)                               &
               - cA(4,k,j,i)*p(k  ,j-1,i) - cA(4,k  ,j+1,i)*p(k  ,j+1,i)&
                                          - cA(5,k+1,j+1,i)*p(k+1,j+1,i)&
               - cA(6,k,j,i)*p(k+1,j,i-1)                               &
               - cA(7,k,j,i)*p(k  ,j,i-1) - cA(7,k  ,j,i+1)*p(k  ,j,i+1)&
                                          - cA(8,k+1,j,i+1)*p(k+1,j,i+1)
          if (cmatrix == 'real') then
          !- Exception for the redefinition of the coef for the bottom level
          r(k,j,i) = r(k,j,i) &
               - cA(5,k,j,i)*p(k,j+1,i-1) - cA(5,k,j-1,i+1)*p(k,j-1,i+1) &
               - cA(8,k,j,i)*p(k,j-1,i-1) - cA(8,k,j+1,i+1)*p(k,j+1,i+1)
          endif

!          res = max(res,abs(r(k,j,i)))
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

!             res = max(res,abs(r(k,j,i)))
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

!          res = max(res,abs(r(k,j,i)))
          res = res+r(k,j,i)*r(k,j,i)
   
       enddo
    enddo

  end subroutine compute_residual_3D_8

  !----------------------------------------
  subroutine compute_residual_3D_27(res,p,b,r,cA,nx,ny,nz)
    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    integer(kind=ip)           :: i,j,k,l
    integer(kind=ip)            :: di,dj,dk
    real(kind=rp)              :: rr

    res = 0._8

    do i = 1,nx
       do j = 1,ny

          k=1 !lower level
          rr = b(k,j,i) - cA(1,k,j,i)*p(k,j,i) - cA(2,k,j,i)*p(k+1,j,i) 
          do l=4,19
             dk=loc(l,1)
             dj=loc(l,2)
             di=loc(l,3)
             rr = rr - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
          enddo
          r(k,j,i) = rr

!          res = max(res,abs(r(k,j,i)))
          res = res+rr*rr

          do k = 2,nz-1 !interior levels
             rr = b(k,j,i)
             do l=1,27
                dk=loc(l,1)
                dj=loc(l,2)
                di=loc(l,3)
                rr = rr - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
             enddo
             r(k,j,i) = rr

!             res = max(res,abs(r(k,j,i)))
             res = res+rr*rr
          enddo

          k=nz !upper level
          rr = b(k,j,i) - cA(1,k,j,i)*p(k,j,i) - cA(3,k,j,i)*p(k-1,j,i)
          do l=12,27
             dk=loc(l,1)
             dj=loc(l,2)
             di=loc(l,3)
             rr = rr - cA(l,k,j,i) * p(k+dk,j+dj,i+di)
          enddo
          r(k,j,i) = rr

!          res = max(res,abs(r(k,j,i)))
          res = res+rr*rr
   
       enddo
    enddo

  end subroutine compute_residual_3D_27

end module mg_relax
