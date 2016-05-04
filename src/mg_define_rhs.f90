module mg_define_rhs

  use mg_mpi
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
!-------------------------------------------------------------------------     
  subroutine rhs_seamount()

    integer(kind=4):: nx,ny,nz
    integer(kind=4):: i,j,k
    integer(kind=4):: pi, pj
    integer(kind=4):: npxg

    real(kind=8) :: bet, x, x1, z1, x2, z2

    real(kind=8), dimension(:,:,:), pointer :: rhs

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    npxg = grid(1)%ngx

    pj = myrank / npxg   
    pi = myrank - pj * npxg

    ! rhs definition
    bet = 1200._8 / (Lx*Lx)
    x1  = Lx   * 0.45_8
    z1  = Htot * (0.75_8 - 1._8)
    x2  = Lx   * 0.75_8
    z2  = Htot * (0.65_8 - 1._8)

    rhs => grid(1)%b

    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1 
          x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(j,i)
          do k = 1,nz
             rhs(k,j,i) = dx(j,i)*dy(j,i)*(zw(k+1,j,i)-zw(k,j,i)) * &
                  (exp(-bet * ((x-x1)**2 + (zr(k,j,i)-z1)**2)) - &
                  exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))
             !           rhs(k,j,i) = (exp(-bet * ((x-x1)**2 + (zr(k,j,i)-z1)**2)) - &
             !                         exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))
             !           rhs(k,j,i) = rhs(k,j,i) * rmask(j,i)
          enddo
       enddo
    enddo

    call fill_halo(1,rhs)

  end subroutine rhs_seamount

  !-------------------------------------------------------------------------     
  subroutine rhs_random()

    integer(kind=4):: nx,ny,nz
    integer(kind=4):: i,j,k

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    call random_number(grid(1)%b)

    grid(1)%b = 2._8 * grid(1)%b - 1._8

    call fill_halo(1,grid(1)%b)

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             grid(1)%b(k,j,i) = &
                  grid(1)%dx(j,i)*grid(1)%dy(j,i)* &
                  ( grid(1)%zw(k+1,j,i) - grid(1)%zw(k,j,i) ) * grid(1)%b(k,j,i)
          enddo
       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='rhs',netcdf_file_name='rhs_random.nc',rank=myrank)
    endif

  end subroutine rhs_random

  !-------------------------------------------------------------------------     
  subroutine setup_random_patches()

    ! define the r.h.s. as a sum of random gaussian patches
    ! sign, width, location are random

    integer(kind=4):: nx,ny,nz
    integer(kind=4):: i,j,k,l
    integer(kind=4):: npxg,npyg,pi,pj
    integer(kind=4):: nbpatch,i0,j0,k0,sign
    real(kind=8)   :: sigh2,sigv2,dh2,dv2,x0,y0,z0

    nbpatch = 100

    sigh2=10.**2
    sigv2=3.**2

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    grid(1)%b = 0._8
    do l=1,nbpatch

       call random_number(x0)
       sign=floor(x0*2)*2-1

       call random_number(x0)
       sigh2=(x0*50)**2
       call random_number(x0)
       sigv2=(x0*20)**2

       call random_number(x0)
       call random_number(y0)
       call random_number(z0)
       ! we use the (i,j,k) global coordinates
       ! instead of (x,y,z)
       i0=1+floor(x0*npxg*nx)
       j0=1+floor(y0*npyg*ny)
       k0 = 1+floor(z0*nz)
       !!if(myrank==0)write(*,10)i0,j0,k0,sign,sqrt(sigh2),sqrt(sigv2)
       do i = 0,nx+1
          do j = 0,ny+1 
             do k = 1,nz
                dh2 = 1.*(i+pi*nx-i0)**2+1.*(j+pj*ny-j0)**2
                dv2 = 1.*(k-k0)**2
                grid(1)%b(k,j,i) = grid(1)%b(k,j,i) +&
                     sign*exp( -dh2 / (2.*sigh2) -dv2/(2.*sigv2))
             enddo
          enddo
       enddo
    enddo

10  format("(i0,j0,k0) =",I5,I5,I4," - sign,sigh,sigv = ",I3,F6.1,F6.1)

    call fill_halo(1,grid(1)%b)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='rhs',netcdf_file_name='rhs.nc',rank=myrank)
    endif

  end subroutine setup_random_patches


end module mg_define_rhs
