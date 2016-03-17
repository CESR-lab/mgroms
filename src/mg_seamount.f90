module mg_seamount


  use mg_mpi
  use mg_grids
  use mg_mpi_exchange
  use mg_define_matrix

  implicit none

    real(kind=8), dimension(:,:), pointer   :: dx, dy
    real(kind=8), dimension(:,:), pointer   :: umask, vmask, rmask
    real(kind=8), dimension(:,:), pointer   :: h
    real(kind=8), dimension(:,:,:), pointer :: zr, zw

    ! MUST BE DEFINED BY THE USER OUTSIDE OF THIS MODULE
    ! before calling setup_realistic_matrix()

    real(kind=8) :: Lx,Ly,Htot
    real(kind=8) :: hc,theta_s,theta_b

contains

  subroutine setup_realistic_matrix()

    integer(kind=4):: npxg,npyg,pi,pj
    integer(kind=4):: nx,ny,nz,nh

    integer(kind=4):: i,j,k
    real(kind=8) :: x,y,z,dist,x0,y0


    npxg = grid(1)%npx
    npyg = grid(1)%npy

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    nh = grid(1)%nh

    ! grid definition
    allocate(h(0:ny+1,0:nx+1))
    allocate(dx(0:ny+1,0:nx+1))
    allocate(dy(0:ny+1,0:nx+1))
    allocate(umask(0:ny+1,0:nx+1))
    allocate(vmask(0:ny+1,0:nx+1))
    allocate(rmask(0:ny+1,0:nx+1))
    allocate(zr(nz,0:ny+1,0:nx+1))
    allocate(zw(nz+1,0:ny+1,0:nx+1))

    dx(:,:) = Lx/real(nx*npxg,kind=rp)
    dy(:,:) = Ly/real(ny*npyg,kind=rp)


    call setup_fine_mask(nx,ny,nz,nh)

    call setup_fine_depth(nx,ny,nz,nh)
    
    call setup_scoord(nx,ny,nz,nh)

    call define_matrices(dx, dy, zr, zw, umask, vmask, rmask)


  end subroutine setup_realistic_matrix

 
  !-------------------------------------------------------------------------     
  subroutine setup_fine_depth(nx,ny,nz,nh)
    integer(kind=4):: nx,ny,nz,nh

    real(kind=8) :: x,y,x0,y0,dist
    integer(kind=4):: i,j,pi,pj,npxg,npyg

    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

    x0 = Lx * 0.5_rp
    y0 = Ly * 0.5_rp
    do i = 0,nx+1 !!!  I need to know my global index range
       do j = 0,ny+1
          x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(i,j)
          y = (real(j+(pj*ny),kind=rp)-0.5_rp) * dy(i,j)
          
          h(j,i) = Htot
          !h(j,i) = Htot * (1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5._rp)**2._rp))
          !h(j,i) = Htot * (1._rp - 0.5_rp * exp(-(x-x0)**2._rp/(Lx/5._rp)**2._rp -(y-y0)**2._rp/(Ly/5._rp)**2._rp))
          ! circular slope
          dist = sqrt((x-x0)**2._rp + (y-y0)**2._rp) / Lx
          !h(j,i) = Htot - Htot* ( (tanh( (dist -0.45)*10. )+1.)*0.5 *0.9 )

       enddo
    enddo

    call write_netcdf(h,vname='h',netcdf_file_name='h.nc',rank=myrank)

  end subroutine setup_fine_depth
  !-------------------------------------------------------------------------     
  subroutine setup_fine_mask(nx,ny,nz,nh)
    integer(kind=4):: nx,ny,nz,nh

    real(kind=8) :: dist,x,y
    integer(kind=4):: i,j,pi,pj,npxg,npyg


    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)
    
    rmask(:,:) = 0.
    do i=0,nx+1
       x=(1._8*i-0.5+pi*nx)/(npxg*nx) -0.5
       do j=0,ny+1
          y=(1._8*j-0.5+pj*ny)/(npyg*ny)-0.5
          dist = sqrt(x*x+y*y) - 0.48

          ! for a circular domain
          !if (dist<0.) then

          ! for a square domain
          if ((x>-0.5).and.(x<0.5).and.(y>-0.5).and.(y<0.5)) then
             rmask(j,i)  = 1.
          endif
       enddo
    enddo

    do i=1,nx+1     
       do j=1,ny+1
          umask(j,i)=rmask(j,i)*rmask(j,i-1)
          vmask(j,i)=rmask(j,i)*rmask(j-1,i)
       enddo
    enddo

    ! unpleasant cooking to fill_halo umask and vmask (no 2D wrapper)
    grid(1)%p(1,:,:)=umask
    grid(1)%p(2,:,:)=vmask
    call fill_halo(1,grid(1)%p)
    umask = grid(1)%p(1,:,:)
    vmask = grid(1)%p(2,:,:)
    grid(1)%p(:,:,:)=0._8
    
    call write_netcdf(umask,vname='umask',netcdf_file_name='umask.nc',rank=myrank)
    call write_netcdf(vmask,vname='vmask',netcdf_file_name='vmask.nc',rank=myrank)
    
    
  end subroutine setup_fine_mask

  !-------------------------------------------------------------------------     
  subroutine setup_rhs(nx,ny,nz,nh,rhs)
    real(kind=8), dimension(:,:,:), pointer,intent(inout)   ::  rhs
    integer(kind=4):: nx,ny,nz,nh

    integer(kind=4):: i,j,k,pi,pj,npxg,npyg
    real(kind=8) :: x,y,bet,x1,z1,x2,z2,y1,y2



    npxg = grid(1)%npx
    npyg = grid(1)%npy
    pj = myrank/npxg
    pi = mod(myrank,npxg)

  ! rhs definition
  bet = 1200._8 / (Lx*Lx)
  x1 = Lx * 0.45_8
  y1 = Ly * 0.45_8
  z1 = Htot * (0.75_8 - 1._8)
  x2 = Lx * 0.75_8
  y2 = Ly * 0.55_8
  z2 = Htot * (0.65_8 - 1._8)

  do i = 0,nx+1 !!!  I need to know my global index range
     do j = 0,ny+1 
        x = (real(i+(pi*nx),kind=rp)-0.5_rp) * dx(i,j)
        y = (real(j+(pj*ny),kind=rp)-0.5_rp) * dy(i,j)
        do k = 1,nz
           rhs(k,j,i) = dx(j,i)*dy(j,i)*(zw(k+1,j,i)-zw(k,j,i)) * &
                
                ! pseudo 2D rhs
                (exp(-bet * ((x-x1)**2 +(zr(k,j,i)-z1)**2)) - &
                 exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))

                ! 3D rhs
!                (exp(-bet * ((x-x1)**2 + (y-y1)**2 + (zr(k,j,i)-z1)**2)) - &
!                 exp(-bet * ((x-x2)**2 + (y-y2)**2 + (zr(k,j,i)-z2)**2)))


!           rhs(k,j,i) = (exp(-bet * ((x-x1)**2 + (zr(k,j,i)-z1)**2)) - &
!                         exp(-bet * ((x-x2)**2 + (zr(k,j,i)-z2)**2)))
!           rhs(k,j,i) = rhs(k,j,i) * rmask(j,i)           
        enddo
     enddo
  enddo

!  call random_number(rhs)
  do i = 0,nx+1 
     do j = 0,ny+1 
        do k = 1,nz
           rhs(k,j,i) = rhs(k,j,i) * rmask(j,i)
        enddo
     enddo
  enddo
  call fill_halo(1,rhs)

end subroutine setup_rhs


  !-------------------------------------------------------------------------     
  subroutine setup_scoord(nx,ny,nz,nh)
    ! compute zr and zw from h(j,i)

    integer(kind=4):: nx,ny,nz,nh


    integer(kind=4):: i,j,k
    real(kind=8) :: cff,cff1,cff2,hinv,sc_w,sc_r,cff_w,cff_r,z_w0,z_r0,cs_r,cs_w

    cff  = 1./float(nz)
    cff1 = 1./sinh(theta_s)
    cff2 = 0.5/tanh(0.5*theta_s)

    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             sc_w = cff*float(k-nz)
             sc_r = cff*(float(k-nz)-0.5)

             cs_r = (1.-theta_b)*cff1*sinh(theta_s*sc_r) &
                  +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5)

             cs_w = (1.-theta_b)*cff1*sinh(theta_s*sc_w) &
                  +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5)
             
             cff_w = hc * sc_w
             cff_r = hc * sc_r

             z_w0 = cff_w + cs_w*h(j,i) 
             z_r0 = cff_r + cs_r*h(j,i)  
 
             hinv = 1. / (h(j,i)+hc)

             ! roms sigma coordinates
!             zw(k,j,i) = z_w0 * (h(j,i)*hinv)
!             zr(k,j,i) = z_r0 * (h(j,i)*hinv)

             ! basic linear sigma coordinates
             zr(k,j,i) = (real(k,kind=rp)-0.5_rp)*h(j,i)/real(nz,kind=rp) - h(j,i)
             zw(k,j,i) = (real(k,kind=rp)-1.0_rp)*h(j,i)/real(nz,kind=rp) - h(j,i)
          enddo
          zw(nz+1,j,i) = 0.0_rp
       enddo
    enddo

    call write_netcdf(zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank)
    
  end subroutine setup_scoord

end module mg_seamount
