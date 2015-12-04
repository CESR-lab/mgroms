module mg_halo

  use tictoc
  use mg_grid

  interface fill_halo
     module procedure   &
          fill_halo_3D, &
          fill_halo_4D
  end interface fill_halo

contains

   !----------------------------------------
  subroutine fill_halo_3D(lev,p)

    integer(kind=is), intent(in):: lev
    real(kind=rl), dimension(:,:,:), pointer, intent(inout)::p
 
    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: nh
    integer(kind=is) :: south, east, north, west
    integer(kind=is) :: southwest, southeast, northeast, northwest

    integer(kind=is) :: etag, wtag, ntag, stag
    integer(kind=is) :: swtag, setag, nwtag, netag
    integer(kind=is) :: ierr,status1 

    real(kind=rl), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rl), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rl), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rl), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    real(kind=rl) :: dumt

    call tic(lev,'fill_halo')

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    sendS => grid(lev)%sendS
    recvS => grid(lev)%recvS
    sendN => grid(lev)%sendN
    recvN => grid(lev)%recvN

    sendE => grid(lev)%sendE
    recvE => grid(lev)%recvE
    sendW => grid(lev)%sendW
    recvW => grid(lev)%recvW

    sendSW => grid(lev)%sendSW
    sendSE => grid(lev)%sendSE
    sendNW => grid(lev)%sendNW
    sendNE => grid(lev)%sendNE

    recvSW => grid(lev)%recvSW
    recvSE => grid(lev)%recvSE
    recvNW => grid(lev)%recvNW
    recvNE => grid(lev)%recvNE

    !-----------
    !    Fill West-East ghost cells
    !
    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh+1:nx)  
    endif

    etag = 1
    call MPI_SendRecv(                                  &
         sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,etag, &
         recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 1
    call MPI_SendRecv(                                  &
         sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,wtag, &
         recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,wtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (west.ne.MPI_PROC_NULL) then
       p(:,1:ny,1-nh:0) = recvW
    else !!Homogenous Neumann  
       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
    endif

    if (east.ne.MPI_PROC_NULL) then
       p(:,1:ny,nx+1:nx+nh) = recvE
    else !!Homogenous Neumann  
       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
    end if
    !

    !-----------
    !    Fill North-South ghost cells
    !
    !     Pack:
    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
    endif
    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh+1:ny,1:nx)  
    endif
    ntag = 1
    call MPI_SendRecv(&
         sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,ntag, &
         recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,ntag, &
         MPI_COMM_WORLD,status1,ierr)
    stag = 1
    call MPI_SendRecv(&
         sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,stag, &
         recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,stag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack:
    if (south.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1:nx)  = recvS
    else !!Homogenous Neumann  
       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
    end if

    if (north.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1:nx)  = recvN
    else!!Homogenous Neumann  
       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
    end if
    !
    !-----------
    !    Fill Corner ghost cells
    !
    ! SW and SE !
    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh+1:nx)  
    endif

    setag = 1
    call MPI_SendRecv(&
         sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,setag, &
         recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,setag, &
         MPI_COMM_WORLD,status1,ierr)

    swtag = 1
    call MPI_SendRecv(&
         sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,swtag, &
         recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,swtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (southwest.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1-nh:0) = recvSW
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,nx+1:nx+nh) = recvSE
    else !!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    end if
    !

    ! NW and NE !
    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh+1:ny,1:nh)  
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx)  
    endif

    netag = 1
    call MPI_SendRecv( &
         sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,netag, &
         recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,netag, &
         MPI_COMM_WORLD,status1,ierr)

    nwtag = 1
    call MPI_SendRecv(&
         sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,nwtag, &
         recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,nwtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (northwest.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1-nh:0) = recvNW
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    end if
    !

    !write(*,*)'rank- p(nz/2,0:1,0)       :', myrank, p(nz/2,0:1,0), p(nz/2,1,1)
    !write(*,*)'rank- p(nz/2,0:1,nx/2)    :', myrank, p(nz/2,0:1,nx/2)
    !write(*,*)'rank- p(nz/2,0:1,nx+1)    :', myrank, p(nz/2,0:1,nx+1), p(nz/2,1,nx)
    !write(*,*)'rank- p(nz/2,ny:ny+1,0)   :', myrank, p(nz/2,ny:ny+1,0), p(nz/2,ny,1)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx/2):', myrank, p(nz/2,ny:ny+1,nx/2)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx+1):', myrank, p(nz/2,ny:ny+1,nx+1), p(nz/2,ny,nx)

    call toc(lev,'fill_halo')

  end subroutine fill_halo_3D

   !----------------------------------------
  subroutine fill_halo_4D(lev,p)

    integer(kind=is), intent(in):: lev
    real(kind=rl), dimension(:,:,:,:), pointer, intent(inout)::p
 
    integer(kind=is) :: nx, ny, nz
    integer(kind=is) :: nh
    integer(kind=is) :: south, east, north, west
    integer(kind=is) :: southwest, southeast, northeast, northwest

    integer(kind=is) :: etag, wtag, ntag, stag
    integer(kind=is) :: swtag, setag, nwtag, netag
    integer(kind=is) :: ierr,status1 

    real(kind=rl), dimension(:,:,:), pointer :: sendN,recvN,sendS,recvS
    real(kind=rl), dimension(:,:,:), pointer :: sendE,recvE,sendW,recvW
    real(kind=rl), dimension(:,:,:), pointer :: sendSW,recvSW,sendSE,recvSE
    real(kind=rl), dimension(:,:,:), pointer :: sendNW,recvNW,sendNE,recvNE

    real(kind=rl) :: dumt

    call tic(lev,'fill_halo_4D')

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nh = grid(lev)%nh

    south     = grid(lev)%neighb(1)
    east      = grid(lev)%neighb(2)
    north     = grid(lev)%neighb(3)
    west      = grid(lev)%neighb(4)
    southwest = grid(lev)%neighb(5)
    southeast = grid(lev)%neighb(6)
    northeast = grid(lev)%neighb(7)
    northwest = grid(lev)%neighb(8)

    sendS => grid(lev)%sendS
    recvS => grid(lev)%recvS
    sendN => grid(lev)%sendN
    recvN => grid(lev)%recvN

    sendE => grid(lev)%sendE
    recvE => grid(lev)%recvE
    sendW => grid(lev)%sendW
    recvW => grid(lev)%recvW

    sendSW => grid(lev)%sendSW
    sendSE => grid(lev)%sendSE
    sendNW => grid(lev)%sendNW
    sendNE => grid(lev)%sendNE

    recvSW => grid(lev)%recvSW
    recvSE => grid(lev)%recvSE
    recvNW => grid(lev)%recvNW
    recvNE => grid(lev)%recvNE

    !-----------
    !    Fill West-East ghost cells
    !
    if (west.ne.MPI_PROC_NULL) then
       sendW = p(:,1:ny,1:nh)  
    endif

    if (east.ne.MPI_PROC_NULL) then
       sendE = p(:,1:ny,nx-nh+1:nx)  
    endif

    etag = 1
    call MPI_SendRecv(                                  &
         sendE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,etag, &
         recvW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,etag, &
         MPI_COMM_WORLD,status1,ierr)

    wtag = 1
    call MPI_SendRecv(                                  &
         sendW,nz*ny*nh,MPI_DOUBLE_PRECISION,west,wtag, &
         recvE,nz*ny*nh,MPI_DOUBLE_PRECISION,east,wtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (west.ne.MPI_PROC_NULL) then
       p(:,1:ny,1-nh:0) = recvW
    else !!Homogenous Neumann  
       p(:,1:ny,1-nh:0) = p(:,1:ny,nh:1:-1)
    endif

    if (east.ne.MPI_PROC_NULL) then
       p(:,1:ny,nx+1:nx+nh) = recvE
    else !!Homogenous Neumann  
       p(:,1:ny,nx+1:nx+nh) = p(:,1:ny,nx:nx-nh+1:-1)
    end if
    !

    !-----------
    !    Fill North-South ghost cells
    !
    !     Pack:
    if (south.ne.MPI_PROC_NULL) then
       sendS = p(:,1:nh,1:nx)  
    endif
    if (north.ne.MPI_PROC_NULL) then
       sendN = p(:,ny-nh+1:ny,1:nx)  
    endif
    ntag = 1
    call MPI_SendRecv(&
         sendN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,ntag, &
         recvS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,ntag, &
         MPI_COMM_WORLD,status1,ierr)
    stag = 1
    call MPI_SendRecv(&
         sendS,nz*nx*nh,MPI_DOUBLE_PRECISION,south,stag, &
         recvN,nz*nx*nh,MPI_DOUBLE_PRECISION,north,stag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack:
    if (south.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1:nx)  = recvS
    else !!Homogenous Neumann  
       p(:,1-nh:0,1:nx) = p(:,nh:1:-1,1:nx)
    end if

    if (north.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1:nx)  = recvN
    else!!Homogenous Neumann  
       p(:,ny+1:ny+nh,1:nx) = p(:,ny:ny-nh+1:-1,1:nx)
    end if
    !
    !-----------
    !    Fill Corner ghost cells
    !
    ! SW and SE !
    if (southwest.ne.MPI_PROC_NULL) then
       sendSW = p(:,1:nh,1:nh)  
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       sendSE = p(:,1:nh,nx-nh+1:nx)  
    endif

    setag = 1
    call MPI_SendRecv(&
         sendSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,setag, &
         recvSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,setag, &
         MPI_COMM_WORLD,status1,ierr)

    swtag = 1
    call MPI_SendRecv(&
         sendSW,nz*nh*nh,MPI_DOUBLE_PRECISION,southwest,swtag, &
         recvSE,nz*nh*nh,MPI_DOUBLE_PRECISION,southeast,swtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (southwest.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,1-nh:0) = recvSW
    else !!Homogenous Neumann  
       p(:,1-nh:0,1-nh:0) = p(:,nh:1:-1,nh:1:-1)
    endif

    if (southeast.ne.MPI_PROC_NULL) then
       p(:,1-nh:0,nx+1:nx+nh) = recvSE
    else !!Homogenous Neumann  
       p(:,1-nh:0,nx+1:nx+nh) = p(:,nh:1:-1,nx:nx-nh+1:-1)
    end if
    !

    ! NW and NE !
    if (northwest.ne.MPI_PROC_NULL) then
       sendNW = p(:,ny-nh+1:ny,1:nh)  
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       sendNE = p(:,ny-nh+1:ny,nx-nh+1:nx)  
    endif

    netag = 1
    call MPI_SendRecv( &
         sendNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,netag, &
         recvNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,netag, &
         MPI_COMM_WORLD,status1,ierr)

    nwtag = 1
    call MPI_SendRecv(&
         sendNW,nz*nh*nh,MPI_DOUBLE_PRECISION,northwest,nwtag, &
         recvNE,nz*nh*nh,MPI_DOUBLE_PRECISION,northeast,nwtag, &
         MPI_COMM_WORLD,status1,ierr)
    !
    !     Unpack: 
    if (northwest.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,1-nh:0) = recvNW
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,1-nh:0) = p(:,ny:ny-nh+1:-1,nh:1:-1)
    endif

    if (northeast.ne.MPI_PROC_NULL) then
       p(:,ny+1:ny+nh,nx+1:nx+nh) = recvNE
    else !!Homogenous Neumann  
       p(:,ny+1:ny+nh,nx+1:nx+nh) = p(:,ny:ny-nh+1:-1,nx:nx-nh+1:-1)
    end if
    !

    !write(*,*)'rank- p(nz/2,0:1,0)       :', myrank, p(nz/2,0:1,0), p(nz/2,1,1)
    !write(*,*)'rank- p(nz/2,0:1,nx/2)    :', myrank, p(nz/2,0:1,nx/2)
    !write(*,*)'rank- p(nz/2,0:1,nx+1)    :', myrank, p(nz/2,0:1,nx+1), p(nz/2,1,nx)
    !write(*,*)'rank- p(nz/2,ny:ny+1,0)   :', myrank, p(nz/2,ny:ny+1,0), p(nz/2,ny,1)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx/2):', myrank, p(nz/2,ny:ny+1,nx/2)
    !write(*,*)'rank- p(nz/2,ny:ny+1,nx+1):', myrank, p(nz/2,ny:ny+1,nx+1), p(nz/2,ny,nx)

    call toc(lev,'fill_halo')

  end subroutine fill_halo


end module mg_halo
