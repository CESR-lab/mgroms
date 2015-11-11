      module mg_mpi

      use mg_grids

      implicit none


      type type_halo
      real*8,dimension(:,:,:,:),allocatable :: rcorners,scorners
      real*8,dimension(:,:,:,:),allocatable :: rxedges,sxedges
      real*8,dimension(:,:,:,:),allocatable :: ryedges,syedges
      integer,dimension(:),allocatable :: recv,send,rstatus,sstatus
      integer,dimension(8):: neighb
      integer:: nb,nc,n1,n2
      end type type_halo

      type type_buffer
      real*8,dimension(:,:,:,:,:),allocatable :: b
      end type type_buffer

      type(type_halo),dimension(40) :: halo

      type(type_buffer),dimension(40):: buffer

! http://www.open-mpi.org/doc/v1.8/man3/MPI_Recv_init.3.php

      contains

      !----------------------------------------
      subroutine init_mpi()

      integer:: lev
      integer:: nx,ny,nz,np,mp

      do lev=1,nlevs
         call init_fill_halo(lev,grid(lev)%neighb)
         !
         nx = grid(lev)%nx
         ny = grid(lev)%ny
         nz = grid(lev)%nz      
         !
         np = grid(lev)%np
         mp = grid(lev)%mp
         !
         if ((grid(lev)%gather).ne.0)then
            ! allocate the buffer for the gathering
            allocate(buffer(lev)%b(nx,ny,nz,np,mp))
            ! define the local communicator for the allgather 
         endif
      enddo
      end subroutine

      !----------------------------------------
      subroutine init_fill_halo(lev,neighb)

      include 'mpif.h'

      integer::lev
      integer,dimension(8)::neighb

      integer :: ierror,comm
      integer :: nb,nc,n1,n2,k,count
      integer :: nx,ny,nz
      
      nx = grid(lev)%nx
      ny = grid(lev)%ny
      nz = grid(lev)%nz      

      nb = 0 ! nb is the number of messages to send: max is 8
      nc = 0 ! corners
      n1 = 0 ! along x
      n2 = 0 ! along y
      do k=1,8
         halo(lev)%neighb(k) = neighb(k)
         if (neighb(k).ne.-1) then
            if ((k.eq.1).or.(k.eq.3).or.(k.eq.6).or.(k.eq.8)) nc = nc+1
            if ((k.eq.2).or.(k.eq.7)) n1 = n1+1
            if ((k.eq.4).or.(k.eq.5)) n2 = n2+1
            nb=nb+1
         endif
      enddo
      halo(lev)%nb = nb
      halo(lev)%nc = nc
      halo(lev)%n1 = n1
      halo(lev)%n2 = n2

      if (nb.gt.0)then
         allocate(halo(lev)%recv(nb))
         allocate(halo(lev)%send(nb))
         allocate(halo(lev)%rstatus(nb))
         allocate(halo(lev)%sstatus(nb))
      endif
      if (nc.gt.0)then
         allocate(halo(lev)%rcorners(nhalo,nhalo,nz,nc))
         allocate(halo(lev)%scorners(nhalo,nhalo,nz,nc))
      endif
      if (n1.gt.0)then
         allocate(halo(lev)%rxedges(nx-2*nhalo,nhalo,nz,n1))
         allocate(halo(lev)%sxedges(nx-2*nhalo,nhalo,nz,n1))
      endif
      if (n2.gt.0)then
         allocate(halo(lev)%ryedges(nhalo,ny-2*nhalo,nz,n2))
         allocate(halo(lev)%syedges(nhalo,ny-2*nhalo,nz,n2))
      endif


! MPI_RECV_INIT(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
! MPI_SSEND_INIT(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST,  IERROR)

      comm = MPI_COMM_WORLD

      nb = 0
      nc = 0      
      n1 = 0
      n2 = 0
      do k=1,8
         if (neighb(k).ne.-1) then
            nb = nb+1
            ! corners
            if ((k.eq.1).or.(k.eq.3).or.(k.eq.6).or.(k.eq.8)) then
               nc = nc+1
               count = nhalo*nhalo*nz
               call MPI_RECV_INIT(halo(lev)%rcorners(1,1,1,nc), count, 
     $              MPI_REAL8, 
     $              neighb(9-k), k, comm, halo(lev)%recv(nb), ierror)
               call MPI_SSEND_INIT(halo(lev)%scorners(1,1,1,nc), count, 
     $              MPI_REAL8, 
     $              neighb(k), k, comm, halo(lev)%recv(nb),  ierror)
            endif
            ! along x edges
            if ((k.eq.2).or.(k.eq.7))then
               n1 = n1+1
               count = (nx-2*nhalo)*nhalo*nz
               call MPI_RECV_INIT(halo(lev)%rxedges(1,1,1,n1), count, 
     $              MPI_REAL8, 
     $              neighb(9-k), k, comm, halo(lev)%recv(nb), ierror)
               call MPI_SSEND_INIT(halo(lev)%sxedges(1,1,1,n1), count, 
     $              MPI_REAL8, 
     $              neighb(k), k, comm, halo(lev)%recv(nb),  ierror)
            endif
            ! along y edges
            if ((k.eq.4).or.(k.eq.5)) then
               n2 = n2+1
               count = nhalo*(ny-2*nhalo)*nz
               call MPI_RECV_INIT(halo(lev)%ryedges(1,1,1,n2), count, 
     $              MPI_REAL8, 
     $              neighb(9-k), k, comm, halo(lev)%recv(nb), ierror)
               call MPI_SSEND_INIT(halo(lev)%syedges(1,1,1,n2), count, 
     $              MPI_REAL8, 
     $              neighb(k), k, comm, halo(lev)%recv(nb),  ierror)
            endif
         endif 
      enddo

      end subroutine init_fill_halo

#if defined FULLSET
      !----------------------------------------
      subroutine fill_halo(lev,x)

      ! halo(lev)%nb: nb of requests=nb of neighbours (max is 8)

      integer:: lev,ierror

      ! recv rbuffers
      call MPI_STARTALL(halo(lev)%nb, halo(lev)%recv, ierror)

      ! copy interior array to sbuffers
      call copy_to_buffers(lev,x)

      ! send sbuffers
      call MPI_STARTALL(halo(lev)%nb, halo(lev)%send, ierror)

      call MPI_WAITALL(halo(lev)%nb, halo(lev)%recv, halo(lev)%rstatus, ierror)
      call MPI_WAITALL(halo(lev)%nb, halo(lev)%send, halo(lev)%sstatus, ierror)

      ! copy rbuffers to halo
      call copy_from_buffers(lev,x)

      end subroutine


      !----------------------------------------
      subroutine fill_halo_matrix(lev,A)

      do kd=1,nd ! loop over diagonals
         dummy3d(:,:,:)=A(kd,:,:,:)
         call fill_halo(lev,dummy3d)
         A(kd,:,:,:)=dummy3d(:,:,:)
      enddo
      end subroutine

      !----------------------------------------
      subroutine copy_to_buffers(lev,x)
      
      nhalo2=nhalo*2
      do k=1,nz
         nc = 0
         n1 = 0
         n2 = 0
         if (halo(lev)%neighb(1).eq.-1)then
            nc=nc+1
            do j=1,nhalo
               do i=1,nhalo
                  halo(lev)%scorners(i,j,k,nc)=x(i+nhalo,j+nhalo,k)
               enddo
            enddo
         endif
         if (halo(lev)%neighb(2).eq.-1)then
            n1=n1+1
            do j=1,nhalo
               do i=nhalo+1,nx-nhalo
                  halo(lev)%sxedges(i-nhalo,j,k,n1)=x(i,j+nhalo,k)
               enddo
            enddo
         endif
         if (halo(lev)%neighb(3).eq.-1)then
            nc=nc+1
            do j=1,nhalo
               do i=1,nhalo
                  halo(lev)%scorners(i,j,k,nc)=x(i+nx-nhalo2,j+nhalo,k)
               enddo
            enddo
         endif
         if (halo(lev)%neighb(4).eq.-1)then
            n2=n2+1
            do j=nhalo+1,ny-nhalo
               do i=1,nhalo
                  halo(lev)%syedges(i,j-nhalo,k,n2)=x(i+nhalo,j,k)
               enddo
            enddo
         endif

            !etc
         enddo
      enddo

      end subroutine

      !----------------------------------------
      subroutine gather(x,y)


      ! 
      integer:: i,j,k,l,ii,jj

      ! b is a preallocated buffer
      call MPI_ALLGATHERV(x,b,localcomm)
      !
      !
      do l=1,mp
         jj=1+(l-1)*(m-2*nh)
         do j=1,m
            do k=1,np
               ii=1+(k-1)*(n-2*nh)
               do i=1,n                  
                  y(jj,ii) = b(l,k,j,i)
                  ii=ii+1
               enddo
            enddo
            jj=jj+1
         enddo
      enddo
      
      end subroutine
#endif     
      end module
