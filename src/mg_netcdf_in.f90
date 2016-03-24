module mg_netcdf_in
  !*******************************************
  ! Netcdf I/O Interface					 
  !*******************************************      
  use netcdf

  implicit none 
  !      
  interface Read_Ncdf_var
     module procedure          &
          Read_Ncdf_var2d_Real, &
          Read_Ncdf_var3d_Real
  end interface Read_Ncdf_var
  !

  integer(kind=4), parameter :: ip= 4
  integer(kind=4), parameter :: rp= 8
  integer(kind=ip), parameter :: lun_error = 0

contains

  !                    
  !*******************************************
  !   subroutine Read_Ncdf_var2d_real
  !*******************************************          
  ! [] = optional argument
  subroutine Read_Ncdf_var2d_Real( &
       varname                   , & ! NetCDF variable name
       tabvar                    , & ! Array to store NetCDF variable
       file                      , & ! [NetCDF file name]
       ncid                      , & ! [NetCDF file Id]
       starts                    , & ! [Start indices in NetCDF variable]
       counts                    )   ! [Shape of array extracted from the NetCDF variable]
    !      
    ! Argument declaration
    character(len=*)                          , intent(in)    :: varname
    real(kind=rp)   , dimension(:,:), pointer , intent(out)   :: tabvar
    character(len=*)                , optional, intent(in)    :: file
    integer(kind=ip)                , optional, intent(inout) :: ncid
    integer(kind=ip),   dimension(:), optional, intent(inout) :: starts
    integer(kind=ip),   dimension(:), optional, intent(inout) :: counts

    !- local declaration
    integer(kind=ip), dimension(10) :: dimIDS
    integer(kind=ip)                :: dim1,dim2
    integer(kind=ip)                :: status
    integer(kind=ip)                :: varid

    ! Open netcdf File or verify ncid argument
    if (present(file)) then
       status = nf90_open(file,NF90_NOWRITE,ncid)      
       if (status/=nf90_noerr) then    
          write(lun_error,*)"unable to open netcdf file : ",file
          stop
       endif
    elseif (present(ncid)) then
       !! just to verify that ncid is set if file is not
    else
       write(lun_error,*)"Please set argument NetCDF file name or NetCDF Id "
       stop
    endif

    ! Get NetCDF variable Id
    status = nf90_inq_varid(ncid,varname,varid)
    if (status /= nf90_noerr) then
       write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var2d_Real'
       write(lun_error,*) 'ERROR: problem to inquire variabe Id !',varname
       stop 
    endif

    if (present(starts).and.present(counts)) then
       dim1 = counts(1)
       dim2 = counts(2)
    else
       status = nf90_inquire_variable(ncid,varid,dimids=dimIDS)
       status = nf90_inquire_dimension(ncid,dimIDS(1),len=dim1)
       status = nf90_inquire_dimension(ncid,dimIDS(2),len=dim2)
       starts = [   1,    1]
       counts = [dim1, dim2]
    endif

    ! Read netcdf 2D array
    status=nf90_get_var(   &
         ncid   = ncid   , &
         varid  = varid  , &
         values = tabvar , &
         start  = starts , &
         count  = counts )

    if (status /= nf90_noerr) then
       write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var2d_Real'
       write(lun_error,*) 'ERROR: problem to read 2D netcdf data !',ncid
       stop 
    endif

    ! Close NetCDF file
    if (present(file)) then
       status = nf90_close(ncid)

       if (status /= nf90_noerr) then
          write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var2d_Real'
          write(lun_error,*) 'ERROR: problem to close netcdf file !',file
          stop 
       endif
    endif

    !     
  end subroutine Read_Ncdf_var2d_Real
  !                    
  !*******************************************
  !   subroutine Read_Ncdf_var3d_real
  !*******************************************          
  ! [] = optional argument
  subroutine Read_Ncdf_var3d_Real( &
       varname                   , & ! NetCDF variable name
       tabvar                    , & ! Array to store NetCDF variable
       file                      , & ! [NetCDF file name]
       ncid                      , & ! [NetCDF file Id]
       starts                    , & ! [Start indices in NetCDF variable]
       counts                    )   ! [Shape of array extracted from the NetCDF variable]
    !      


    ! Argument declaration
    character(len=*)                            , intent(in)    :: varname
    real(kind=rp)   , dimension(:,:,:), pointer , intent(out)   :: tabvar
    character(len=*)                  , optional, intent(in)    :: file
    integer(kind=ip)                  , optional, intent(inout) :: ncid
    integer(kind=ip),     dimension(:), optional, intent(inout) :: starts
    integer(kind=ip),     dimension(:), optional, intent(inout) :: counts

    !- local declaration
    integer(kind=ip), dimension(10) :: dimIDS
    integer(kind=ip)                :: dim1,dim2,dim3
    integer(kind=ip)                :: status
    integer(kind=ip)                :: varid

    ! Open netcdf File or verify ncid argument
    if (present(file)) then
       status = nf90_open(file,NF90_NOWRITE,ncid)      
       if (status/=nf90_noerr) then    
          write(lun_error,*)"unable to open netcdf file : ",file
          stop
       endif
    elseif (present(ncid)) then
       !! just to verify that ncid is set if file is not
    else
       write(lun_error,*)"Please set argument NetCDF file name or NetCDF Id "
       stop
    endif

    ! Get NetCDF variable Id
    status = nf90_inq_varid(ncid,varname,varid)
    if (status /= nf90_noerr) then
       write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var3d_Real'
       write(lun_error,*) 'ERROR: problem to inquire variabe Id !',varname
       stop 
    endif

    if (present(starts).and.present(counts)) then
       dim1 = counts(1)
       dim2 = counts(2)
       dim3 = counts(3)
    else
       status = nf90_inquire_variable(ncid,varid,dimids=dimIDS)
       status = nf90_inquire_dimension(ncid,dimIDS(1),len=dim1)
       status = nf90_inquire_dimension(ncid,dimIDS(2),len=dim2)
       status = nf90_inquire_dimension(ncid,dimIDS(3),len=dim3)
       starts = [   1,    1,    1]
       counts = [dim1, dim2, dim3]
    endif

    ! Read netcdf 3D array
    status=nf90_get_var(   &
         ncid   = ncid   , &
         varid  = varid  , &
         values = tabvar , &
         start  = starts , &
         count  = counts )

    if (status /= nf90_noerr) then
       write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var2d_Real'
       write(lun_error,*) 'ERROR: problem to read 3D netcdf data !',ncid
       stop 
    endif

    ! Close NetCDF file
    if (present(file)) then
       status = nf90_close(ncid)

       if (status /= nf90_noerr) then
          write(lun_error,*) 'ERROR: io_netcdf: Read_Ncdf_var3d_Real'
          write(lun_error,*) 'ERROR: problem to close netcdf file !',file
          stop 
       endif
    endif

    !     
  end subroutine Read_Ncdf_var3d_Real

  !******************************************* 
end module mg_netcdf_in
!******************************************* 
