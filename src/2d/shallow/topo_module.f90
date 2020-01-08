! ============================================================================
!  Module for topography data
! ============================================================================
module topo_module

    use amr_module, only: xlower,xupper,ylower,yupper
    
    implicit none
    save

    logical, private :: module_setup = .false.

    ! General topography file type
    type topo_file_type

        ! Topography data
        integer :: topo_type
        real(kind=8), allocatable :: z(:, :)

        ! Geometry in space-time
        integer :: num_cells(2)
        real(kind=8) :: x_lower, x_upper, y_lower, y_upper, t_lower, t_upper
        real(kind=8) :: dx, dy

        ! Level information
        integer :: min_level, max_level

        ! Originating file
        character(len=256) :: file_path

        ! dtopo support
        logical :: save_topo = .false.

    end type topo_file_type

    ! Array of topo file data
    integer, private :: num_topo_files = 0
    integer, private, allocatable :: topo_order(:)
    type(topo_file_type), pointer, private :: topo_data(:)

    ! Analytical topography
    real(kind=8), private :: topo_left, topo_right, topo_location
    real(kind=8), private :: topo_x0,topo_x1,topo_x2,topo_basin_depth
    real(kind=8), private :: topo_shelf_depth,topo_shelf_slope,topo_beach_slope

    ! NetCDF4 support
    real(kind=4), parameter :: CONVENTION_REQUIREMENT = 1.0

    ! Missing value support
    real(kind=8) topo_missing

contains

    ! ========================================================================
    ! Read topography files as specified in topography.data
    !
    ! Each topography file has a type stored in topotype(i).
    !   topotype = 1:  standard GIS format: 3 columns: lon,lat,height(m)
    !   topotype = 2:  Header as in DEM file, height(m) one value per line
    !   topotype = 3:  Header as in DEM file, height(m) one row per line
    !   topotype = 4:  NetCDf file
    ! For other formats modify readtopo routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Finest value of topography in a given region will be used for
    ! computation
    ! ========================================================================
    subroutine read_topo_settings(file_name)

        use geoclaw_module, only: GEO_PARM_UNIT
        use dtopo_module, only: num_dtopo_files, 

        implicit none

        ! Input arguments
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: iunit = 7
        integer :: i, j, itopo, finer_than, rank
        real(kind=8) :: area_i, area_j, x_junk, y_junk
        real(kind=8) :: area, area_domain

        if (.not.module_setup) then

            ! Open and begin parameter file output
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SETTOPO:'
            write(GEO_PARM_UNIT,*) '---------'


            if (present(file_name)) then
                call opendatafile(iunit, file_name)
            else
                call opendatafile(iunit, 'topo.data')
            endif

            ! Read in value to use in place of no_data_value in topofile
            read(iunit,*) topo_missing

            ! Read in topography specification type
            read(iunit,"(i1)") test_topography

            ! Primary topography type, read in topography files specified
            if (test_topography == 0) then
                read(iunit,*) num_topo_files

                if (num_topo_files == 0) then
                    write(GEO_PARM_UNIT,*) '   num_topo_files = 0'
                    write(GEO_PARM_UNIT,*) '   No topo files specified, '
                    write(GEO_PARM_UNIT,*) '      will set B(x,y) = 0 in setaux'
                    return
                endif

                write(GEO_PARM_UNIT,*) '   num_topo_files = ', num_topo_files
                
                ! Read and allocate paramters for topo files - note that here we
                ! also allocate space for the work arrays needed for moving 
                ! topography but do not fill them in until a bit later
                allocate(topo_data(num_topo_files + num_dtopo_files))
                do i = 1, num_topo_files
                    read(iunit, *) topo_data(i)%file_path
                    read(iunit, *) topo_data(i)%topo_type,     &
                                   topo_data(i)%min_level,     &
                                   topo_data(i)%max_level,     &
                                   topo_data(i)%t_low,         &
                                   topo_data(i)%t_high

                    allocate(topo_data(i)%z(topo_data(i)%num_cells(1),   &
                                            topo_data(i)%num_cells(2)))

                    write(GEO_PARM_UNIT,*) '   '
                    write(GEO_PARM_UNIT,*) '   ',topo_data(i)%file_path
                    write(GEO_PARM_UNIT,*) '  topo_type = ',          &
                                    topo_data(i)%topo_type
                    write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                    topo_data(i)%min_level,           &
                                    topo_data(i)%max_level
                    write(GEO_PARM_UNIT,*) '  tlow, thi = ',          &
                                    topo_data(i)%t_low, topo_data(i)%t_high

                    if (abs(topo_data(i)%topo_type) == 1) then
                        print *, 'WARNING: topotype 1 has been deprecated'
                        print *, 'converting to topotype > 1 is encouraged'
                        print *, 'python tools for converting files are provided'
                    end if

                    ! Read in header data
                    call read_topo_header(topo_data(i))
                end do


                ! Check that topo arrays cover full domain:
                ! TODO:  Need to fix this
                call topoarea(xlower,xupper,ylower,yupper,1,area)
                area_domain = (yupper-ylower)*(xupper-xlower)
                if (abs(area - area_domain) > 1d-2 * area_domain) then
                    print *, '**** topo arrays do not cover domain'
                    print *, '**** area of overlap = ', area
                    print *, '**** area of domain  = ', area_domain
                    stop
                else if (abs(area - area_domain) > 1d-12 * area_domain) then
                    print *, '**** WARNING'
                    print *, '**** topo arrays do not quite cover domain'
                    print *, '**** area of overlap = ', area
                    print *, '**** area of domain  = ', area_domain
                    print *, '**** error is less than 1% so proceeding...'
                endif

                ! Fill in data for dtopo files - basically copied from the dtopo
                ! module data
                do i = num_topo_files + 1, num_dtopo_files
                    topo_data(i)%topo_type = dtopo_data%dtopo_type
                    topo_data(i)%num_cells = dtopo_data%num_cells
                    topo_data(i)%x_lower = dtopo_data%x_lower
                    topo_data(i)%y_lower = dtopo_data%y_lower
                    topo_data(i)%x_upper = dtopo_data%x_upper
                    topo_data(i)%y_upper = dtopo_data%y_upper
                    topo_data(i)%dx = dtopo_data%dx
                    topo_data(i)%dy = dtopo_data%dy
                    topo_data(i)%min_level = dtopo_data%min_level
                    topo_data(i)%max_level = dtopo_data%max_level
                    topo_data(i)%t_lower = dtopo_data%t_lower
                    topo_data(i)%t_upper = dtopo_data%t_upper

                    allocate(topo_data(i)%topo(topo_data(i)%num_cells(1),   &
                                               topo_data(i)%num_cells(2)))
                end do

                ! Read topography and allocate space for each file
                do i = 1, num_topo_files
                    call read_topo_file(topo_data(i))
                    
                    ! At this point see if this topography file intersects with
                    ! a dtopo file.  If so mark it as needing saving
                    do j = num_topo_files + 1, num_dtopo_files
                        if ((topo_files(i)%x_lower < topo_files(j)%x_lower) .or. &
                            (topo_files(i)%y_lower < topo_files(j)%y_lower) .or. &
                            (topo_files(i)%x_higher > topo_files(j)%x_higher) .or. &
                            (topo_files(i)%y_higher > topo_files(j)%y_higher) ) then

                            topo_files(i)%save_topo = .true.
                        end if
                    end do
                end do

                ! Rank topography in order of resolution - does not include
                ! work arrays from dtopo
                do i = 1, num_topo_files
                    finer_than = 0
                    do j = 1, num_topo_files
                        if (j /= i) then
                            area_i = topo_data(i)%dx * topo_data(i)%dy
                            area_j = topo_data(i)%dx * topo_data(i)%dy

                            ! if two files have the same resolution, order is
                            ! arbitrarily chosen
                            if ((area_i < area_j) .or.  &
                                (area_i == area_j).and.(j < i)) then
                                finer_than = finer_than + 1
                            end if
                        endif
                    enddo
                    rank = num_topo_files - finer_than
                    topo_order(rank) = i
                enddo

                ! This only displays the ranks of actual topography files
                write(GEO_PARM_UNIT,*) ' '
                write(GEO_PARM_UNIT,*) '  Ranking of topography files', &
                    '  finest to coarsest: ', &
                    (topo_order(rank), rank = 1, num_topo_files)
                write(GEO_PARM_UNIT,*) ' '

                ! Set the values in the topo work arrays that deal with dtopo
                do i = num_topo_files + 1, num_dtopo_files
                    call set_topo_for_dtopo(topo_data(i))
                end do

                ! Allocate extra topo data work space for the topography files
                ! that need it
                topo_finalized = .true.
                aux_finalized = 2
                if (num_dtopo_files > 0) then

                    topo_finalized = .false.
                    aux_finalized = 0



                end if

                !create topo0work array for finest arrays covering dtopo
                !arrays to be saved are indicated in topo0save
                topo_finalized = .true.
                aux_finalized = 2   !# indicates aux arrays properly set with dtopo
                if (num_dtopo>0) then
                   topo_finalized = .false.
                   aux_finalized = 0  !# will be incremented each time level 1 goes
                   i0topo0(1) = 1
                   mtopo0size = dot_product(mtopo,topo0save)
                   allocate(topo0work(mtopo0size))
                   do i = 2,mtopofiles
                      i0topo0(i)= i0topo0(i-1) + mtopo(i-1)*topo0save(i-1)
                   enddo

                   do i = 1,mtopofiles
                      if (topo0save(i)>0) then
                         topo0work(i0topo0(i):i0topo0(i)+mtopo(i)-1) = &
                            topowork(i0topo(i):i0topo(i)+mtopo(i)-1)
                      endif
                   enddo
                endif

            !-------------- tests for synthetic bathymetry ------------------
            ! Simple jump discontinuity in bathymetry
            else if (test_topography == 1) then
                topo_finalized = .true.
                read(iunit,"(d16.8)") topo_location
                read(iunit,"(d16.8)") topo_left
                read(iunit,"(d16.8)") topo_right

            ! Idealized ocean shelf
            else if (test_topography == 2 .or. test_topography == 3) then
                topo_finalized = .true.
                read(iunit,"(d16.8)") topo_x0
                read(iunit,"(d16.8)") topo_x1
                read(iunit,"(d16.8)") topo_x2
                read(iunit,"(d16.8)") topo_basin_depth
                read(iunit,"(d16.8)") topo_shelf_depth
                read(iunit,"(d16.8)") topo_beach_slope
                topo_shelf_slope = (topo_basin_depth - topo_shelf_depth) &
                                            / (topo_x0 - topo_x1)
            else
                print *,"Error:  Unknown test topography type ",test_topography
                stop
            endif

            module_setup = .true.
        end if

    end subroutine read_topo_settings

    ! ========================================================================
    !  read_topo_file(mx,my,topo_type,fname,xll,yll,topo)
    !
    !  Read topo file.
    ! ========================================================================

    subroutine read_topo_file(topo_data)

#ifdef NETCDF
        use netcdf
#endif

        ! use geoclaw_module, only:
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Arguments
        type(topo_file_type), intent(inout) :: topo_data

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,topo_start,n
        real(kind=8) :: no_data_value,x,y,z,topo_temp
        real(kind=8) :: values(10)
        character(len=80) :: str
        integer(kind=4) :: row_index

        ! NetCDF Support
        character(len=10) :: direction, x_dim_name, x_var_name, y_dim_name, &
            y_var_name, z_var_name, var_name, dim_name_tmp
        real(kind=8), allocatable :: nc_buffer(:, :), xlocs(:), ylocs(:)
        integer(kind=4) :: x_var_id, y_var_id, z_var_id, x_dim_id, y_dim_id
        integer(kind=4) :: xstart(1), ystart(1), mx_tot, my_tot, m_tmp
        integer(kind=4) :: ios, nc_file, num_values, dim_ids(2), num_dims, &
            var_type, var_ids(2), num_vars, num_dims_tot, z_dim_ids(2)

        ! Format statements


        print *, ' '
        print *, 'Reading topography file  ', topo_data%file_path

        select case(abs(topo_type))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                i = 0
                status = 0
                do i=1,mx*my
                    read(iunit,fmt=*,iostat=status) x,y,topo_temp
                    if ((i > mx * my) .and. (status == 0)) then
                        print *,'*** Error: i > mx*my = ',mx*my
                        print *,'*** i, mx, my: ',i,mx,my
                        print *,'*** status = ',status
                        stop
                    endif

                    if (status /= 0) then
                        print *,"Error reading topography file, reached EOF."
                        print *,"  File = ",fname
                        stop
                    else
                        ! TODO: Translate to actual 2D array
                        topo_data%z = topo_temp
                    endif
                enddo
                
                close(unit=iunit)

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if topo_type=2 or
            ! mx values per line if topo_type=3
            ! ================================================================
            case(2:3)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo
                
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                no_data_value = values(1)

                ! Read in data
                missing = 0
                select case(abs(topo_type))
                    ! TODO: modify this to also be 2D always
                    case(2)
                        do i=1,mx*my
                            read(iunit,*) topo_data%z(i)
                            if (topo(i) == no_data_value) then
                                missing = missing + 1
                                topo_data%z(i) = topo_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (topo((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (topo((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    topo((j-1)*mx + i) = topo_missing
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print ('WARNING... ', i6, ' missing data values in this ', &
                           'topofile.'), missing
                    print ('   These values have been set to topo_missing = ', &
                          f13.3, ' in read_topo_file.'), topo_missing
                    if (topo_missing > 9999.d0) then
                        print *, 'ERROR... do not use this default value'
                        print *, 'Fix your topofile or set'
                        print *, '  rundata.topo_data.topo_missing in setrun.py'
                        stop
                    endif
                    print *, ' '
                endif

                close(unit=iunit)
            
            ! NetCDF
            case(4)
#ifdef NETCDF
                ! Open file    
                call check_netcdf_error(nf90_open(fname, nf90_nowrite, nc_file))
                
                ! Get number of dimensions and vars
                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, &
                    num_vars))
                
                ! get x and y dimension info
                call get_dim_info(nc_file, num_dims_tot, x_dim_id, x_dim_name, &
                    mx_tot, y_dim_id, y_dim_name, my_tot)
                    
                allocate(xlocs(mx_tot),ylocs(my_tot))
                
                call check_netcdf_error(nf90_get_var(nc_file, x_dim_id, xlocs, start=(/ 1 /), count=(/ mx_tot /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_dim_id, ylocs, start=(/ 1 /), count=(/ my_tot /)))
                xstart = minloc(xlocs, mask=(xlocs.eq.xll))
                ystart = minloc(ylocs, mask=(ylocs.eq.yll))
                deallocate(xlocs,ylocs)
                
                z_var_id = -1
                do n=1, num_vars
                    call check_netcdf_error(nf90_inquire_variable(nc_file, n, var_name, var_type, num_dims, dim_ids))
                    
                    ! Assume dim names are same as var ids
                    if (var_name == x_dim_name) then
                        x_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, x_var_name, x_var_id))
                    else if (var_name == y_dim_name) then
                        y_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, y_var_name, y_var_id))
                    else if (num_dims == 2) then
                        z_var_name = var_name
                        z_dim_ids = dim_ids
                        call check_netcdf_error(nf90_inq_varid(nc_file, z_var_name, z_var_id))
                    end if

                end do
                if (z_var_id == -1) then
                    stop "Unable to find topography data!"
                end if

                ! Load in data
                ! TODO: Provide striding into data if need be
                
                ! only try to access data if theres overlap with the domain
                if ((mx > 0) .and. (my > 0)) then
                    if (z_dim_ids(1) == x_dim_id) then
                        allocate(nc_buffer(mx, my))
                    else if (z_dim_ids(1) == y_dim_id) then
                        allocate(nc_buffer(my, mx))
                    else 
                        stop " NetCDF z variable has dimensions that don't align with x and y"
                    end if

                    call check_netcdf_error(nf90_get_var(nc_file, z_var_id, nc_buffer, &
                                                start = (/ xstart(1), ystart(1) /), &
                                                count = (/ mx, my /)))

                    ! check for lon/lat ordering of z variable
                    ! TODO: Modify to remain 2D
                    if (z_dim_ids(1) == x_dim_id) then
                        do j = 0, my - 1
                            topo(j * mx + 1:(j + 1) * mx) = nc_buffer(:, my - j)
                        end do
                    else if (z_dim_ids(1) == y_dim_id) then
                        do j = 0, my - 1
                            topo(j * mx + 1:(j + 1) * mx) = nc_buffer(my - j, :)
                        end do
                    end if
                    deallocate(nc_buffer)

                    ! Check if the topography was defined positive down and flip the
                    ! sign if need be.  Just in case this is true but topo_type < 0
                    ! we do not do anything here on this to avoid doing it twice.
                    ios = nf90_get_att(nc_file, z_var_id, 'positive', direction)
                    call check_netcdf_error(nf90_close(nc_file))
                    if (ios == NF90_NOERR) then
                        if (to_lower(direction) == "down") then
                            if (topo_type < 0) then
                                topo = -topo
                            endif
                        end if
                    end if
                end if
#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of GeoClaw."
                stop
#endif
        end select

        ! Handle negative topo types
        if (topo_data%topo_type < 0) then
            topo_data%z = -topo_data%z
        endif

        ! ====================================================================

    end subroutine read_topo_file

    ! ========================================================================
    ! subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read topo file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - topo_type - (int) Type of topography file (-3 < topo_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_topo_header(topo_data)
#ifdef NETCDF
        use netcdf
#endif

        ! use geoclaw_module
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Input and Output
        type(topo_file_type), intent(inout) :: topo_data

        ! Local
        integer, parameter :: iunit = 19
        integer :: topo_size, status, n, i
        real(kind=8) :: x,y,z,nodata_value
        logical :: found_file
        real(kind=8) :: values(10)
        character(len=80) :: str
        logical :: verbose
        logical :: xll_registered, yll_registered

        ! NetCDF Support
        ! character(len=1) :: axis_string
        ! character(len=6) :: convention_string
        ! integer(kind=4) :: convention_version
        integer(kind=4) :: ios, nc_file, num_values
        real(kind=8), allocatable :: xlocs(:),ylocs(:)
        logical, allocatable :: x_in_dom(:),y_in_dom(:)
        integer(kind=4) :: dim_ids(2), num_dims, var_type, var_ids(2), num_vars, num_dims_tot
        character(len=10) :: var_name, x_var_name, y_var_name, z_var_name
        character(len=10) :: x_dim_name, y_dim_name
        integer(kind=4) :: x_var_id, y_var_id, z_var_id, x_dim_id, y_dim_id

        verbose = .false.

        inquire(file=topo_data%file_path, exist=found_file)
        if (.not. found_file) then
            print *, 'Missing topography file:'
            print *, '   ', topo_data%file_path
            stop
        endif

        select case(abs(topo_data%topo_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                open(unit=iunit, file=fname, status='unknown',form='formatted')

                ! Initial size variables
                topo_size = 0
                topo_data%num_cells(1) = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) topo_data%x_lower, topo_data%y_upper
                topo_size = topo_size + 1
                topo_data%num_cells(1) = topo_data%num_cells(1) + 1

                ! Go through first row figuring out mx, continue to count
                y = topo_data%y_upper
                do while (topo_data%y_upper == y)
                    read(iunit,*) x, y, z
                    topo_size = topo_size + 1
                    topo_data%num_cells(1) = topo_data%num_cells(1) + 1
                enddo
                topo_data%num_cells(1) = topo_data%num_cells(1) - 1
                ! Continue to count the rest of the lines
                do
                    read(iunit,fmt=*,iostat=status) x, y, z
                    if (status /= 0) exit
                    topo_size = topo_size + 1
                enddo
                if (status > 0) then
                    print *,"ERROR:  Error reading header of topography file ", topo_data%file_path
                    stop
                endif

                ! Calculate remaining values
                topo_data%num_cells(2) = topo_size / topo_data%num_cells(1)
                opo_data%x_upper = x
                topo_data%y_lower = y
                topo_data%dx = (opo_data%x_upper-topo_data%x_lower) / (topo_data%num_cells(1)-1)
                topo_data%dy = (topo_data%y_upper-topo_data%y_lower) / (topo_data%num_cells(2)-1)
                
            ! ASCII file with header followed by z data
            case(2:3)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                topo_data%num_cells(1) = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                topo_data%num_cells(2) = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                topo_data%x_lower = values(1)
                str = to_lower(str)  ! convert to lower case
                xll_registered = (index(str, 'xllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                topo_data%y_lower = values(1)
                str = to_lower(str)  ! convert to lower case
                yll_registered = (index(str, 'yllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                topo_data%dx = values(1)
                if (n == 2) then
                    topo_data%dy = values(2)
                else
                    topo_data%dy = topo_data%dx
                endif

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                nodata_value = values(1)

                if (xll_registered) then
                    topo_data%x_lower = topo_data%x_lower + 0.5d0 * topo_data%dx
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting xllcorner by 0.5*dx to cell center'
                endif 

                if (yll_registered) then
                    topo_data%y_lower = topo_data%y_lower + 0.5d0 * topo_data%dy
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting yllcorner by 0.5*dy to cell center'
                endif 


                topo_data%x_upper = topo_data%x_lower +         &
                                (topo_data%num_cells(1) - 1) * topo_data%dx
                topo_data%y_upper = topo_data%y_lower +         &
                                (topo_data%num_cells(2) - 1) * topo_data%dy
                
            ! NetCDF
            case(4)
#ifdef NETCDF

                ! Open file    
                call check_netcdf_error(nf90_open(fname, nf90_nowrite, nc_file))

                ! get x and y dimension info
                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, &
                    num_vars))
                call get_dim_info(nc_file, num_dims_tot, x_dim_id, x_dim_name, &
                    mx, y_dim_id, y_dim_name, my)
                
                ! allocate vector to hold lon and lat vals
                allocate(xlocs(mx),ylocs(my),x_in_dom(mx),y_in_dom(my))

                if (verbose) then
                    print *, "Names = (", x_dim_name, ", ", y_dim_name, ")"
                    print *, "Size = (", mx, ", ", my, ")"
                end if

                ! Read in variables
                call check_netcdf_error(nf90_inq_varids(nc_file, num_vars, var_ids))
                if (verbose) then
                    print *, "n, var_name, var_type, num_dims, dim_ids"
                end if

                do n=1, num_vars
                    call check_netcdf_error(nf90_inquire_variable(nc_file, n, var_name, var_type, num_dims, dim_ids))
                    if (verbose) then
                        print *, n, ": ", var_name, "|", var_type, "|", num_dims, "|", dim_ids
                    end if

                    ! Assume dim names are same as var ids
                    if (var_name == x_dim_name) then
                        x_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, x_var_name, x_var_id))
                        ! x_var_id = n
                        ! x_var_name = var_name
                    else if (var_name == y_dim_name) then
                        y_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, y_var_name, y_var_id))
                        ! y_var_id = n
                        ! y_var_name = var_name
                    ! Assume that this is the topography data
                    else if (num_dims == 2) then
                        z_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, z_var_name, z_var_id))
                        ! z_var_id = n
                        ! z_var_name = var_name
                    else
                        if (verbose) then
                            print *, "Not using var_id ", n
                        end if
                    end if

                end do

                if (verbose) then
                    print *, "x_var_name, x_var_id = ", x_var_name, x_var_id
                    print *, "y_var_name, y_var_id = ", y_var_name, y_var_id
                    print *, "z_var_name, z_var_id = ", z_var_name, z_var_id
                end if

                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xlocs, start=(/ 1 /), count=(/ mx /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, ylocs, start=(/ 1 /), count=(/ my /)))
                
                topo_data%dx = xlocs(2) - xlocs(1)
                topo_data%dy = ylocs(2) - ylocs(1)
                
                ! find which locs are within domain (with a dx/dy buffer around domain)
                x_in_dom = (xlocs.gt.(xlower-dx)) .and. (xlocs.lt.(xupper+dx))
                y_in_dom = (ylocs.gt.(ylower-dy)) .and. (ylocs.lt.(yupper+dy))
                
                topo_data%x_lower = minval(xlocs, mask=x_in_dom)
                topo_data%y_lower = minval(ylocs, mask=y_in_dom)
                topo_data%x_upper = maxval(xlocs, mask=x_in_dom)
                topo_data%y_upper = maxval(ylocs, mask=y_in_dom)
                
                ! adjust mx, my
                topo_data%num_cells(1) = count(x_in_dom)
                topo_data%num_cells(2) = count(y_in_dom)

                call check_netcdf_error(nf90_close(nc_file))
                deallocate(xlocs,ylocs,x_in_dom,y_in_dom)
#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of GeoClaw."
                stop
#endif

            case default
                print *, 'ERROR:  Unrecognized topo_type'
                print *, '    topo_type = ', topo_data%topo_type
                print *, '  for topography file:'
                print *, '   ', topo_data%file_path
                stop
        end select

        close(iunit)
        write(GEO_PARM_UNIT,*) '  mx = ',topo_data%num_cells(1),'  x = (',   &
                                    topo_data%x_lower,',',topo_data%x_upper,')'
        write(GEO_PARM_UNIT,*) '  my = ',topo_data%num_cells(2),'  y = (',   &
                                    topo_data%y_lower,',',topo_data%y_upper,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', &
                                    topo_data%dx, topo_data%dy

    end subroutine read_topo_header

        ! ========================================================================
    !  set_topo_for_dtopo()
    !
    !  Set topography values in new topo arrays that correspond to dtopo spatialy
    !  array values come from the finest topography already in topowork
    ! ========================================================================

    subroutine set_topo_for_dtopo(mx,my,dx,dy,xlow,ylow,xhi,yhi,newtopo)

        !arguments
        integer, intent(in) :: mx,my
        real(kind=8), intent(in) :: dx,dy,xlow,xhi,ylow,yhi
        real(kind=8), intent(inout) :: newtopo(1:mx*my)

        !locals
        integer :: i,j,k,ij,id,irank,itopo1,itopo2,jtopo1,jtopo2
        integer :: ijll,ijlr,ijul,ijur
        real(kind=8) :: x,y,xl,xr,yu,yl,zll,zlr,zul,zur,z,dxdy

        do j=1,my
               y = yhi - (j-1)*dy
            do i=1,mx
               x = xlow + (i-1)*dx
               ij = (j-1)*mx + i
               !find intersection starting from finest topo
               !all points must lie in some topo file therefore the
               !finest topo file for all dtopo points will be saved in topo0
               do irank = 1,mtopofiles
                  id = mtopoorder(irank)
                  if (id.gt.mtopofiles-num_dtopo) then
                     !this is another dtopo ==> topo file: skip
                     cycle
                  elseif ( (x>xhitopo(id)).or.(x<xlowtopo(id)).or. &
                          (y>yhitopo(id)).or.(y<ylowtopo(id))) then
                     !no intersection
                     cycle
                  else !lies in this topofile

                     ! Old way of setting topo0save up to v5.4.1.
                     ! This assumed topo file did not
                     ! intersect dtopo file if this point was never reached.
                     ! Not true if topofile has such small extent that it lies between
                     ! dtopo points.  Now instead we set topo0save earlier.
                     !topo0save(id) = 1

                     !find indices for bilinear cell in topo
                     !arrays are in form of DEM...high y values first
                     !note for xy points lying on nodes all indices will be equal
                     itopo1 = int(floor((x-xlowtopo(id))/dxtopo(id)))+1
                     itopo2 = int(ceiling((x-xlowtopo(id))/dxtopo(id)))+1
                     jtopo1 = int(floor((yhitopo(id)-y)/dytopo(id))) + 1
                     jtopo2 = int(ceiling((yhitopo(id)-y)/dytopo(id))) + 1
                     !indices for work array
                     ijll = i0topo(id) + (jtopo2-1)*mxtopo(id) + itopo1 -1
                     ijlr = i0topo(id) + (jtopo2-1)*mxtopo(id) + itopo2 -1
                     ijul = i0topo(id) + (jtopo1-1)*mxtopo(id) + itopo1 -1
                     ijur = i0topo(id) + (jtopo1-1)*mxtopo(id) + itopo2 -1
                     !find x,y,z values for bilinear
                     !z may be from only 1 or 2 nodes for aligned grids
                     !bilinear should still evaluate correctly
                     zll = topowork(ijll)
                     zlr = topowork(ijlr)
                     zul = topowork(ijul)
                     zur = topowork(ijur)
                     xl = xlowtopo(id) + real(itopo1-1,kind=8)*dxtopo(id)
                     xr = xl + dxtopo(id)
                     yu = yhitopo(id) - real(jtopo1-1,kind=8)*dytopo(id)
                     yl = yu - dytopo(id)
                     dxdy = dxtopo(id)*dytopo(id)
                     z = zll*(xr-x)*(yu-y) + zlr*(x-xl)*(yu-y) + zul*(xr-x)*(y-yl) + zur*(x-xl)*(y-yl)
                     newtopo(ij) = z/dxdy
                     !this was the finest topo file, move to next point
                     exit
                  endif
               enddo
            enddo
        enddo

    end subroutine set_topo_for_dtopo

    real(kind=8) pure function test_topo(x,y) result(topography)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x,y

        if (test_topography == 1) then
            if (x < topo_location) then
                topography = topo_left
            else
                topography = topo_right
            endif
        else if (test_topography == 2) then
            if (x < topo_x0) then
                topography = topo_basin_depth
            else if (topo_x0 <= x .and. x < topo_x1) then
                topography = topo_shelf_slope * (x-topo_x0) + topo_basin_depth
            else if (topo_x1 <= x .and. x < topo_x2) then
                topography = topo_shelf_depth
            else
                topography = topo_beach_slope * (x-topo_x2) + topo_shelf_depth
            endif
        endif

    end function test_topo


    recursive subroutine topoarea(x1, x2, y1, y2, m, area)

        ! Compute the area of overlap of topo with the rectangle (x1,x2) x (y1,y2)
        ! using topo arrays indexed mtopoorder(mtopofiles) through mtopoorder(m) 
        ! (coarse to fine).

        ! The main call to this subroutine has corners of a physical domain for
        ! the rectangle and m = 1 in order to compute the area of overlap of
        ! domain by all topo arrays.  Used to check inputs and insure topo
        ! covers domain.

        ! The recursive strategy is to first compute the area using only topo 
        ! arrays mtopoorder(mtopofiles) to mtopoorder(m+1), 
        ! and then apply corrections due to adding topo array mtopoorder(m).
         
        ! Corrections are needed if the new topo array intersects the grid cell.
        ! Let the intersection be (x1m,x2m) x (y1m,y2m).
        ! Two corrections are needed, first to subtract out the area over
        ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
        ! topo arrays mtopoorder(mtopofiles) to mtopoorder(m+1),
        ! and then adding in the area over this same region using 
        ! topo array mtopoorder(m).

        ! Based on the recursive routine rectintegral that integrates
        ! topo over grid cells using a similar strategy.

        implicit none

        ! arguments
        real (kind=8), intent(in) :: x1,x2,y1,y2
        integer, intent(in) :: m
        real (kind=8), intent(out) :: area

        ! local
        real(kind=8) :: xmlo,xmhi,ymlo,ymhi,x1m,x2m, &
            y1m,y2m, area1,area2,area_m
        integer :: mfid, indicator, i0
        real(kind=8), external :: topointegral  


        mfid = mtopoorder(m)
        i0=i0topo(mfid)

        if (m == mtopofiles) then
             ! innermost step of recursion reaches this point.
             ! only using coarsest topo grid -- compute directly...
             call intersection(indicator,area,xmlo,xmhi, &
                 ymlo,ymhi, x1,x2,y1,y2, &
                 xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

        else
            ! recursive call to compute area using one fewer topo grids:
            call topoarea(x1,x2,y1,y2,m+1,area1)

            ! region of intersection of cell with new topo grid:
            call intersection(indicator,area_m,x1m,x2m, &
                 y1m,y2m, x1,x2,y1,y2, &
                 xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

            
            if (area_m > 0) then
            
                ! correction to subtract out from previous set of topo grids:
                call topoarea(x1m,x2m,y1m,y2m,m+1,area2)
        
                ! adjust integral due to corrections for new topo grid:
                area = area1 - area2 + area_m
            else
                area = area1
            endif
        endif

    end subroutine topoarea

    recursive subroutine rectintegral(x1,x2,y1,y2,m,integral)

        ! Compute the integral of topo over the rectangle (x1,x2) x (y1,y2)
        ! using topo arrays indexed mtopoorder(mtopofiles) through mtopoorder(m) 
        ! (coarse to fine).

        ! The main call to this subroutine has corners of a grid cell for the 
        ! rectangle and m = 1 in order to compute the integral over the cell 
        ! using all topo arrays.

        ! The recursive strategy is to first compute the integral using only topo 
        ! arrays mtopoorder(mtopofiles) to mtopoorder(m+1), 
        ! and then apply corrections due to adding topo array mtopoorder(m).
         
        ! Corrections are needed if the new topo array intersects the grid cell.
        ! Let the intersection be (x1m,x2m) x (y1m,y2m).
        ! Two corrections are needed, first to subtract out the integral over
        ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
        ! topo arrays mtopoorder(mtopofiles) to mtopoorder(m+1),
        ! and then adding in the integral over this same region using 
        ! topo array mtopoorder(m).

        ! Note that the function topointegral returns the integral over the 
        ! rectangle based on a single topo array, and that routine calls
        ! bilinearintegral.


        implicit none

        ! arguments
        real (kind=8), intent(in) :: x1,x2,y1,y2
        integer, intent(in) :: m
        real (kind=8), intent(out) :: integral

        ! local
        real(kind=8) :: xmlo,xmhi,ymlo,ymhi,area,x1m,x2m, &
            y1m,y2m, int1,int2,int3
        integer :: mfid, indicator, mp1fid, i0
        real(kind=8), external :: topointegral  


        mfid = mtopoorder(m)
        i0=i0topo(mfid)

        if (m == mtopofiles) then
             ! innermost step of recursion reaches this point.
             ! only using coarsest topo grid -- compute directly...
             call intersection(indicator,area,xmlo,xmhi, &
                 ymlo,ymhi, x1,x2,y1,y2, &
                 xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

             if (indicator.eq.1) then
                ! cell overlaps the file
                ! integrate surface over intersection of grid and cell
                integral = topointegral( xmlo,xmhi,ymlo, &
                        ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                        dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
             else
                integral = 0.d0
             endif

        else
            ! recursive call to compute area using one fewer topo grids:
            call rectintegral(x1,x2,y1,y2,m+1,int1)

            ! region of intersection of cell with new topo grid:
            call intersection(indicator,area,x1m,x2m, &
                 y1m,y2m, x1,x2,y1,y2, &
                 xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

            
            if (area > 0) then
            
                ! correction to subtract out from previous set of topo grids:
                call rectintegral(x1m,x2m,y1m,y2m,m+1,int2)
        
                ! correction to add in for new topo grid:
                int3 = topointegral(x1m,x2m, y1m,y2m, &
                            xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                            dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
        
                ! adjust integral due to corrections for new topo grid:
                integral = int1 - int2 + int3
            else
                integral = int1
            endif
        endif

    end subroutine rectintegral

    subroutine intersection(indicator,area,xintlo,xinthi, &
               yintlo,yinthi,x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi)

        ! find the intersection of two rectangles, return the intersection
        ! and it's area, and indicator =1
        ! if there is no intersection, indicator =0

          implicit none

          integer, intent(out) :: indicator

          real(kind=8), intent(in) ::  x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi
          real(kind=8), intent(out) :: area,xintlo,xinthi,yintlo,yinthi

          xintlo=dmax1(x1lo,x2lo)
          xinthi=dmin1(x1hi,x2hi)
          yintlo=dmax1(y1lo,y2lo)
          yinthi=dmin1(y1hi,y2hi)


          if (xinthi.gt.xintlo.and.yinthi.gt.yintlo) then
             area = (xinthi-xintlo)*(yinthi-yintlo)
             indicator = 1
          else
             area = 0.d0
             indicator = 0
          endif

    end subroutine intersection


#ifdef NETCDF
    subroutine check_netcdf_error(ios)

        use netcdf

        implicit none

        integer, intent(in) :: ios

        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if

    end subroutine check_netcdf_error

    subroutine get_dim_info(nc_file, ndims, x_dim_id, x_dim_name, mx, &
        y_dim_id, y_dim_name, my)
        use netcdf
        implicit none
        integer, intent(in) :: nc_file, ndims
        integer, intent(out) :: x_dim_id, y_dim_id, mx, my
        character (len = *), intent(out) :: x_dim_name, y_dim_name
        integer :: m_tmp, n
        character(20) :: dim_name_tmp

        ! get indices to start at for reading netcdf within domain
        do n=1, ndims
            call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                n, dim_name_tmp, m_tmp))
            if (ANY((/ 'LON      ','LONGITUDE','X        ' /) == Upper(dim_name_tmp))) then
                x_dim_name = dim_name_tmp
                mx = m_tmp
                x_dim_id = n
            else if (ANY((/ 'LAT     ','LATITUDE','Y       ' /) == Upper(dim_name_tmp))) then
                y_dim_name = dim_name_tmp
                my = m_tmp
                y_dim_id = n
            end if
        end do
    end subroutine get_dim_info
#endif

    function Upper(s1)  RESULT (s2)
        CHARACTER(*)       :: s1
        CHARACTER(LEN(s1)) :: s2
        CHARACTER          :: ch
        INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
        INTEGER            :: i

        DO i = 1,LEN(s1)
           ch = s1(i:i)
           IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
           s2(i:i) = ch
        END DO
    end function Upper

end module topo_module