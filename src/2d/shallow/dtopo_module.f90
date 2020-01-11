! ==============================================================================
!  Variable Topography module
!
!    This module contains storage and routines for dealing with variable 
!    topography in GeoClaw.
!    
! ==============================================================================

module dtopo_module

    implicit none

    save
    
    logical, private :: module_setup = .false.

    ! Moving topography
    type dtopo_file_type

        ! dTopography data
        integer :: dtopo_type
        real(kind=8), allocatable :: dz(:, :, :)

        ! Geometry in space-time
        integer :: num_cells(3)
        real(kind=8) :: x_lower, x_upper, y_lower, y_upper, t_lower, t_upper
        real(kind=8) :: dx, dy, dt

        ! Level information
        integer :: min_level, max_level

        ! Originating file
        character(len=256) :: file_path

    end type dtopo_file_type

    integer, private :: num_dtopo_files = 0
    integer, private, allocatable :: dtopo_order(:)
    type(dtopo_file_type), pointer, private :: dtopo_data(:)
    real(kind=8) :: max_dt_dtopo

contains

    ! ========================================================================
    !  set_dtopo(fname)
    ! ========================================================================
    ! Read moving topography info from setdtopo.data
    ! Time-dependend topography is used to initiate tsunami, for example.
    !
    ! If num_dtopo = 0, no movement of topography specified.
    !
    ! If num_dtopo = 1, the topo changes dynamically.
    ! the topofile can then be formated in the following ways
    ! dtopotype > 1: file contains a header, with a line for each record
    ! my; mx; mt; xlower; ylower; t0; dx; dy; dt;
    ! dtopotype = 1:
    ! Then the dtopofile should have 4 columns:
    ! time,longitude,latitude,vertical displacement(m) since t0
    !
    ! Longitude and latitude advance in the standard GIS way from
    ! upper left corner across in x and then down in y.
    ! Time column advances most slowly.
    ! ========================================================================
    subroutine read_dtopo_settings(file_name)

        use geoclaw_module, only: GEO_PARM_UNIT

        implicit none

        ! Input arguments
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: iunit = 79
        integer :: itopo,finer_than,rank
        real(kind=8) :: area_i,area_j
        real(kind=8) :: xcell, xim, xip, ycell, yjm, yjp, ztopoij
        real(kind=8) :: capac_area
        integer :: i,j,m,ib,jb,ij,ijdtopo,jbr

        if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SETDTOPO:'
            write(GEO_PARM_UNIT,*) '-------------'

            if (present(file_name)) then
                call opendatafile(iunit,file_name)
            else
                call opendatafile(iunit,'dtopo.data')
            endif

            read(iunit,*) num_dtopo_files
            write(GEO_PARM_UNIT,*) '   num dtopo files = ', num_dtopo_files
            if (num_dtopo_files == 0) then
                return
            endif

            ! Allocate and read in dtopo info
            allocate(dtopo_data(num_dtopo_files))
            allocate(dtopo_order(num_dtopo_files))

            do i = 1, num_dtopo_files
                read(iunit,*) dtopo_data(i)%file_path
                read(iunit,*) dtopo_data(i)%dtopo_type,        &
                read(inunt,*) dtopo_data(i)%minleveldtopo,     &
                read(iunit,*) dtopo_data(i)%maxleveldtopo

                allocate(dtopo_data(i)%dz(:, :, :))

                write(GEO_PARM_UNIT,*) '   fname:', dtopo_data(i)%file_path
                write(GEO_PARM_UNIT,*) '   topo type:',dtopo_data(i)%dtopo_type
                write(GEO_PARM_UNIT,*) '   minlevel, maxlevel:'
                write(GEO_PARM_UNIT,*)       dtopo_data(i)%min_level,   &
                                             dtopo_data(i)%max_level

                ! Read in header data
                call read_dtopo_header(dtopo_data(i))
            enddo

            ! Largest allowable dt while dtopo is moving
            read(iunit,*) max_dt_dtopo
            
            ! dtopo order..for updating topo from finest dtopo model
            ! The finest topography will be given priority in any region
            ! mtopoorder(rank) = i means that i'th topography file has rank rank,
            ! where the file with rank=1 is the finest and considered first.
            do i = 1, num_dtopo_files
                finer_than = 0
                do j = 1, num_dtopo_files
                    if (j /= i) then
                        area_i = dtopo_data(i)%dx * dtopo_data(i)%dy
                        area_j = dtopo_data(j)%dx * dtopo_data(j)%dy

                        ! If two files have the same resolution, order is
                        ! arbitrarily chosen
                        if ((area_i < area_j) .or.          &
                            ((area_i == area_j) .and. (j < i)) ) then
                            
                            finer_than = finer_than + 1
                        end if
                    endif
                enddo
                ! ifinerthan tells how many other files, file i is finer than
                rank = num_dtopo_files - finer_than
                dtopo_order(rank) = i
            enddo

            ! Read in dtopo data
            do i = 1, num_dtopo_files
                call read_dtopo(dtopo_data(i))
            enddo

            module_setup = .true.
        end if

    end subroutine read_dtopo_settings
    ! ========================================================================

    ! ========================================================================
    !  subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi,
    !                               yhi,tf,dx,dy,dt)
    ! ========================================================================

    !  Read in dtopo file header and either read or calculate the grid info
    !
    !  :Input:
    !   - fname - (char) Name of the dtopo file
    !   - topo_type - (int) Topography file type (1-3 are valid)
    !
    !  :Output:
    !   - mx,my,mt - (int) Number of grid point in space (mx,my) and time (mt)
    !   - xlow,ylow - (dp) Lower corner spatial coordinate of grid
    !   - xhi,yhi - (dp) Upper corner spatial coodinate of grid
    !   - t0,tf - (dp) Beginning and end times for the file
    !   - dx,dy,dt - (dp) Distance between space (dx,dy) and time (dt) points
    ! ========================================================================
    subroutine read_dtopo_header(dtopo_data)

        implicit none

        ! Locals
        integer, parameter :: iunit = 7
        integer :: status
        real(kind=8) :: x,y,t,y_old,t_old
        logical :: found_file

        ! Open file
        inquire(file=dtopo_data%file_path, exist=found_file)
        if (.not.found_file) then
            print *, 'Missing dtopo file:'
            print *, '    ', dtopo_data%file_path
            stop
        endif
        open(unit=iunit,file=fname,status='unknown',form='formatted')

        select case(topo_type)
            ! Old style ASCII dtopo files
            case(1)
                ! Initial size variables
                topo_size = 0
                dtopo_data%num_cells(2) = 1
                dtopo_data%num_cells(3) = 1

                ! Read in first values, determines xlow, yhi and t0
                read(iunit,*) t0, dtopo_data%x_lower, dtopo_data%y_upper
                topo_size = topo_size + 1
                t = t0
                y_old = dtopo_data%y_upper
                ! Go through entire file figuring out my, mt and topo_size
                status = 0
                do while (status == 0 .and. abs(t - t0) < 1e-15)
                    read(iunit,fmt=*,iostat=status) t, x, y
                    topo_size = topo_size + 1
                    if (y /= y_old .and. abs(t - t0) < 1e-15) then
                        dtopo_data%num_cells(2) = dtopo_data%num_cells(2) + 1
                        y_old = y
                    endif
                enddo
                dtopo_data%num_cells(1) = (topo_size - 1) / dtopo_data%num_cells(2)
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) t, x, y
                    topo_size = topo_size + 1
                enddo


                if (status > 0) then
                    print *, "IO error occured in ", dtopo_data%file_path, ", aborting!"
                    stop
                endif

                ! Calculate remaining values
                dtopo_data%num_cells(3) = (topo_size - 1) / (dtopo_data%num_cells(1) * dtopo_data%num_cells(2))
                dtopo_data%x_upper = x
                dtopo_data%y_lower = y
                dtopo_data%t_upper = t
                dtopo_data%dx = (dtopo_data%x_upper - dtopo_data%x_lower) / (dtopo_data%num_cells(1) - 1)
                dtopo_data%dy = (dtopo_data%y_upper - dtopo_data%y_lower) / (dtopo_data%num_cells(2) - 1)
                dtopo_data%dt = (dtopo_data%t_upper - dtopo_data%t_lower) / (dtopo_data%num_cells(3) - 1)

            ! New ASCII headered dtopo files, similar to topography files type
            ! 2 and 3
            case(2:3)
                ! Read in header directly
                read(iunit,*) dtopo_data%num_cells(1)
                read(iunit,*) dtopo_data%num_cells(2)
                read(iunit,*) dtopo_data%num_cells(3)
                read(iunit,*) dtopo_data%x_lower
                read(iunit,*) dtopo_data%y_lower
                read(iunit,*) dtopo_data%t_lower
                read(iunit,*) dtopo_data%dx
                read(iunit,*) dtopo_data%dy
                read(iunit,*) dtopo_data%dt

                dtopo_data%x_upper = dtopo_data%x_lower + dtopo_data%dx * (dtopo_data%num_cells(1) - 1)
                dtopo_data%y_upper = dtopo_data%y_lower + dtopo_data%dy * (dtopo_data%num_cells(2) - 1)
                dtopo_data%t_upper = dtopo_data%t_lower + dtopo_data%dt * (dtopo_data%num_cells(3) - 1)
            case default
                print *, 'ERROR:  Unrecognized topography type'
                print *, '    topo_type = ', dtopo_data%dtopo_type
                print *, '  for dtopo file:'
                print *, '   ', dtopo_data%file_path
                stop
        end select

         close(iunit)
    end subroutine read_dtopo_header

    ! ========================================================================
    !  read_dtopo(fname)
    ! ========================================================================
    subroutine read_dtopo(mx,my,mt,dtopo_type,fname,dtopo)

      implicit none

      ! Arguments
      integer, intent(in) :: mx,my,mt,dtopo_type
      character*150, intent(in) :: fname
      real(kind=8), intent(inout) :: dtopo(1:mx*my*mt)

      ! Local
      integer, parameter :: iunit = 29
      integer :: i,j,k,dtopo_size,status
      real(kind=8) :: t,x,y

      open(unit=iunit, file=fname, status = 'unknown',form='formatted')

      select case(abs(dtopo_type))
         case(1)
            ! ASCII file with 4 columns
            do i = 1,mx*my*mt
               read(iunit,fmt=*,iostat=status) t,x,y, dtopo(i)
            enddo

         case(2)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            ! read the data
            do i = 1,mx*my*mt
               read(iunit,*) dtopo(i)
            enddo
         case(3)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            do k = 1,mt
               do j = 1,my
                  read(iunit,*) (dtopo((k-1)*mx*my + (j-1)*mx + i) , i=1,mx)
               enddo
            enddo
      end select

    end subroutine read_dtopo

end module dtopo_module