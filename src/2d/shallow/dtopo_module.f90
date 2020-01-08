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
        integer :: num_cells(2)
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
    type(topo_file_type), pointer, private :: topo_files_0(:)
    real(kind=8) :: max_dt_dtopo
    
    logical :: topo_finalized
    integer :: imovetopo, aux_finalized

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
            allocate(dtopo_files(num_dtopo_files))
            allocate(dtopo_order(num_dtopo_files))

            do i = 1, num_dtopo_files
                read(iunit,*) dtopo_files(i)%file_path
                read(iunit,*) dtopo_files(i)%dtopo_type,        &
                read(inunt,*) dtopo_files(i)%minleveldtopo,     &
                read(iunit,*) dtopo_files(i)%maxleveldtopo

                allocate(dtopo_data(i)%dz(:, :, :))

                write(GEO_PARM_UNIT,*) '   fname:', dtopo_files(i)%file_path
                write(GEO_PARM_UNIT,*) '   topo type:',dtopo_files(i)%dtopo_type
                write(GEO_PARM_UNIT,*) '   minlevel, maxlevel:'
                write(GEO_PARM_UNIT,*)  dtopo_files(i)%minleveldtopo,   &
                                        dtopo_files(i)%maxleveldtopo

                ! Read in header data
                call read_dtopo_header(dtopo_files(i))
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
                        area_i = dxdtopo(i) * dydtopo(i)
                        area_j = dxdtopo(j) * dydtopo(j)

                        ! If two files have the same resolution, order is
                        ! arbitrarily chosen
                        if ((area_i < area_j) .or.          &
                            ((area_i == area_j) .and. (j < i)) ) then
                            
                            finer_than = finer_than + 1
                        end if
                    endif
                enddo
                ! ifinerthan tells how many other files, file i is finer than
                rank = num_dtopo - finer_than
                dtopo_order(rank) = i
            enddo

            ! Read in dtopo data
            do i = 1, num_dtopo_files
                call read_dtopo(dtopo_files(i))
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
    subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi, &
        yhi,tf,dx,dy,dt)

        implicit none

        ! Input Arguments
        character*150, intent(in) :: fname
        integer, intent(in) :: topo_type

        ! Output Arguments
        integer, intent(out) :: mx,my,mt
        real(kind=8), intent(out) :: xlow,ylow,t0,xhi,yhi,tf,dx,dy,dt

        ! Locals
        integer, parameter :: iunit = 7
        integer :: topo_size,status
        real(kind=8) :: x,y,t,y_old,t_old
        logical :: found_file

        ! Open file
        inquire(file=fname,exist=found_file)
        if (.not.found_file) then
            print *, 'Missing dtopo file:'
            print *, '    ', fname
            stop
        endif
        open(unit=iunit,file=fname,status='unknown',form='formatted')

        select case(topo_type)
            ! Old style ASCII dtopo files
            case(1)
                ! Initial size variables
                topo_size = 0
                my = 1
                mt = 1

                ! Read in first values, determines xlow, yhi and t0
                read(iunit,*) t0,xlow,yhi
                topo_size = topo_size + 1
                t = t0
                y_old = yhi
                ! Go through entire file figuring out my, mt and topo_size
                status = 0
                do while (status == 0.and. t.eq.t0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                    if (y /= y_old .and. t.eq.t0 ) then
                        my = my + 1
                        y_old = y
                    endif
                enddo
                mx = (topo_size-1)/my
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                enddo


                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                mt = (topo_size-1)/ (my*mx)
                xhi = x
                ylow = y
                tf = t
                dx = (xhi-xlow) / (mx-1)
                dy = (yhi-ylow) / (my-1)
                dt = (tf - t0) / (mt-1)

            ! New ASCII headered dtopo files, similar to topography files type
            ! 2 and 3
            case(2:3)
                ! Read in header directly
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) mt
                read(iunit,*) xlow
                read(iunit,*) ylow
                read(iunit,*) t0
                read(iunit,*) dx
                read(iunit,*) dy
                read(iunit,*) dt

                xhi = xlow + dx*(mx-1)
                yhi = ylow + dy*(my-1)
                tf = t0 + dt*(mt-1)
            case default
                print *, 'ERROR:  Unrecognized topography type'
                print *, '    topo_type = ',topo_type
                print *, '  for dtopo file:'
                print *, '   ', fname
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