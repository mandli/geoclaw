! GeoClaw specific output....  add eta to q array before printing
!
! Write the results to the file fort.q<iframe>
! Use format required by matlab script  plotclaw2.m or Python tools
!
! set outaux = .true. to also output the aux arrays to fort.a<iframe>
!
!
subroutine valout (lst, lend, time, num_eqn, num_aux)

    use amr_module, only: alloc, store1, storeaux, levelptr, node, rnode
    use amr_module, only: ndilo, ndihi, ndjlo, cornxlo, ndjhi, cornylo
    use amr_module, only: hxposs, hyposs, lstart
    use amr_module, only: matlabu
    use amr_module, only: num_ghost => nghost
    use amr_module, only: output_format, output_aux_onlyonce, output_aux_components
    use amr_module, only: t0
    use amr_module, only: timeValout

    use geoclaw_module, only: coordinate_system

#ifdef NETCDF
    use netcdf
#endif
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: lst, lend, num_eqn, num_aux
    real(kind=8), intent(in) :: time

    ! Local
    logical :: output_aux
    integer :: output_aux_num, ios, m, i, j, loc, locaux
    integer :: grid_pointer, num_grids
    real(kind=8) :: h, hu, hv, eta
    real(kind=8), allocatable :: qeta(:)
    character(len=10) :: fname(5)

    integer :: level_pointer, level, num_cells(2), total_cells(2), num_dim
    real(kind=8) :: lower(2)

    ! Output units
    integer, parameter :: Q_UNIT = 50
    integer, parameter :: T_UNIT = 60
    integer, parameter :: A_UNIT = 70
    integer, parameter :: B_UNIT = 71

    ! Output formats
    character(len=*), parameter :: q_header_format = &
                    "(i5,'                 grid_number',/," // &
                     "i5,'                 AMR_level',/,"   // &
                     "i5,'                 mx',/,"          // &
                     "i5,'                 my',/,"          // &
                     "e18.8,'    xlow', /,"                 // &
                     "e18.8,'    ylow', /,"                 // &
                     "e18.8,'    dx', /,"                   // &
                     "e18.8,'    dy',/)"

    character(len=*), parameter :: q_header_format_1D = &
                    "(i5,'                 grid_number',/," // &
                     "i5,'                 AMR_level',/,"   // &
                     "i5,'                 mx',/,"          // &
                     "e18.8,'    xlow', /,"                 // &
                     "e18.8,'    dx', /)"

    character(len=*), parameter :: t_format = &
                    "(e18.8,'    time', /,"              // &
                     "i5,'                 meqn'/,"      // &
                     "i5,'                 ngrids'/,"    // &
                     "i5,'                 naux'/,"      // &
                     "i5,'                 ndim'/,"      // &
                     "i5,'                 nghost'/,/)"

    character(len=*), parameter :: console_format = &
                     "('AMRCLAW: Frame ',i4, "   // &
                     "' output files done at time t = ', d12.6,/)"

    ! Timing
    integer :: clock_start, clock_finish, clock_rate

    call system_clock(clock_start,clock_rate)

    ! How many aux components requested for output
    output_aux_num = 0
    do i = 1, num_aux
        output_aux_num = output_aux_num + output_aux_components(i)
    end do
        
    ! Currently outputs all aux components if any are requested
    ! TODO: Only output requested fields
    output_aux = ((output_aux_num > 0) .and.            &
                 ((.not. output_aux_onlyonce) .or. (time == t0)))

    ! Create file names for output
    fname(1) = 'fort.q' // zfill(matlabu)
    fname(2) = 'fort.t' // zfill(matlabu)
    fname(3) = 'fort.a' // zfill(matlabu)
    fname(4) = 'claw' // zfill(matlabu) // ".nc"
    fname(5) = 'fort.b' // zfill(matlabu)

    ! Open q file, note that even if we output in binary we still provide header
    ! information for each grid through the q file
    open(unit=Q_UNIT, file=fname(1), iostat=ios, status="unknown",  &
         action="write", form='formatted')
    if ( ios /= 0 ) then
        print *, "Error opening file ", fname(1),"!"
        stop
    end if

    if (output_format == 2) then

#ifdef NETCDF
        ! Open and set the relevant meta-data for NetCDF 4 files
        check_netcdf_error( nf90_create(fname(4), NF90_NETCDF4, ncid) )

        if (coordinate_system==2) then
            coord_sys = 'Latitude-Longitude'
            coord_units = 'degrees Lat-Lon'
        else
            coord_sys = 'Cartesian'
            coord_units = 'meters'
        endif

        ! Prepare NetCDF4 dataset (for NetCDF3 change NF90_NETCDF4 to NF90_NOCLOBBER)
        check_netcdf_error( nf90_create(fname(4), NF90_NETCDF4, ncid) )

        ! Define dimension for number of output variables
        check_netcdf_error( nf90_def_dim(ncid, 'num_vars', num_vars, num_vars_id) )

        ! Define time variable
        check_netcdf_error( nf90_def_var(ncid, 'time', NF90_FLOAT, time_id) )
        check_netcdf_error( nf90_put_att(ncid, time_id, 'units', 'seconds') )

        ! Assign known global attributes
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'time', time) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'num_dim', 2) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'amr_levels', mxnest) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'num_vars', num_vars) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'q_components', output_q_components) )
        if (output_aux) then
            check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'aux_components', output_aux_components) )
        else
            check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'aux_components', 0 * output_aux_components) )
        endif

        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'coord_sys', trim(coord_sys)) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'coord_units', trim(coord_units)) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'coord_xbounds', (/xlower, xupper/)) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'coord_ybounds', (/ylower, yupper/)) )
        check_netcdf_error( nf90_enddef(ncid) )

        ! Assign time
        check_netcdf_error( nf90_put_var(ncid,time_id,time) )

#else
        print *, "GeoClaw was not compiled with support for NetCDF IO.  Please"
        print *, "recompile and include the appropriate compiled flags to "
        print *, "enable NetCDF support."
        stop
#endif

    else if (output_format == 3) then
        open(unit=B_UNIT, file=fname(4), iostat=ios, status="unknown",  &
             action="write", access='stream')
        if ( ios /= 0 ) then
            print *, "Error opening file ", fname(4),"!"
            stop
        end if
    end if

    ! Output q fields - loop over grids at each level starting with the coarsest
    level = 0
    num_grids = 0
    do while (level <= lend)
        level = level + 1
        grid_pointer = lstart(level)
        do
            ! Check to see if there are any more grids and if not cycle back to
            ! level loop
            if (grid_pointer == 0) exit

            num_grids = num_grids + 1
            num_cells(1) = node(ndihi, grid_pointer) - node(ndilo, grid_pointer) + 1
            num_cells(2) = node(ndjhi, grid_pointer) - node(ndjlo, grid_pointer) + 1
            loc = node(store1, grid_pointer)
            locaux = node(storeaux, grid_pointer)
            total_cells(1) = num_cells(1) + 2 * num_ghost
            total_cells(2) = num_cells(2) + 2 * num_ghost

            lower(1) = rnode(cornxlo, grid_pointer)
            lower(2) = rnode(cornylo, grid_pointer)

            ! Output header info, the special condition is for 1D AMR
            if (num_cells(2) == 1) then
                write(Q_UNIT, q_header_format_1D) grid_pointer, level,          &
                                                  num_cells(1), lower(1),       &
                                                  hxposs(level)
            else
                write(Q_UNIT, q_header_format) grid_pointer, level, num_cells(1), &
                                               num_cells(2), lower(1), lower(2),  &
                                               hxposs(level), hyposs(level)
            end if

            ! ASCII output only
            if (output_format == 1) then

                do j = num_ghost + 1, total_cells(2) - num_ghost
                    do i = num_ghost + 1, total_cells(1) - num_ghost

                        do m = 1, num_eqn
                            ! Zero out small values
                            if (abs(alloc(iadd(m, i, j))) < 1d-90) then
                                alloc(iadd(m, i, j)) = 0.d0
                            end if
                        end do

                        ! GeoClaw specific output
                        ! Extract depth and momenta
                        h = alloc(iadd(1, i, j))
                        hu = alloc(iadd(2, i, j))
                        hv = alloc(iadd(3, i, j))

                        ! Calculate surfaces
                        eta = h + alloc(iaddaux(1, i, j))
                        if (abs(eta) < 1d-90) then
                            eta = 0.d0
                        end if

                        write(Q_UNIT, "(50e26.16)") h, hu, hv, eta
                    end do
                    write(Q_UNIT,*) ' '
                end do

            ! NetCDF Output
            else if (output_format == 2) then
#ifdef NETCDF
            check_netcdf_error( NF90_REDEF(ncid) )

            write(gridstr,'(I4.4)') grid_pointer

            ! Define dimensions for grid
            check_netcdf_error( nf90_def_dim(ncid, 'dimx_' // trim(gridstr), &
                                                         num_cells(1), dimxid) )
            check_netcdf_error( nf90_def_dim(ncid, 'dimy_' // trim(gridstr), &
                                                         num_cells(2), dimyid) )

            ! Define grid variables
            rcode = nf90_def_var(ncid, 'grid_' // trim(gridstr), NF90_FLOAT, &
                                 (/num_vars_id, dimxid, dimyid/), gridid)
            ! Assign grid attributes
            check_netcdf_error( nf90_put_att(ncid, gridid, 'gridno', grid_pointer) )
            check_netcdf_error( nf90_put_att(ncid, gridid, 'level', level) )
            dim_names = "['num_vars','dimx','dimy']"
            check_netcdf_error( nf90_put_att(ncid,gridid,'dim_names', TRIM(dim_names)) )
            check_netcdf_error( nf90_put_att(ncid,gridid,'num_vars' ,num_vars) )
            check_netcdf_error( nf90_put_att(ncid,gridid,'dimx.m', num_cells(1)) )
            check_netcdf_error( nf90_put_att(ncid,gridid,'dimy.m', num_cells(2)) )
            check_netcdf_error( nf90_put_att(ncid,gridid,'dimx.lower', lower(1)) )
            check_netcdf_error( nf90_put_att(ncid, gridid, 'dimy.lower', lower(2)) )
            check_netcdf_error( nf90_put_att(ncid, gridid, 'dimx.d', hxposs(level)) )
            check_netcdf_error( nf90_put_att(ncid, gridid, 'dimy.d', hyposs(level)) )
            
            ! Grid is defined, leave define mode
            rcode = nf90_enddef(ncid)


            ! Create and fill grid array
            allocate(grid(num_eqn, num_cells(1), num_cells(2)))
            do j = num_ghost + 1, total_cells(2) - num_ghost
                do i = num_ghost + 1, total_cells(1) - num_ghost
                    iq_store = 0
                    do m=1, num_eqn
                        if (abs(alloc(iadd(m, i, j))) < 1d-90) then
                            alloc(iadd(m, i, j)) = 0.d0
                        endif
                        if (output_q_components(m) > 0) then
                            iq_store = iq_store + 1
                            grid(iq_store, i - num_ghost, j - num_ghost) =  &
                                        alloc(iadd(m, i, j))
                        endif
                    enddo
                    iaux_store = output_q_num
                    do m = 1, num_aux
                        if (output_aux .and. output_aux_components(m) > 0) then
                            iaux_store = iaux_store + 1
                            grid(iaux_store, i - num_ghost, j - num_ghost) =  &
                                        alloc(iaddaux(m, i, j))
                        endif
                    enddo
                enddo
            enddo

            ! Save grid
            check_netcdf_error( nf90_put_var(ncid, gridid, grid, (/1,1,1/), (/num_vars, nx, ny/)) )
            if (allocated(grid)) deallocate(grid)
#endif
            ! Fortran direct binary format
            else if (output_format == 3) then

                ! Allocate array to hold eta
                allocate(qeta(4 * total_cells(1) * total_cells(2)))
                do j = 1, total_cells(2)
                    do i = 1, total_cells(1)
                        do m = 1, 3
                            qeta(iaddqeta(m, i, j)) = alloc(iadd(m, i, j))
                        end do
                        eta = alloc(iadd(1, i, j)) + alloc(iaddaux(1, i, j))
                        qeta(iaddqeta(4, i, j)) = eta
                    end do
                end do

                write(B_UNIT) qeta

                deallocate(qeta)

            else
                print *,"Unknown output format ", output_format,"."
                stop
            end if 

            ! Fetch next grid
            grid_pointer = node(levelptr, grid_pointer)

        end do ! Grid loop
    end do ! Level loop

    close(Q_UNIT)
    if (output_format == 3) then
        close(B_UNIT)
    end if

    ! Assign global attribute num_grids and variable time
#ifdef NETCDF
    if (output_format == 3) then
        check_netcdf_error( nf90_put_var(ncid, time_id, time) )
        check_netcdf_error( nf90_put_att(ncid, NF90_GLOBAL, 'num_grids', num_grids) )
        check_netcdf_error( nf90_close(ncid) )
    end if
#endif

    ! Output fort.t file
    open(unit=T_UNIT, file=fname(2), iostat=ios, status="unknown",   &
        action="write", form='formatted')
    if ( ios /= 0 ) then
        print *, "Error opening file ", fname(2),"!"
        stop
    end if

    if (num_cells(2) > 1) then
        num_dim = 2
    else
        num_dim = 1
    end if
    ! Need to include num_ghost in order to strip the ghost cells from the 
    ! binary output.  Also note the total fields is num_eqn + 1 due to the 
    ! addition of the eta field
    write(T_UNIT, t_format) time, num_eqn + 1, num_grids, num_aux, num_dim, &
                            num_ghost

    print console_format, matlabu, time

    close(T_UNIT)

    ! Finally output aux arrays
!     if (out_aux) then

!         level = lst
!         if (level > lend)

!     end if


! c       -------------------
! c       # output aux arrays
! c       -------------------

!         if (outaux) then
! c        # output aux array to fort.aXXXX

!          level = lst
!  165     if (level .gt. lend) go to 190
!             mptr = lstart(level)
!  170        if (mptr .eq. 0) go to 180
!               nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
!               ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
!               locaux  = node(storeaux,mptr)
!               mitot   = nx + 2*nghost
!               mjtot   = ny + 2*nghost


!           if (output_format == 1) then
!              open(unit=matunit3,file=fname3,status='unknown',
!      .            form='formatted')
!               if (ny.gt.1) then
!                   write(matunit3,1001) mptr, level, nx, ny
!                 else
! c                 # output in 1d format if ny=1:
!                   write(matunit3,1003) mptr, level, nx
!                 endif
!               xlow = rnode(cornxlo,mptr)
!               ylow = rnode(cornylo,mptr)
!               if (ny.gt.1) then
!                   write(matunit3,1002)
!      &              xlow,ylow,hxposs(level),hyposs(level)
!                 else
!                   write(matunit3,1004)
!      &              xlow,hxposs(level)
!                 endif

!              do j = nghost+1, mjtot-nghost
!                 do i = nghost+1, mitot-nghost
!                    do ivar=1,naux
!                       if (abs(alloc(iaddaux(ivar,i,j))) .lt. 1d-90) 
!      &                   alloc(iaddaux(ivar,i,j)) = 0.d0
!                    enddo
!                    write(matunit3,109) (alloc(iaddaux(ivar,i,j)), 
!      &                              ivar=1,naux)
!                 enddo
!                 write(matunit3,*) ' '
!              enddo
!             endif
            
!          if (output_format == 3) then
! c            # binary output          
!              open(unit=matunit3,file=fname3,status='unknown',
!      &               access='stream')
!              i1 = iaddaux(1,1,1)
!              i2 = iaddaux(naux,mitot,mjtot)
! c            # NOTE: we are writing out ghost cell data also, unlike ascii
!              write(matunit3) alloc(i1:i2)
!              endif


!             mptr = node(levelptr, mptr)
!             go to 170
!  180     level = level + 1
!          go to 165

!  190    continue
!         close(unit=matunit3)
!         endif !# end outputting aux array



!       matlabu = matlabu + 1

!       close(unit=matunit1)
!       close(unit=matunit2)
!       if (output_format == 3) then
!           close(unit=matunit4)
!           endif

    ! Increment frame counter
    matlabu = matlabu + 1

    ! Timing support
    call system_clock(clock_finish,clock_rate)
    timeValout = timeValout + clock_finish - clock_start

contains

    integer pure function iadd(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar, i ,j
        iadd = loc + ivar - 1 + num_eqn * ((j-1) * total_cells(1) + i - 1)
    end function iadd

    integer pure function iaddaux(iaux,i,j)
        implicit none
        integer, intent(in) :: iaux, i, j
        iaddaux = locaux + iaux - 1 + num_aux * (i - 1)             &
                         + num_aux * total_cells(1) * (j - 1)
    end function iaddaux

    integer pure function iaddqeta(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar, i, j
        iaddqeta = 1 + ivar - 1 + 4 * ((j - 1) * total_cells(1) + i - 1)
    end function iaddqeta

    ! Given an integer frame fill in 0s to the left up to length 4
    ! TODO: Make the length of the string with zeros variable
    character(len=4) function zfill(frame) result(output)

        implicit none
        integer, intent(in) :: frame

        integer :: step, digit, position

        output = "xxxx"

        step = frame
        do position = 4, 1, -1
            digit = mod(step, 10)
            output(position:position) = char(ichar('0') + digit)
            step = step / 10
        end do

    end function zfill

#ifdef NETCDF
    subroutine check_netcdf_error(ios)

        implicit none

        integer, intent(in) :: ios

        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            stop
        end if

    end subroutine check_netcdf_error
#endif

end subroutine valout