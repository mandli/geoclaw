! This module provides utility functions for working with NetCDF files using the
! netcdf-fortran library.
module netcdf_utils_module
#ifdef NETCDF
    use netcdf

    implicit none
    private

    public :: nc_check
    public :: nc_open_readonly, nc_close
    public :: nc_has_att, nc_get_att_str, nc_get_att_int, nc_get_att_real
    public :: nc_get_var_units
    public :: nc_varid_by_name
    public :: nc_find_var_by_names
    public :: nc_find_var_by_att_str

contains

    subroutine nc_check(ierr, context)
        integer, intent(in) :: ierr
        character(len=*), intent(in), optional :: context
        if (ierr /= nf90_noerr) then
            if (present(context)) then
                write(*,*) "NetCDF error in: ", trim(context)
            else
                write(*,*) "NetCDF error:"
            end if
            write(*,*) "  ", trim(nf90_strerror(ierr))
            stop 1
        end if
    end subroutine nc_check

  subroutine nc_open_readonly(path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    call nc_check(nf90_open(trim(path), nf90_nowrite, ncid), "nf90_open("//trim(path)//")")
  end subroutine nc_open_readonly

  subroutine nc_close(ncid)
    integer, intent(inout) :: ncid
    if (ncid > 0) then
      call nc_check(nf90_close(ncid), "nf90_close")
      ncid = -1
    end if
  end subroutine nc_close

  logical function nc_has_att(ncid, varid, attname)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    integer :: ierr
    ierr = nf90_inquire_attribute(ncid, varid, trim(attname))
    nc_has_att = (ierr == nf90_noerr)
  end function nc_has_att

  subroutine nc_get_att_str(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    character(len=*), intent(out) :: value
    logical, intent(out), optional :: found
    integer :: ierr

    value = ""
    ierr = nf90_get_att(ncid, varid, trim(attname), value)
    if (present(found)) found = (ierr == nf90_noerr)
    if (.not. present(found)) call nc_check(ierr, "nf90_get_att("//trim(attname)//")")
  end subroutine nc_get_att_str

  subroutine nc_get_att_int(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    integer, intent(out) :: value
    logical, intent(out), optional :: found
    integer :: ierr

    value = 0
    ierr = nf90_get_att(ncid, varid, trim(attname), value)
    if (present(found)) found = (ierr == nf90_noerr)
    if (.not. present(found)) call nc_check(ierr, "nf90_get_att("//trim(attname)//")")
  end subroutine nc_get_att_int

  subroutine nc_get_att_real(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    real(kind=8), intent(out) :: value
    logical, intent(out), optional :: found
    integer :: ierr

    value = 0.0d0
    ierr = nf90_get_att(ncid, varid, trim(attname), value)
    if (present(found)) found = (ierr == nf90_noerr)
    if (.not. present(found)) call nc_check(ierr, "nf90_get_att("//trim(attname)//")")
  end subroutine nc_get_att_real

  subroutine nc_get_var_units(ncid, varid, units, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(out) :: units
    logical, intent(out), optional :: found
    call nc_get_att_str(ncid, varid, "units", units, found)
  end subroutine nc_get_var_units

  integer function nc_varid_by_name(ncid, name, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    logical, intent(out), optional :: found
    integer :: ierr, vid

    ierr = nf90_inq_varid(ncid, trim(name), vid)
    if (present(found)) then
      found = (ierr == nf90_noerr)
      if (found) nc_varid_by_name = vid
      if (.not. found) nc_varid_by_name = -1
    else
      call nc_check(ierr, "nf90_inq_varid("//trim(name)//")")
      nc_varid_by_name = vid
    end if
  end function nc_varid_by_name

  integer function nc_find_var_by_names(ncid, names, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: names(:)
    logical, intent(out) :: found
    integer :: k, vid
    logical :: ok

    found = .false.
    nc_find_var_by_names = -1
    do k = 1, size(names)
      vid = nc_varid_by_name(ncid, trim(names(k)), ok)
      if (ok) then
        found = .true.
        nc_find_var_by_names = vid
        return
      end if
    end do
  end function nc_find_var_by_names

  integer function nc_find_var_by_att_str(ncid, attname, attvalue, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: attname, attvalue
    logical, intent(out) :: found

    integer :: nvars, ndims, ngatts, unlimdim
    integer :: varid, ierr
    character(len=256) :: v
    logical :: has

    call nc_check(nf90_inquire(ncid, ndims, nvars, ngatts, unlimdim), "nf90_inquire")

    found = .false.
    nc_find_var_by_att_str = -1

    do varid = 1, nvars
      ! netcdf-fortran uses 1..nvars in inquiry? careful: varids are 0..nvars-1.
      ! We'll query by 0-based varid.
    end do

    do varid = 0, nvars-1
      has = nc_has_att(ncid, varid, trim(attname))
      if (.not. has) cycle
      call nc_get_att_str(ncid, varid, trim(attname), v, has)
      if (has) then
        if (trim(v) == trim(attvalue)) then
          found = .true.
          nc_find_var_by_att_str = varid
          return
        end if
      end if
    end do
  end function nc_find_var_by_att_str

#else
    implicit none
    write(*,*) "This program was built without NETCDF support."
    write(*,*) "Rebuild with -DNETCDF (or equivalent) and link netcdf-fortran."
    stop 2
#endif
end module netcdf_utils_module
