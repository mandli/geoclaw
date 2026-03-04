! Stub version of netcdf_utils_module, for use when NETCDF support is not available.
module netcdf_utils_module

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
    end subroutine nc_check

     subroutine nc_open_readonly(path, ncid)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ncid
    end subroutine nc_open_readonly

    subroutine nc_close(ncid)
        integer, intent(inout) :: ncid
    end subroutine nc_close

    logical function nc_has_att(ncid, varid, attname)
        integer, intent(in) :: ncid, varid
        character(len=*), intent(in) :: attname
    end function nc_has_att

  subroutine nc_get_att_str(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    character(len=*), intent(out) :: value
    logical, intent(out), optional :: found
  end subroutine nc_get_att_str

  subroutine nc_get_att_int(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    integer, intent(out) :: value
    logical, intent(out), optional :: found
  end subroutine nc_get_att_int

  subroutine nc_get_att_real(ncid, varid, attname, value, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: attname
    real(kind=8), intent(out) :: value
    logical, intent(out), optional :: found
  end subroutine nc_get_att_real

  subroutine nc_get_var_units(ncid, varid, units, found)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(out) :: units
    logical, intent(out), optional :: found
  end subroutine nc_get_var_units

  integer function nc_varid_by_name(ncid, name, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    logical, intent(out), optional :: found
  end function nc_varid_by_name

  integer function nc_find_var_by_names(ncid, names, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: names(:)
    logical, intent(out) :: found
  end function nc_find_var_by_names

  integer function nc_find_var_by_att_str(ncid, attname, attvalue, found)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: attname, attvalue
    logical, intent(out) :: found
  end function nc_find_var_by_att_str

end module netcdf_utils_module
