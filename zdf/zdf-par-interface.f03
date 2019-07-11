#include "os-preprocess.fpp"
#include "os-config.h"

module zdf_parallel

use zdf

integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte   = selected_int_kind(2)
integer, private, parameter :: p_int64  = selected_int_kind(10)

integer(c_int64_t), parameter :: ZDF_GET_OFFSET = -1

enum, bind(c)
  enumerator :: ZDF_MPI, ZDF_INDEPENDENT, ZDF_MPIIO_INDEPENDENT, ZDF_MPIIO_COLLECTIVE
endenum

type, bind(C) :: t_zdf_par_file
    type(t_zdf_file)   :: zdf_file
    integer(c_int)     :: iomode
    integer(c_int32_t) :: comm
    integer(c_int32_t) :: fh
    integer(c_int64_t) :: fpos
end type

interface zdf_mpi_type
    module procedure zdf_mpi_type
end interface

interface zdf_open_file
  module procedure zdf_par_open_file
end interface

interface zdf_close_file
  subroutine zdf_par_close_file_f( file ) bind(C, name = "zdf_par_close_file_f")
    import t_zdf_par_file
    type(t_zdf_par_file)    :: file
  end subroutine zdf_par_close_file_f
end interface

interface zdf_add
  module procedure zdf_par_add_string
  module procedure zdf_par_add_int32
  module procedure zdf_par_add_double
  module procedure zdf_par_add_iteration
  module procedure zdf_par_add_grid_info
  module procedure zdf_par_add_part_info
end interface

interface zdf_start_cdset
    module procedure zdf_par_start_cdset
end interface

interface zdf_par_write_cdset
  subroutine zdf_par_write_cdset_f( file, cdset, chunk, offset ) bind(C, name="zdf_par_write_cdset_f")
    import t_zdf_par_file, t_zdf_dataset, t_zdf_chunk, c_int64_t
    type(t_zdf_par_file)   :: file
    type(t_zdf_dataset)    :: cdset
    type(t_zdf_chunk)   :: chunk
    integer(c_int64_t), value :: offset
  end subroutine zdf_par_write_cdset_f
end interface

interface zdf_end_cdset
  subroutine zdf_par_end_cdset_f( file, dataset ) bind(C, name="zdf_par_end_cdset_f")
      import t_zdf_par_file, t_zdf_dataset
      type(t_zdf_par_file)   :: file
      type(t_zdf_dataset)    :: dataset
  end subroutine zdf_par_end_cdset_f
end interface

interface zdf_iomode_name
  module procedure zdf_iomode_name
end interface

contains

function zdf_iomode_name( iomode )

  implicit none

  integer, intent(in) :: iomode
  character(len=18) :: zdf_iomode_name

  select case(iomode)
  case (ZDF_MPI)
    zdf_iomode_name = "mpi"
  case (ZDF_INDEPENDENT)
    zdf_iomode_name = "independent"
  case (ZDF_MPIIO_INDEPENDENT)
    zdf_iomode_name = "mpi-io independent"
  case (ZDF_MPIIO_COLLECTIVE)
    zdf_iomode_name = "mpi-io collective"
  case default
    zdf_iomode_name = "invalid"
  end select

end function zdf_iomode_name

function zdf_mpi_type( data_type )

  use mpi

  implicit none

  integer, intent(in) :: data_type
  integer :: zdf_mpi_type

  select case ( data_type )
  case( zdf_null )
      zdf_mpi_type = MPI_DATATYPE_NULL
  case( zdf_int8 )
      zdf_mpi_type = MPI_BYTE
  case( zdf_uint8 )
      zdf_mpi_type = MPI_BYTE
  case( zdf_int16 )
      zdf_mpi_type = MPI_INTEGER2
  case( zdf_uint16 )
      zdf_mpi_type = MPI_INTEGER2
  case( zdf_int32 )
      zdf_mpi_type = MPI_INTEGER4
  case( zdf_uint32 )
      zdf_mpi_type = MPI_INTEGER4
  case( zdf_float32 )
      zdf_mpi_type = MPI_REAL
  case( zdf_int64 )
      zdf_mpi_type = MPI_INTEGER8
  case( zdf_uint64 )
      zdf_mpi_type = MPI_INTEGER8
  case( zdf_float64 )
      zdf_mpi_type = MPI_DOUBLE_PRECISION
  case default
      zdf_mpi_type = MPI_DATATYPE_NULL
  end select

end function


subroutine zdf_par_open_file( file, name, amode, comm, iomode )

  implicit none

  type(t_zdf_par_file), intent(inout) :: file
  character(len=*), intent(in) :: name
  integer, intent(in) :: amode
  integer, intent(in) :: comm
  integer, intent(in) :: iomode

  interface
    subroutine zdf_par_open_file_f( file, name, amode, comm, iomode ) bind(C, name = "zdf_par_open_file_f")
      import t_zdf_par_file, c_char, c_int32_t
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      integer(c_int32_t), value      :: amode
      integer(c_int32_t), value      :: comm
      integer(c_int32_t), value      :: iomode
    end subroutine zdf_par_open_file_f
  end interface

  call zdf_par_open_file_f( file, trim(name) // C_NULL_CHAR, amode, comm, iomode )

end subroutine zdf_par_open_file

subroutine zdf_par_add_string( file, name, str )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  character(len = *), intent(in) :: str

  interface
    subroutine zdf_par_add_string_f( file, name, str ) bind(C, name = "zdf_par_add_string_f")
      import t_zdf_par_file, c_char
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      character(kind=c_char)  :: str(*)
    end subroutine zdf_par_add_string_f
  end interface

  call zdf_par_add_string_f( file, trim(name) // C_NULL_CHAR, trim(str) // C_NULL_CHAR )

end subroutine zdf_par_add_string


subroutine zdf_par_add_int32( file, name, value )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  integer, intent(in) :: value

  interface
    subroutine zdf_par_add_int32_f( file, name, val ) bind(C, name = "zdf_par_add_int32_f")
      import t_zdf_par_file, c_char, c_int32_t
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      integer(c_int32_t), value  :: val
    end subroutine zdf_par_add_int32_f
  end interface

  call zdf_par_add_int32_f( file, trim(name) // C_NULL_CHAR, value )

end subroutine zdf_par_add_int32


subroutine zdf_par_add_double( file, name, value )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  real(p_double), intent(in) :: value

  interface
    subroutine zdf_par_add_double_f( file, name, val ) bind(C, name = "zdf_par_add_double_f")
      import t_zdf_par_file, c_char, c_double
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      real(c_double), value  :: val
    end subroutine zdf_par_add_double_f
  end interface

  call zdf_par_add_double_f( file, trim(name) // C_NULL_CHAR, value )

end subroutine zdf_par_add_double


subroutine zdf_par_add_iteration( file, name, iter )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  type(t_zdf_iteration), intent(in) :: iter

  interface
    subroutine zdf_par_add_iteration_f( file, name, iter ) bind(C, name = "zdf_par_add_iteration_f")
      import t_zdf_par_file, c_char, t_zdf_iteration_c
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      type(t_zdf_iteration_c) :: iter
    end subroutine zdf_par_add_iteration_f
  end interface

  type( t_zdf_iteration_c ) :: iter_c
  character(len = len( iter % time_units ) + 1 ), target :: time_units

  ! Prepare iteration info
  time_units = trim( iter % time_units ) // C_NULL_CHAR
  iter_c % n = iter % n
  iter_c % t = iter % t
  iter_c % time_units = c_loc( time_units )

  call zdf_par_add_iteration_f( file, trim(name) // C_NULL_CHAR, iter_c )

end subroutine zdf_par_add_iteration


subroutine zdf_par_add_grid_info( file, name, info )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  type(t_zdf_grid_info), intent(in) :: info

  interface
    subroutine zdf_par_add_grid_info_f( file, name, info ) bind(C, name = "zdf_par_add_grid_info_f")
      import t_zdf_par_file, c_char, t_zdf_grid_info_c
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      type(t_zdf_grid_info_c) :: info
    end subroutine zdf_par_add_grid_info_f
  end interface

  type(t_zdf_grid_info_c) :: info_c
  character(len = p_str_maxlen ), target :: label, units

  type(t_zdf_grid_axis_c), dimension(info % ndims), target :: axis
  character(len = p_str_maxlen ), dimension(info % ndims), target :: axis_label
  character(len = p_str_maxlen ), dimension(info % ndims), target :: axis_units

  integer :: i

  ! Prepare grid info
  info_c % ndims = info % ndims
  do i = 1, info % ndims
      info_c % count(i) = info % count(i)
  enddo

  ! Label and units
  label = trim( info % label ) // C_NULL_CHAR
  units = trim( info % units ) // C_NULL_CHAR
  info_c % label = c_loc( label )
  info_c % units = c_loc( units )

  ! Axis
  do i = 1, info % ndims
      axis(i) % type = info % axis(i) % type
      axis(i) % min  = info % axis(i) % min
      axis(i) % max  = info % axis(i) % max
      axis_label(i) = trim(info % axis(i) % label) // C_NULL_CHAR
      axis_units(i) = trim(info % axis(i) % units) // C_NULL_CHAR
      axis(i) % label = c_loc( axis_label(i) )
      axis(i) % units = c_loc( axis_units(i) )
  enddo

  info_c % axis = c_loc( axis )

  call zdf_par_add_grid_info_f( file, trim(name) // C_NULL_CHAR, info_c )

end subroutine zdf_par_add_grid_info


subroutine zdf_par_add_part_info( file, name, info )

  implicit none

  type(t_zdf_par_file), intent(in) :: file
  character(len = *), intent(in) :: name
  type(t_zdf_part_info), intent(in) :: info

  interface
    subroutine zdf_par_add_part_info_f( file, name, info ) bind(C, name = "zdf_par_add_part_info_f")
      import t_zdf_par_file, c_char, t_zdf_part_info_c
      type(t_zdf_par_file)    :: file
      character(kind=c_char)  :: name(*)
      type(t_zdf_part_info_c) :: info
    end subroutine zdf_par_add_part_info_f
  end interface

  type( t_zdf_part_info_c ) :: info_c
  character(len = p_str_maxlen ), target :: lname
  character(len = p_str_maxlen ), dimension(:), pointer :: quants, labels, units
  type(c_ptr), dimension(64), target :: pquants, plabels, punits

  integer :: i

  ! Prepare particles info
  lname = trim(info % name) // C_NULL_CHAR
  info_c % name = c_loc(lname)
  info_c % np = info % np
  info_c % nquants = info % nquants

  allocate( quants( info % nquants ))
  allocate( labels( info % nquants ))
  allocate( units( info % nquants ))

  do i = 1, info % nquants
      quants(i) = trim(info % quants(i)) // C_NULL_CHAR
      labels(i) = trim(info % labels(i)) // C_NULL_CHAR
      units(i)  = trim(info % units(i)) // C_NULL_CHAR

      pquants(i) = c_loc( quants(i) )
      plabels(i) = c_loc( labels(i) )
      punits(i)  = c_loc( units(i) )
  enddo

  info_c % quants = c_loc( pquants )
  info_c % labels = c_loc( plabels )
  info_c % units  = c_loc( punits )
  call zdf_par_add_part_info_f( file, trim(name) // C_NULL_CHAR, info_c )

end subroutine zdf_par_add_part_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_par_start_cdset( file, name, dataset )

    implicit none

    type(t_zdf_par_file), intent(in) :: file
    character(len=*), intent(in)     :: name
    type(t_zdf_dataset), intent(in)  :: dataset

    interface
        subroutine zdf_par_start_cdset_f( file, name, dataset ) bind(C, name="zdf_par_start_cdset_f")
            import t_zdf_par_file, c_char, t_zdf_dataset
            type(t_zdf_par_file)   :: file
            character(kind=c_char) :: name(*)
            type(t_zdf_dataset)    :: dataset
        end subroutine zdf_par_start_cdset_f
    end interface

    call zdf_par_start_cdset_f( file, trim(name) // C_NULL_CHAR, dataset )

end subroutine zdf_par_start_cdset


end module zdf_parallel