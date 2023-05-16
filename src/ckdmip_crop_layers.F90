! ckdmip_crop_layers.F90 - CKDMIP tool for cropping layers
! Author: Robin Hogan <r.j.hogan@ecmwf.int>
! Copyright (C) 2023- ECMWF
!
! This program reads line-by-line absorption spectra in the form of
! layer optical depths, and converts to another format such as
! mass-absorption coefficient

program ckdmip_crop_layers

  use parkind1,        only : jprb, jpib
  use easy_netcdf,     only : netcdf_file

  implicit none

#include "ckdmip_version.h"

  real(jprb), parameter :: AccelDueToGravity  = 9.80665_jprb ! m s-2
  real(jprb), parameter :: MolarMassAir       = 28.970_jprb  ! g mol-1

  type(netcdf_file)  :: in_file, out_file
  character(len=512) :: in_file_name, out_file_name, argument_str
  character(len=512) :: constituent_id, profile_type, experiment, experiment_id
  character(len=512) :: new_experiment, new_experiment_id
  character(len=4000) :: history_str

  ! Number and index of command-line arguments, and number of bands
  integer :: narg, iarg_in, iarg_out, istatus, ncol, jcol, jlev, nwav, nlev, nlevnew

  ! Layer optical depth of a single column, dimensioned (nwav,nlev)
  real(jprb), allocatable :: layer_od(:,:)

  ! Chunk sizes in compression of "something else"
  integer :: chunksizes(3)

  real(jprb), allocatable :: pressure_hl(:,:), temperature_hl(:,:)
  real(jprb), allocatable :: pressure_fl(:,:), temperature_fl(:,:)
  real(jprb), allocatable :: mole_fraction_fl(:,:)

  ! Temporary arrays
  integer,    allocatable :: int_vector(:)
  real(jprb), allocatable :: real_matrix(:,:)
  
  ! Molar mass of gas (g mol-1)
  real(jprb) :: molar_mass

  ! Number of layers to crop from top and bottom of atmosphere
  integer :: ncroptop, ncropbottom

  ! Argument number
  integer :: iarg

  ncroptop    = 0
  ncropbottom = 0
  
  ! Loop over arguments
  narg = command_argument_count()

  iarg = 1
  
  if (narg < 4) then
    call print_usage_and_exit()
  end if

  do while (iarg <= narg)
    
    call get_command_argument(iarg, argument_str)
    if (trim(argument_str) == '--ncropbottom') then
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str,*) ncropbottom
    else
      exit
    end if
    iarg = iarg+1

  end do

  iarg_in  = iarg
  iarg_out = iarg+1

  call get_command_argument(iarg_in, in_file_name, status=istatus)
  if (istatus /= 0) then
    error stop '*** Error: failed to read name of input file as string of length < 512'
  end if
  
  call get_command_argument(iarg_out, out_file_name, status=istatus)
  if (istatus /= 0) then
    error stop '*** Error: failed to read name of output file as string of length < 512'
  end if

  ! Open input file  
  call in_file%open(trim(in_file_name), iverbose=3)

  ! Get global attributes
  call in_file%get_global_attribute('history', history_str)
  call in_file%get_global_attribute('profile_type', profile_type)
  call in_file%get_global_attribute('experiment', experiment)
  call in_file%get_global_attribute('experiment_id', experiment_id)
  call in_file%get_global_attribute('constituent_id', constituent_id)

  ! Create output file
  call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=5)

  ncol = in_file%get_outer_dimension('optical_depth')
  nlev = in_file%get_outer_dimension('level')
  nwav = in_file%get_outer_dimension('wavenumber')
  nlevnew = nlev-ncroptop-ncropbottom
  !call out_file%copy_dimensions(in_file)
  call out_file%define_dimension("column", ncol)
  call out_file%define_dimension("level", nlevnew)
  call out_file%define_dimension("half_level", nlevnew+1)
  call out_file%define_dimension("wavenumber", nwav)
    
  call out_file%copy_variable_definition(in_file, 'level')
  call out_file%copy_variable_definition(in_file, 'half_level')
  call out_file%copy_variable_definition(in_file, 'pressure_hl')
  call out_file%copy_variable_definition(in_file, 'temperature_hl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_hl')
  call out_file%copy_variable_definition(in_file, 'pressure_fl')
  call out_file%copy_variable_definition(in_file, 'temperature_fl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_fl')
  if (in_file%exists('reference_surface_mole_fraction')) then
    call out_file%copy_variable_definition(in_file, 'reference_surface_mole_fraction')
  end if
  call out_file%copy_variable_definition(in_file, 'wavenumber')

  ! Put global attributes
  call out_file%put_global_attributes( &
       &  title_str = "Spectral optical depth profiles of " // trim(to_upper(constituent_id)), &
       &  inst_str  = "European Centre for Medium-Range Weather Forecasts (ECMWF)", &
       &  source_str= "Line-By-Line Radiative Transfer Model (LBLRTM)", &
       &  creator_name="Robin Hogan", &
       &  creator_email_str="r.j.hogan@ecmwf.int", &
       &  contributor_name="Marco Matricardi", &
       &  project_str="CKDMIP", &
       &  comment_str="LBLRTM computes layer optical depth from half-level pressure, temperature and mole fraction." &
       &  // new_line('a') // "Layer-mean temperature and mole fractions are computed from the " &
       &  //"half-level values assuming a linear variation with pressure.", &
       &  conventions_str="CF-1.7", prior_history_str=history_str)

  call out_file%put_global_attribute('constituent_id', trim(constituent_id))
  call out_file%put_global_attribute('profile_type', trim(profile_type))

  write(new_experiment,'(a,a,i0,a)') trim(experiment), ', cropping ', ncropbottom, ' layers from bottom of atmosphere'
  write(new_experiment_id,'(a,a,i0)') trim(experiment_id), '-cropsurf', ncropbottom

  call out_file%put_global_attribute('experiment',    trim(new_experiment))
  call out_file%put_global_attribute('experiment_id', trim(new_experiment_id))

  ! Choose chunk sizes for compression, allowing for the limit on the
  ! chunk size being 2**32 (which is apparently less)
  chunksizes = [nwav,nlevnew,1]
  if (int(nwav,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
    chunksizes(2) = 1
  end if

  call out_file%define_variable('optical_depth', &
       &  dim3_name='column', dim2_name='level', dim1_name='wavenumber', &
       &  long_name='Layer optical depth', &
       &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)

  call  in_file%get('level', int_vector)
  call out_file%put('level', int_vector(1+ncroptop:nlev-ncropbottom))
  deallocate(int_vector)
  
  call  in_file%get('half_level', int_vector)
  call out_file%put('half_level', int_vector(1+ncroptop:nlev+1-ncropbottom))
  deallocate(int_vector)

  call  in_file%get('pressure_hl', real_matrix)
  call out_file%put('pressure_hl', real_matrix(1+ncroptop:nlev+1-ncropbottom,:))
  deallocate(real_matrix)

  call  in_file%get('temperature_hl', real_matrix)
  call out_file%put('temperature_hl', real_matrix(1+ncroptop:nlev+1-ncropbottom,:))
  deallocate(real_matrix)

  call  in_file%get('mole_fraction_hl', real_matrix)
  call out_file%put('mole_fraction_hl', real_matrix(1+ncroptop:nlev+1-ncropbottom,:))
  deallocate(real_matrix)

  call  in_file%get('pressure_fl', real_matrix)
  call out_file%put('pressure_fl', real_matrix(1+ncroptop:nlev-ncropbottom,:))
  deallocate(real_matrix)

  call  in_file%get('temperature_fl', real_matrix)
  call out_file%put('temperature_fl', real_matrix(1+ncroptop:nlev-ncropbottom,:))
  deallocate(real_matrix)

  call  in_file%get('mole_fraction_fl', real_matrix)
  call out_file%put('mole_fraction_fl', real_matrix(1+ncroptop:nlev-ncropbottom,:))
  deallocate(real_matrix)

  if (in_file%exists('reference_surface_mole_fraction')) then
    call out_file%copy_variable(in_file, 'reference_surface_mole_fraction')
  end if
  call out_file%copy_variable(in_file, 'wavenumber')

  do jcol = 1,ncol
    call  in_file%get('optical_depth', real_matrix, jcol)
    call out_file%put('optical_depth', real_matrix(:,1+ncroptop:nlev-ncropbottom), jcol)
  end do

  call  in_file%close()
  call out_file%close()

contains

  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: ckdmip_crop_layers [--ncroptop <ntopcrop>] [--ncropbottom <ncropbottom>] <input.h5> <output.h5>'
    write(*,'(a)') '  Copies the file optionally cropping layers from the top and/or bottom'
    stop

  end subroutine print_usage_and_exit

  ! Convert lower case string to upper case.
  ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
  ! Original author: Clive Page
  function to_upper(str_in) result(str_out)
    
    implicit none
    
    character(len=*), intent(in) :: str_in
    character(len=len(str_in)) :: str_out
    integer :: i,j
    
    do i = 1, len(str_in)
      j = iachar(str_in(i:i))
      if (j>= iachar("a") .and. j<=iachar("z") ) then
        str_out(i:i) = achar(iachar(str_in(i:i))-32)
      else
        str_out(i:i) = str_in(i:i)
      end if
    end do
    
  end function to_upper
  
end program ckdmip_crop_layers
