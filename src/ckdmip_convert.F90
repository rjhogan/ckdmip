! ckdmip_convert.F90 - CKDMIP tool for converting units
! Author: Robin Hogan <r.j.hogan@ecmwf.int>
! Copyright (C) 2020 ECMWF
!
! This program reads line-by-line absorption spectra in the form of
! layer optical depths, and converts to another format such as
! mass-absorption coefficient

program ckdmip_convert

  use parkind1,        only : jprb, jpib
  use easy_netcdf,     only : netcdf_file

  implicit none

#include "ckdmip_version.h"

  real(jprb), parameter :: AccelDueToGravity  = 9.80665_jprb ! m s-2
  real(jprb), parameter :: MolarMassAir       = 28.970_jprb  ! g mol-1

  type(netcdf_file)  :: in_file, out_file
  character(len=512) :: in_file_name, out_file_name, argument_str
  character(len=512) :: constituent_id, profile_type, experiment, experiment_id
  character(len=4000) :: history_str

  ! Number and index of command-line arguments, and number of bands
  integer :: narg, iarg_in, iarg_out, istatus, ncol, jcol, jlev, nspec, nlev

  ! Optical depth, to be converted to something else
  real(jprb), allocatable :: data(:,:)

  ! Chunk sizes in compression of "something else"
  integer :: chunksizes(3)

  real(jprb), allocatable :: pressure_hl(:,:), mole_fraction_fl(:,:)

  ! Molar mass of gas (g mol-1)
  real(jprb) :: molar_mass

  real(jprb) :: scaling

  ! Convert optical depth to (1) mass extinction coefficient or (2)
  ! molar extinction coefficient
  integer :: convert_type 

  ! Loop over arguments
  narg = command_argument_count()

  if (narg < 3) then
    call print_usage_and_exit()
  end if

  call get_command_argument(1, argument_str)
  if (trim(argument_str) == '--mass') then
    convert_type = 1
  else if (trim(argument_str) == '--molar') then
    convert_type = 2
  else
    error stop 'First argument must be --mass or --molar'
  end if

  iarg_in = 2
  iarg_out = 3

  call get_command_argument(iarg_in, in_file_name, status=istatus)
  if (istatus /= 0) then
    error stop 'Error: failed to read name of input file as string of length < 512'
  end if
  
  call get_command_argument(iarg_out, out_file_name, status=istatus)
  if (istatus /= 0) then
    error stop 'Error: failed to read name of output file as string of length < 512'
  end if
  
  call in_file%open(trim(in_file_name), iverbose=5)

  call in_file%get_global_attribute('constituent_id', constituent_id)
  if (trim(constituent_id) == 'h2o') then
    molar_mass = 18.0152833_jprb
  else if (trim(constituent_id) == 'co2') then
    molar_mass = 44.011_jprb
  else if (trim(constituent_id) == 'o3') then
    molar_mass = 47.9982_jprb
  else if (trim(constituent_id) == 'n2o') then
    molar_mass = 44.013_jprb
  else if (trim(constituent_id) == 'co') then
    molar_mass = 28.0101_jprb
  else if (trim(constituent_id) == 'ch4') then
    molar_mass = 16.043_jprb
  else if (trim(constituent_id) == 'o2') then
    molar_mass = 31.9988_jprb
  else if (trim(constituent_id) == 'cfc11') then
    molar_mass = 137.3686_jprb
  else if (trim(constituent_id) == 'cfc12') then
    molar_mass = 120.914_jprb
  else if (trim(constituent_id) == 'hcfc22') then
    molar_mass = 86.469_jprb
  else if (trim(constituent_id) == 'ccl4') then
    molar_mass = 153.823_jprb
  else if (trim(constituent_id) == 'no2') then
    molar_mass = 46.0055_jprb
  else
    write(*,'(a)') 'constituent_id "', trim(constituent_id), '" not understood'
    error stop
  end if

  ! Get global attributes
  call in_file%get_global_attribute('history', history_str)
  call in_file%get_global_attribute('profile_type', profile_type)
  call in_file%get_global_attribute('experiment', experiment)
  call in_file%get_global_attribute('experiment_id', experiment_id)

  ! Create output file
  call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=5)
  call out_file%copy_dimensions(in_file)
  call out_file%copy_variable_definition(in_file, 'level')
  call out_file%copy_variable_definition(in_file, 'half_level')
  call out_file%copy_variable_definition(in_file, 'pressure_hl')
  call out_file%copy_variable_definition(in_file, 'temperature_hl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_hl')
  call out_file%copy_variable_definition(in_file, 'pressure_fl')
  call out_file%copy_variable_definition(in_file, 'temperature_fl')
  call out_file%copy_variable_definition(in_file, 'mole_fraction_fl')
  call out_file%copy_variable_definition(in_file, 'reference_surface_mole_fraction')
  call out_file%copy_variable_definition(in_file, 'wavenumber')

  ! Put global attributes
  call out_file%put_global_attributes( &
       &  title_str = "Spectral absorption profiles of " // trim(to_upper(constituent_id)), &
       &  source_str= "Line-By-Line Radiative Transfer Model (LBLRTM)", &
       &  project_str="CKDMIP", &
       &  comment_str="LBLRTM computes layer optical depth from half-level pressure, temperature and mole fraction." &
       &  // new_line('a') // "Layer-mean temperature and mole fractions are computed from the " &
       &  //"half-level values assuming a linear variation with pressure.", &
       &  conventions_str="CF-1.7", prior_history_str=history_str)

  call out_file%put_global_attribute('constituent_id', trim(constituent_id))
  call out_file%put_global_attribute('profile_type', trim(profile_type))
  call out_file%put_global_attribute('experiment', trim(experiment))
  call out_file%put_global_attribute('experiment_id', trim(experiment_id))

  ! Choose chunk sizes for compression, allowing for the limit on the
  ! chunk size being 2**32 (which is apparently less)
  nspec = in_file%get_outer_dimension('wavenumber')
  nlev  = in_file%get_outer_dimension('level')
  chunksizes = [nspec,nlev,1]
  if (int(nspec,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
    chunksizes(2) = 1
  end if

  if (convert_type == 1) then
    call out_file%define_variable('mass_extinction_coefft', &
         &  dim3_name='column', dim2_name='level', dim1_name='wavenumber', &
         &  units_str='m2 kg-1', long_name='Mass extinction coefficient', &
             &   deflate_level=2, shuffle=.true., chunksizes=chunksizes)
  else
    call out_file%define_variable('molar_extinction_coefft', &
         &  dim3_name='column', dim2_name='level', dim1_name='wavenumber', &
         &  units_str='m2 mol-1', long_name='Molar extinction coefficient', &
             &   deflate_level=2, shuffle=.true., chunksizes=chunksizes)
  end if

  call out_file%copy_variable(in_file, 'level')
  call out_file%copy_variable(in_file, 'half_level')
  call out_file%copy_variable(in_file, 'pressure_hl')
  call out_file%copy_variable(in_file, 'temperature_hl')
  call out_file%copy_variable(in_file, 'mole_fraction_hl')
  call out_file%copy_variable(in_file, 'pressure_fl')
  call out_file%copy_variable(in_file, 'temperature_fl')
  call out_file%copy_variable(in_file, 'mole_fraction_fl')
  call out_file%copy_variable(in_file, 'reference_surface_mole_fraction')
  call out_file%copy_variable(in_file, 'wavenumber')

  call in_file%get('pressure_hl', pressure_hl)
  call in_file%get('mole_fraction_fl', mole_fraction_fl)

  ncol = in_file%get_outer_dimension('optical_depth')
  do jcol = 1,ncol
    call in_file%get('optical_depth', data, jcol)
    do jlev = 1,nlev
      ! "scaling" here has units of m2 mol-1 (per unit optical
      ! depth). The factor of 0.001 converts molar mass from g mol-1
      ! to kg mol-1.
      scaling = AccelDueToGravity * MolarMassAir * 0.001_jprb &
           &  / (mole_fraction_fl(jlev,jcol) * (pressure_hl(jlev+1,jcol) - pressure_hl(jlev,jcol)))
      if (convert_type == 2) then
        ! Convert "scaling" to units of m2 kg-1 (per unit optical
        ! depth)
        scaling = scaling / (molar_mass * 0.001_jprb)
      end if
      data(:,jlev) = data(:,jlev) * scaling
    end do
    if (convert_type == 1) then
      call out_file%put('mass_extinction_coefft', data, jcol)
    else
      call out_file%put('molar_extinction_coefft', data, jcol)
    end if
  end do

  call in_file%close()
  call out_file%close()

contains

  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: ckdmip_convert [--molar|--mass] <input.h5> <output.h5>'
    write(*,'(a)') '  Converts the the single-gas layer optical depth in the input file to either'
    write(*,'(a)') '  mass-extinction coefficient or molar extinction coefficient in the output file.'
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
  
end program ckdmip_convert
