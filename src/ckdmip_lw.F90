! ckdmip_lw.F90 - longwave radiative transfer for CKDMIP
!
! Copyright (C) 2019- ECMWF
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author: Robin Hogan <r.j.hogan@ecmwf.int>
!
! This program performs radiative transfer calculations on
! line-by-line absorption spectra, or on CKD optical depths for a
! number of g points. Run with no arguments to see usage information.

program ckdmip_lw

  use parkind1,        only : jprb
  use easy_netcdf,     only : netcdf_file
  use longwave_config, only : longwave_config_type
  use planck_function, only : calc_planck_function, calc_planck_function1
  use longwave_fluxes, only : calc_longwave_fluxes_n, MAX_ANGLES

  implicit none

#include "ckdmip_version.h"

  real(jprb), parameter :: PI = 3.14159265358979323846_jprb

  integer,      parameter :: N_MAX_FILES  = 32
  character(5), parameter :: LW_UNITS_STR = 'W m-2'

  integer, parameter :: I_AS_IS = 0
  integer, parameter :: I_SCALE = 1
  integer, parameter :: I_CONC  = 2
  integer, parameter :: I_CONST = 3

  ! How to treat concentrations (0-4 above), a function of FILE
  integer :: conc_operator(N_MAX_FILES)

  ! Concentration factor, interpreted according to conc_operator, a function of FILE
  real(jprb) :: conc_factor(N_MAX_FILES)

  ! Is a particular file itself a merge of several gases? A function of FILE
  integer :: ngas_in_file(N_MAX_FILES)

  ! Input files
  type(netcdf_file) :: in_file(N_MAX_FILES)

  ! Reference surface concentration, a function of GAS
  real(jprb) :: ref_conc(N_MAX_FILES)

  ! Target surface concentration, a function of GAS
  real(jprb) :: target_conc(N_MAX_FILES)

  type(netcdf_file) :: out_file

  type(longwave_config_type) :: config

  ! Number of levels, spectral intervals, gases, profiles and files
  integer :: nlev, nspecin, ngas, ncol, nspec, ninspec, nfile, ngaslocal

  ! Number of blocks into which to divide the calculation
  integer :: nblock

  ! Loop indices
  integer :: jlev, jspec, igas, jcol, jblock, jarg, jband, jfile

  ! Range of spectral indices to process
  integer :: istartspec,iendspec

  ! Range of columns to process
  integer :: istartcol, iendcol

  ! Number and index to output column
  integer ::ncolout, icolout, iendcolout
  
  ! Pressure and temperature at half levels (nlev+1,ncol)
  real(jprb), allocatable, dimension(:,:) :: pressure_hl, temperature_hl
  ! Planck function at half levels (nspec,nlev+1)
  real(jprb), allocatable, dimension(:,:) :: planck_hl
  ! Wavenumber in cm-1 (nspec)
  real(jprb), allocatable, dimension(:)   :: wavenumber_cm1
  ! Optical depth of 1 gas (nspecin,nlev) and all gases (nspec,nlev)
  real(jprb), allocatable, dimension(:,:) :: od_1gas
  real(jprb), allocatable, dimension(:,:) :: od
  ! Longwave surface emission in W m-2 (nspec) and emissivity
  real(jprb), allocatable, dimension(:) :: surf_emission
  real(jprb), allocatable, dimension(:) :: surf_emissivity

  ! Surface temperature may be provided (K)
  real(jprb), allocatable, dimension(:) :: surf_temperature

  ! Fluxes (W m-2)
  real(jprb), allocatable, dimension(:,:) :: flux_dn, flux_up

  ! Broadband fluxes (W m-2)
  real(jprb), allocatable, dimension(:) :: bb_flux_dn, bb_flux_up

  ! Mole fraction (mol mol-1), dimensioned (ngas,nlev)
  real(jprb), allocatable, dimension(:,:) :: mole_fraction_fl

  ! Mole fraction (mol mol-1) of a single gas (nlev)
  real(jprb), allocatable, dimension(:) :: mole_fraction_fl_1gas

  ! Mole fractions (mol mol-1) of multiple gases in a single file (ngas_in_file,nlev)
  real(jprb), allocatable, dimension(:,:) :: mole_fraction_fl_ngas

  ! Concentration scaling as a function of pressure
  real(jprb), allocatable, dimension(:) :: mole_fraction_scaling

  ! Band boundaries
  real(jprb), allocatable, dimension(:) :: band_wavenumber1, band_wavenumber2

  ! Number of band to which each wavenumber belongs
  integer, allocatable, dimension(:) :: band_number

  ! Fluxes in bands (W m-2), dimensioned (nband, nlev+1)
  real(jprb), allocatable, dimension(:,:) :: band_flux_dn, band_flux_up

  ! Reference surface concentrations for all the gases in a merge file
  real(jprb), allocatable, dimension(:) :: ref_conc_merge

  ! Number and index of command-line arguments, and number of bands
  integer :: narg, iarg, nbands

  integer :: istatus, itoken

  ! File name
  character(len=512) :: argument_str, file_name, out_file_name

  character(len=512) :: molecule_str, molecule_list, scenario_str
  character(len=64)  :: spectral_dim_name, profile_type

  ! Do we output anything a function of wavenumber or g-point?
  logical :: do_save_spectrum

  ! Do we simply merge optical depths without any radiative transfer?
  logical :: do_merge_only

  ! Is there a scenario string?
  logical :: is_scenario

  ! Is the input in the form of optical depths and Planck functions in
  ! the g-points of a CKD model?
  logical :: is_ckd

  ! Initialize configuration variables
  do_merge_only = .false.
  is_scenario   = .false.
  is_ckd        = .false.
  istartcol     = 1
  iendcol       = 0
  nbands        = 0

  ! Initialize concentration scaling variables
  conc_operator = I_AS_IS   ! No scaling by default
  conc_factor   = -1.0_jprb ! No scaling factor by default
  ref_conc      = -1.0_jprb 
  target_conc   = -1.0_jprb
  ngas_in_file  = 1

  ! Initialize gas and argument counters
  nfile = 0
  ngas  = 0
  iarg  = 1

  ! Loop over arguments
  narg = command_argument_count()

  do while (iarg <= narg)
    
    call get_command_argument(iarg, argument_str, status=istatus)
    if (istatus /= 0) then
      write(*,'(a)') 'Failed to read argument as string of length < 512'
      call print_usage_and_exit()
    end if

    if (trim(argument_str) == "-c" .or. trim(argument_str) == "--config") then

      ! Read configuration namelist
      if (iarg == narg) then
        write(*,'(a)') '"-c" or "--config" argument must be followed by the name of a namelist file'
        call print_usage_and_exit()
      end if
      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        error stop 'Error: failed to read name of namelist file as string of length < 512'
      end if
      if (config%iverbose >= 2) then
        write(*,'(a,a)') 'Reading namelist file ', trim(file_name)
      end if
      call config%read(trim(file_name))

      ! Extract band structure, if provided
      nbands = config%get_bands(band_wavenumber1, band_wavenumber2)

    else if (argument_str == '-o' .or. argument_str == '--output') then

      ! Read output file name
      if (iarg == narg) then
        write(*,'(a)') '"-o" or "--output" argument must be followed by the name of the output file'
        call print_usage_and_exit()
      end if
      iarg = iarg+1
      call get_command_argument(iarg, out_file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of output file as string of length < 512'
        call print_usage_and_exit()
      end if
      ! Open output file
      call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=config%iverbose)

    else if (argument_str == '-m' .or. argument_str == '--merge-only') then

      do_merge_only = .true.
  
    else if (argument_str == '--scenario') then

      is_scenario = .true.
      iarg = iarg+1
      call get_command_argument(iarg, scenario_str)

    else if (argument_str == '--column-range') then

      if (iarg >= narg-1) then
        write(*,'(a)') '"--column-range" must be followed by two arguments'
        call print_usage_and_exit()
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) istartcol

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) iendcol

      if (istartcol < 1 .or. iendcol < istartcol) then
        write(*,'(a,i0,a,i0,a)') 'Error: column range of ', istartcol, ' to ', iendcol, ' is invalid'
        error stop
      end if

      if (config%iverbose >= 3) then
        write(*,'(a,i0,a,i0)') 'Processing columns ', istartcol, ' to ', iendcol
      end if

    else if (argument_str == '--scale' .or. &
         &   argument_str == '--const' .or. &
         &   argument_str == '--conc') then
      ! Open file containing the spectral optical depth of a single
      ! gas

      if (iarg >= narg-1) then
        write(*,'(a)') 'Error: "--scale", "--const" and "--conc" must be followed by two arguments'
        call print_usage_and_exit()
      else if (is_ckd) then
        write(*,'(a)') 'Error: cannot read spectral optical depth files if already using CKD optical depth file'
        error stop
      end if

      nfile = nfile+1
      ngas  = ngas+1
      if (argument_str == '--scale') then
        conc_operator(nfile) = I_SCALE
      else if (argument_str == '--const') then
        conc_operator(nfile) = I_CONST
      else if (argument_str == '--conc') then
        conc_operator(nfile) = I_CONC
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) conc_factor(nfile)

      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of input file as string of length < 512'
        call print_usage_and_exit()
      end if

      call in_file(nfile)%open(trim(file_name), iverbose=config%iverbose)

      call in_file(nfile)%get_global_attribute("profile_type", profile_type)
      if (trim(profile_type) == "merge") then
        error stop 'Error: cannot scale a file containing a merge of different gases'
      end if

    else if (argument_str == '--ckd') then
      ! Read a Correlated K-Distribution file containing optical depth
      ! and Planck function in a number of g-points
      if (nfile /= 0) then
        write(*,'(a)') 'Error: cannot read CKD optical depth file if already using spectral optical depth file(s)'
        error stop
      else if (iarg == narg) then
        write(*,'(a)') 'Error: "--ckd" must be followed by a file name'
        error stop
      end if

      is_ckd = .true.

      nfile = nfile+1
      ngas  = ngas+1
      iarg  = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of CKD input file as string of length < 512'
        call print_usage_and_exit()
      end if
      call in_file(nfile)%open(trim(file_name), iverbose=config%iverbose)

    else if (argument_str(1:1) == '-') then
      
      write(*,'(a,a,a)') 'Argument "', trim(argument_str), '" not understood'
      call print_usage_and_exit()

    else
      ! Open file containing the spectral optical depth of a single
      ! gas (or could be a merge)
      if (is_ckd) then
        write(*,'(a)') 'Error: cannot read spectral optical depth files if already using CKD optical depth file'
        error stop
      end if

      nfile = nfile+1

      file_name = trim(argument_str)

      call in_file(nfile)%open(trim(file_name), iverbose=config%iverbose)

      conc_operator(nfile) = I_AS_IS
      conc_factor(nfile)   = 1.0_jprb

      call in_file(nfile)%get_global_attribute("profile_type", profile_type)

      if (trim(profile_type) == 'merge') then
        call in_file(nfile)%get("reference_surface_mole_fraction", ref_conc_merge)
        ngas_in_file(nfile) = size(ref_conc_merge)
        ngas = ngas+size(ref_conc_merge)
      else
        ngas = ngas+1
      end if

    end if

    iarg = iarg + 1

  end do

  if (nfile == 0 .or. narg == 0) then
    call print_usage_and_exit()
  end if

  if (is_ckd .and. do_merge_only) then
    error stop 'Error: "--ckd" and "--merge-only" arguments are incompatible'
  end if

  if (is_ckd) then
    spectral_dim_name = 'gpoint_lw'
  else
    spectral_dim_name = 'wavenumber'    
  end if

  if (.not. is_ckd) then
    ! Report what gases will be merged

    igas = 1
    do jfile = 1,nfile
      if (in_file(jfile)%global_attribute_exists("constituent_id")) then
        call in_file(jfile)%get_global_attribute("constituent_id", molecule_str)
      else
        write(*,'(a,i0)') 'Error: constituent_id not present for file ', jfile
        error stop
      end if
      if (jfile == 1) then
        molecule_list = molecule_str
      else
        molecule_list = trim(molecule_list) // ' ' // trim(molecule_str)
      end if
      
      if (ngas_in_file(jfile) == 1) then
        ! Single gas
        ngaslocal = 1
        ! Load reference concentration, if present
        if (in_file(jfile)%exists("reference_surface_mole_fraction")) then
          call in_file(jfile)%get("reference_surface_mole_fraction", ref_conc(igas))
          target_conc(igas) = ref_conc(igas)
        else if (conc_operator(jfile) == I_CONC .or. conc_operator(jfile) == I_CONST) then
          write(*,'(a,a,a)') 'Error: cannot set concentration of ', to_upper(trim(molecule_str)), &
               &         '  as no reference surface concentration is present'
          error stop
        end if
      else
        ! Multiple gases
        call in_file(jfile)%get("reference_surface_mole_fraction", ref_conc_merge)
        ref_conc   (igas:igas+ngas_in_file(jfile)-1) = ref_conc_merge
        target_conc(igas:igas+ngas_in_file(jfile)-1) = ref_conc_merge
      end if
      
      if (config%iverbose >= 3) then
        write(*,'(a,i0,a,a)') '  File ', jfile, ': ', to_upper(trim(molecule_str))
        if (conc_operator(jfile) == I_AS_IS) then
          write(*,'(a)') '    No scaling to be applied'
        else if (conc_operator(jfile) == I_SCALE) then
          write(*,'(a,e11.5)') '    Scaling by a factor of ', conc_factor(jfile)
        else if (conc_operator(jfile) == I_CONC) then
          write(*,'(a,e11.5,a,e11.5,a)') '    Setting surface mole fraction to ', conc_factor(jfile), ' mol mol-1 ' &
               &             //'while reference value is ', ref_conc(jfile), ' mol mol-1:'
          write(*,'(a,f8.3)') '    Scaling by a factor of ', conc_factor(jfile) / ref_conc(igas)
        else if (conc_operator(jfile) == I_CONST) then
          write(*,'(a,e11.5,a)') '    Setting mole fraction to a constant value of ', conc_factor(jfile), ' mol mol-1'
        end if
      end if
      
      ! If a concentration has been specified then change it to a
      ! scaling
      if (conc_operator(jfile) == I_CONC) then
        conc_factor(jfile) = conc_factor(jfile) / ref_conc(igas)
        conc_operator(jfile) = I_SCALE
      end if
      if (conc_operator(jfile) == I_SCALE) then
        target_conc(igas) = ref_conc(igas) * conc_factor(jfile)
      else if (conc_operator(jfile) == I_CONST) then
        target_conc(igas) = conc_factor(jfile)
      end if
      
      igas = igas+ngaslocal

    end do
    
  end if

  ! Read basic properties from the first input file

  call in_file(1)%get(trim(config%pressure_name), pressure_hl)
  if (config%pressure_scaling /= 1.0_jprb) then
    pressure_hl = pressure_hl * config%pressure_scaling
  end if

  ncol = size(pressure_hl,2)
  nlev = size(pressure_hl,1)-1

  if (is_ckd) then
    config%do_read_planck = .true.
  end if

  ! Don't output bands if we are only merging spectra
  if (do_merge_only) then
    nbands = 0
  end if

  if (.not. config%do_read_planck) then
    ! Spectral files
    call in_file(1)%get(trim(config%wavenumber_name), wavenumber_cm1)
    if (config%nspectralstride == 1) then
      nspec = size(wavenumber_cm1)
      ninspec = nspec
    else
      nspec = (size(wavenumber_cm1)+config%nspectralstride-1) / config%nspectralstride
      ninspec = 1+(nspec-1)*config%nspectralstride
      write(*,'(a,i0)') 'Wavenumber stride: ', config%nspectralstride
      do jspec = 2,nspec
        wavenumber_cm1(jspec) = wavenumber_cm1(1+(jspec-1)*config%nspectralstride)
      end do
    end if
      
    call in_file(1)%get(trim(config%temperature_name), temperature_hl)

    if (config%surf_temperature >= 0.0_jprb) then
      write(*,'(a,f6.2,a)') 'Surface temperature fixed at ', config%surf_temperature, ' K'
    else if (in_file(1)%exists(trim(config%surf_temperature_name))) then
      call in_file(1)%get(trim(config%surf_temperature_name), surf_temperature)
      write(*,'(a,a)') 'Surface temperature read from variable ', trim(config%surf_temperature_name)
    else
      write(*,'(a)') 'Surface temperature assumed to equal lowest air temperature'
    end if

    ! Find locations of the requested bands
    if (nbands > 1) then
      allocate(band_number(nspec))
      band_number = 0
      allocate(band_flux_dn(nbands,nlev+1))
      allocate(band_flux_up(nbands,nlev+1))
      do jband = 1,nbands
        where (wavenumber_cm1 > band_wavenumber1(jband) .and. wavenumber_cm1 <= band_wavenumber2(jband))
          band_number = jband
        end where
      end do
    end if

  else
    ! CKD file
    call in_file(1)%get(trim(config%planck_name), planck_hl, 1)
    nspec = size(planck_hl,1)
    deallocate(planck_hl)
    ninspec = nspec
  end if

  allocate(surf_emission(nspec))
  allocate(surf_emissivity(nspec))
  surf_emissivity = 1.0_jprb

  if (.not. do_merge_only) then
    allocate(flux_dn(nspec,nlev+1))
    allocate(flux_up(nspec,nlev+1))

    allocate(bb_flux_dn(nlev+1))
    allocate(bb_flux_up(nlev+1))

    allocate(planck_hl(nspec,nlev+1))
  end if

  allocate(mole_fraction_fl(nlev,ngas))
  allocate(mole_fraction_scaling(nlev))

  if (iendcol <= 0) then
    iendcol = ncol
  else
    iendcol = iendcol
  end if
  ncolout = iendcol + 1 - istartcol

  if (do_merge_only) then
    config%do_write_optical_depth = .true.
    config%do_write_spectral_fluxes = .false.
    config%do_write_planck = .false.
  end if

  if (config%do_write_planck &
       &  .or. config%do_write_optical_depth &
       &  .or. config%do_write_spectral_fluxes &
       &  .or. config%do_write_spectral_boundary_fluxes) then
    do_save_spectrum = .true.
  else
    do_save_spectrum = .false.
  end if

  ! Define dimensions of output file
   if (config%do_write_planck &
       &  .or. config%do_write_optical_depth &
       &  .or. config%do_write_spectral_fluxes) then
    ! File will be large: don't use an unlimited dimension, enabling
    ! compression
    call out_file%define_dimension("column", ncolout)
  else
    ! File will be small: use an unlimited dimension so that files can
    ! be easily concatenated along this dimension
    call out_file%define_dimension("column", 0)
  end if

  call out_file%define_dimension("gas", ngas)
  call out_file%define_dimension("level", nlev)
  call out_file%define_dimension("half_level", nlev+1)
  if (nbands > 0) then
    call out_file%define_dimension("band_lw", nbands)
  end if
  if (do_save_spectrum) then
    call out_file%define_dimension(trim(spectral_dim_name), nspec)
  end if

  if (.not. do_merge_only) then
    call out_file%define_variable("radiative_transfer_mode", &
         &  long_name="Radiative transfer mode", &
         &  data_type_name='short', &
         &  comment_str='If zero, a single angle of 53 degrees was used per hemisphere (diffusivity factor of 1.66);' &
         &  // new_line('a') // 'otherwise, this is the number of angles per hemisphere using Gauss-Legendre quadrature.')
  end if

  ! Define coordinate variables
  call out_file%define_variable("pressure_hl", &
       &   dim2_name="column", dim1_name="half_level", &
       &   units_str="Pa", long_name="Pressure", &
       &   standard_name="air_pressure")
  if (.not. config%do_read_planck) then
    call out_file%define_variable("temperature_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="K", long_name="Temperature", &
         &   standard_name="air_temperature")
  end if

  if (.not. is_ckd) then
    call out_file%define_variable("reference_surface_mole_fraction", &
         &  dim1_name="gas", &
         &  units_str="1", long_name="Reference surface mole fraction", &
         &  fill_value=-1.0_jprb)
  end if

  ! Define other
  if (.not. is_ckd) then
    call out_file%define_variable("mole_fraction_fl", &
         &  dim3_name="column", dim2_name="gas", dim1_name="level", &
         &  units_str="1", long_name="Mole fraction")
  end if

  if (.not. do_merge_only) then
    call out_file%define_variable("flux_up_lw", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str=LW_UNITS_STR, long_name="Upwelling longwave flux", &
         &   standard_name="upwelling_longwave_flux_in_air")
    call out_file%define_variable("flux_dn_lw", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str=LW_UNITS_STR, long_name="Downwelling longwave flux", &
         &   standard_name="downwelling_longwave_flux_in_air")
  end if

  if (nbands > 0) then
    call out_file%define_variable("band_wavenumber1_lw", &
         &   dim1_name="band_lw", long_name="Lower bound wavenumber for longwave band", &
         &   units_str="cm-1")
    call out_file%define_variable("band_wavenumber2_lw", &
         &   dim1_name="band_lw", long_name="Upper bound wavenumber for longwave band", &
         &   units_str="cm-1")
    call out_file%define_variable("band_flux_up_lw", &
         &   dim3_name="column", dim2_name="half_level", dim1_name="band_lw", &
         &   units_str=LW_UNITS_STR, long_name="Upwelling longwave flux in bands")
    call out_file%define_variable("band_flux_dn_lw", &
         &   dim3_name="column", dim2_name="half_level", dim1_name="band_lw", &
         &   units_str=LW_UNITS_STR, long_name="Downwelling longwave flux in bands")
  end if

  if (do_save_spectrum .and. .not. is_ckd) then
    call out_file%define_variable(trim(spectral_dim_name), &
         &   dim1_name=trim(spectral_dim_name), units_str="cm-1", &
         &   long_name=trim(spectral_dim_name), is_double=.true., &
         &   deflate_level=2, shuffle=.true.)
  end if

  if (config%do_write_optical_depth) then
    call out_file%define_variable("optical_depth", &
         &   dim3_name="column", dim2_name="level", dim1_name=trim(spectral_dim_name), &
         &   long_name="Layer optical depth")
  end if
  if (config%do_write_planck) then
    call out_file%define_variable("planck_hl", &
         &   dim3_name="column", dim2_name="half_level", dim1_name=trim(spectral_dim_name), &
         &   units_str=LW_UNITS_STR, long_name="Planck function")
  end if
  if (config%do_write_spectral_fluxes) then
    call out_file%define_variable("spectral_flux_up_lw", &
         &   dim3_name="column", dim2_name="half_level", dim1_name=trim(spectral_dim_name), &
         &   units_str=LW_UNITS_STR, long_name="Upwelling spectral longwave flux", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
    call out_file%define_variable("spectral_flux_dn_lw", &
         &   dim3_name="column", dim2_name="half_level", dim1_name=trim(spectral_dim_name), &
         &   units_str=LW_UNITS_STR, long_name="Downwelling spectral longwave flux", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
  end if


  if (config%do_write_spectral_boundary_fluxes) then
    call out_file%define_variable("spectral_flux_up_toa_lw", &
         &   dim2_name="column", dim1_name=trim(spectral_dim_name), &
         &   units_str=LW_UNITS_STR, long_name="Upwelling spectral longwave flux at top-of-atmosphere", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
    call out_file%define_variable("spectral_flux_dn_surf_lw", &
         &   dim2_name="column", dim1_name=trim(spectral_dim_name), &
         &   units_str=LW_UNITS_STR, long_name="Downwelling spectral longwave flux at surface", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
  end if

  ! Write global attributes
  if (is_ckd) then
    call out_file%put_global_attributes(title_str="Longwave fluxes from a CKD model", &
         &  project_str="CKDMIP", conventions_str="CF-1.7")
  else
    if (do_merge_only) then
      call out_file%put_global_attributes(title_str="Spectral optical depth profiles from multiple gases", &
           &  project_str="CKDMIP", conventions_str="CF-1.7")
      call out_file%put_global_attribute("profile_type", "merge")
    else
      call out_file%put_global_attributes(title_str="Line-by-line longwave fluxes", &
           &  project_str="CKDMIP", conventions_str="CF-1.7")
    end if
    call out_file%put_global_attribute("constituent_id", molecule_list)
  end if
  if (is_scenario) then
    call out_file%put_global_attribute("scenario", trim(scenario_str))
  end if
  call out_file%put_global_attribute("software_version", CKDMIP_VERSION)

  if (.not. do_merge_only) then
    call out_file%put("radiative_transfer_mode", real(min(MAX_ANGLES, config%nangle),jprb))
  end if

  ! Write coordinate variables
  call out_file%put("pressure_hl", pressure_hl(:,istartcol:iendcol))
  if (.not. config%do_read_planck) then
    call out_file%put("temperature_hl", temperature_hl(:,istartcol:iendcol))
  end if

  if (nbands > 0) then
    call out_file%put("band_wavenumber1_lw", band_wavenumber1(1:nbands))
    call out_file%put("band_wavenumber2_lw", band_wavenumber2(1:nbands))
  end if

  if (do_save_spectrum .and. .not. is_ckd) then
    call out_file%put(trim(spectral_dim_name), wavenumber_cm1(1:nspec))
  end if

  if (.not. is_ckd) then
    call out_file%put("reference_surface_mole_fraction", target_conc(1:ngas))
  end if

  allocate(od(nspec,nlev))

  do jcol = istartcol,iendcol

    if (config%iverbose >= 3) then
      write(*,'(a,i0)') 'Column ', jcol
    end if

    ! Read and sum optical depths
    od = 0.0_jprb

    if (is_ckd) then
      ! Read from a single CKD file
      call in_file(1)%get(trim(config%planck_name), planck_hl, jcol)

      if (config%input_planck_per_sterad) then
        if (config%iverbose >= 3) then
          write(*,'(a)') '    Scaling input Planck function by pi'
        end if
        planck_hl = planck_hl * PI
      else
        if (config%iverbose >= 3) then
          write(*,'(a)') '    Assuming that Planck function is power into a horizontal surface'
        end if   
      end if

      if (in_file(1)%exists(trim(config%surf_emissivity_name))) then
        call in_file(1)%get(trim(config%surf_emissivity_name), surf_emissivity, jcol)
      else
        surf_emissivity = 1.0_jprb
      end if

      if (in_file(1)%exists(trim(config%surf_emission_name))) then
        call in_file(1)%get(trim(config%surf_emission_name), surf_emission, jcol)
      else
        surf_emission = planck_hl(:,nlev+1) * surf_emissivity
      end if

      call in_file(1)%get(trim(config%optical_depth_name), od, jcol)

      call calc_longwave_fluxes_n(config%nangle, nlev, nspec, 1, nspec, od, planck_hl, &
           &  surf_emissivity, surf_emission, flux_dn, flux_up)

    else
      ! Read from multiple gas files
      igas = 1
      do jfile = 1,nfile
        mole_fraction_scaling = 1.0_jprb
        if (ngas_in_file(jfile) == 1) then
          ! Single gas
          ! First the mole fraction
          call in_file(jfile)%get("mole_fraction_fl", mole_fraction_fl_1gas, jcol)
          if (size(mole_fraction_fl_1gas) /= nlev) then
            error stop 'Error: incorrect size for mole_fraction_fl'
          end if
          if (conc_operator(jfile) == I_SCALE) then
            mole_fraction_scaling = conc_factor(jfile)
            if (config%iverbose >= 3) then
              write(*,'(a,i0,a,f8.3)') '  Reading gas ', igas, ': scaling by ', conc_factor(jfile)
            end if
          else if (conc_operator(jfile) == I_CONST) then
            mole_fraction_scaling = conc_factor(jfile) / mole_fraction_fl_1gas
            if (config%iverbose >= 3) then
              write(*,'(a,i0,a)') '  Reading gas ', igas, ': scaling by variable amount with pressure'
            end if
          else
            if (config%iverbose >= 3) then
              write(*,'(a,i0)') '  Reading gas ', igas
            end if
          end if
          mole_fraction_fl_1gas = mole_fraction_fl_1gas * mole_fraction_scaling
          mole_fraction_fl(:,igas) = mole_fraction_fl_1gas
        else
          ! Multiple gases
          ! First the mole fraction
          call in_file(jfile)%get("mole_fraction_fl", mole_fraction_fl_ngas, jcol)
          if (config%iverbose >= 3) then
            write(*,'(a,i0,a,i0)') '  Reading gases ', igas, '-', igas+ngas_in_file(jfile)-1
          end if
          mole_fraction_fl(:,igas:igas+ngas_in_file(jfile)-1) = mole_fraction_fl_ngas
        end if

        call in_file(jfile)%get(trim(config%optical_depth_name), od_1gas, jcol)

        write(*,*) '  raw optical depth range: ', minval(od_1gas), '-', maxval(od_1gas)

        if (size(od,2) /= nlev) then
          error stop 'Error: incorrect size for optical depth'
        end if

        do jlev = 1,nlev
          if (mole_fraction_scaling(jlev) /= 1.0_jprb) then
            od_1gas(:,jlev) = od_1gas(:,jlev) * mole_fraction_scaling(jlev)
          end if
        end do

        write(*,*) '  optical depth range: ', minval(od_1gas), '-', maxval(od_1gas)

        if (config%nspectralstride /= 1) then
          od = od + od_1gas(1:ninspec:config%nspectralstride,:)
        else
          od = od + od_1gas
        end if

        igas = igas + ngas_in_file(jfile)
      end do

      deallocate(od_1gas)
      if (allocated(mole_fraction_fl_1gas))  deallocate(mole_fraction_fl_1gas)
      if (allocated(mole_fraction_fl_ngas))  deallocate(mole_fraction_fl_ngas)

      if (config%do_read_planck) then
        call in_file(1)%get(trim(config%planck_name), planck_hl, jcol)
        call in_file(1)%get(trim(config%surf_emission_name), surf_emission, jcol)
        call in_file(1)%get(trim(config%surf_emissivity_name), surf_emissivity, jcol)
      end if

      nblock = (nspec - 1 + config%nblocksize) / config%nblocksize

      if (.not. do_merge_only) then

        if (config%iverbose >= 3) then
          write(*,'(a)') '  Radiative transfer'
        end if

        !$OMP PARALLEL DO PRIVATE(istartspec,iendspec) SCHEDULE(RUNTIME)
        do jblock = 1,nblock
          istartspec = (jblock-1) * config%nblocksize + 1
          iendspec = min(istartspec + config%nblocksize - 1, nspec)

          if (.not. config%do_read_planck) then
            call calc_planck_function(nspec,nlev+1,istartspec,iendspec, &
                 &  temperature_hl(:,jcol), wavenumber_cm1, planck_hl)
            if (config%surf_temperature >= 0.0_jprb) then
              call calc_planck_function1(nspec,istartspec,iendspec, &
                   &  config%surf_temperature, wavenumber_cm1, surf_emission)
              surf_emission(istartspec:iendspec) = surf_emission(istartspec:iendspec) &
                   &  * surf_emissivity(istartspec:iendspec)
            else if (allocated(surf_temperature)) then
              call calc_planck_function1(nspec,istartspec,iendspec, &
                   &  surf_temperature(jcol), wavenumber_cm1, surf_emission)
              surf_emission(istartspec:iendspec) = surf_emission(istartspec:iendspec) &
                   &  * surf_emissivity(istartspec:iendspec)
            else
              surf_emission(istartspec:iendspec) &
                   &  = surf_emissivity(istartspec:iendspec)*planck_hl(istartspec:iendspec,nlev+1)
            end if
          end if

          call calc_longwave_fluxes_n(config%nangle, nlev, nspec, istartspec, iendspec, &
               &  od,planck_hl, surf_emissivity, surf_emission, flux_dn, flux_up)
        end do
        !$OMP END PARALLEL DO

      end if

    end if

    icolout = jcol + 1 - istartcol

    if (.not. is_ckd) then
      ! Write a slice of mole fraction to this column
      call out_file%put("mole_fraction_fl", mole_fraction_fl, icolout)
    end if

    if (.not. do_merge_only) then

      bb_flux_dn = sum(flux_dn, 1)
      bb_flux_up = sum(flux_up, 1)

      ! Write slice of broadband flux corresponding to this column
      call out_file%put("flux_up_lw", bb_flux_up, icolout)
      call out_file%put("flux_dn_lw", bb_flux_dn, icolout)

      ! Write band fluxes
      if (nbands > 0) then
        do jband = 1,nbands
          do jlev = 1,nlev+1
            band_flux_dn(jband,jlev) = sum(flux_dn(:,jlev), mask=(band_number==jband))
            band_flux_up(jband,jlev) = sum(flux_up(:,jlev), mask=(band_number==jband))
          end do
        end do
        call out_file%put("band_flux_up_lw", band_flux_up, icolout)
        call out_file%put("band_flux_dn_lw", band_flux_dn, icolout)
      end if

    end if

    ! Write optical depth
    if (config%do_write_optical_depth) then
      call out_file%put("optical_depth", od, icolout)
    end if

    ! Write slice of spectral flux corresponding to this column
    if (config%do_write_planck) then
      call out_file%put("planck_hl", planck_hl, icolout)
    end if
    if (config%do_write_spectral_fluxes) then
      call out_file%put("spectral_flux_up_lw", flux_up, icolout)
      call out_file%put("spectral_flux_dn_lw", flux_dn, icolout)
    end if

    if (config%do_write_spectral_boundary_fluxes) then
      call out_file%put("spectral_flux_up_toa_lw", &
           &  flux_up(:,1), icolout)
      call out_file%put("spectral_flux_dn_surf_lw", &
           &  flux_dn(:,nlev), icolout)
    end if

  end do
    
  do jfile = 1,nfile
    call in_file(jfile)%close()
  end do

  call out_file%close()

contains

  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: ckdmip_lw [arguments]'
    write(*,'(a)') '  where the possible arguments are:'
    write(*,'(a)') '              input.h5  : Add unscaled spectral optical depth from this file'
    write(*,'(a)') '    --scale X input.h5  : Add spectral optical depths scaled by X'
    write(*,'(a)') '    --conc  X input.h5  : Add optical depths scaled so surface mole fraction is X'
    write(*,'(a)') '    --const X input.h5  : Add optical depths with constant mole fraction of X'
    write(*,'(a)') '    --ckd     input.h5  : Read optical depth and Planck function from CKD model'
    write(*,'(a)') '           Note: this option is not compatible with spectral optical depth inputs'
    write(*,'(a)') ' -o|--output  output.h5 : Output file'
    write(*,'(a)') ' -c|--config  conf.nam  : Configure using namelist file'
    write(*,'(a)') ' -m|--merge-only        : Only merge optical depths'
    write(*,'(a)') '    --scenario str      : Add a "scenario" global attribute to the output'
    write(*,'(a)') '    --column-range M N  : Only process columns M to N'
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
  
end program ckdmip_lw
