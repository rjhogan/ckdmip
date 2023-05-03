! ckdmip_sw.F90 - shortwave radiative transfer for CKDMIP
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

program ckdmip_sw

  use parkind1,         only : jprb
  use easy_netcdf,      only : netcdf_file
  use shortwave_config, only : shortwave_config_type
  use shortwave_fluxes, only : calc_shortwave_fluxes

  implicit none

#include "ckdmip_version.h"

  integer,      parameter :: N_MAX_FILES  = 32
  character(5), parameter :: SW_UNITS_STR = 'W m-2'

  integer, parameter :: I_AS_IS = 0
  integer, parameter :: I_SCALE = 1
  integer, parameter :: I_CONC  = 2
  integer, parameter :: I_CONST = 3

  real, parameter :: LAYER_OD_TOL = 1.0e-15

  ! How to treat concentrations (0-4 above)
  integer :: conc_operator(N_MAX_FILES)

  ! Concentration factor, interpreted according to conc_operator
  real(jprb) :: conc_factor(N_MAX_FILES)

  ! Input files
  type(netcdf_file) :: in_file(N_MAX_FILES)

  ! Solar spectral irradiance file
  type(netcdf_file) :: ssi_file

  ! Reference surface concentration
  real(jprb) :: ref_conc(N_MAX_FILES)

  ! Target surface concentration
  real(jprb) :: target_conc(N_MAX_FILES)

  type(netcdf_file) :: out_file

  type(shortwave_config_type) :: config

  ! Number of levels, spectral intervals, gases and profiles
  integer :: nlev, ngas, ncol, nspec

  ! Number of levels in input files, may be modified
  integer :: nlevorig

  ! Number of blocks into which to divide the calculation
  integer :: nblock

  ! Loop indices
  integer :: jlev, jspec, jgas, jcol, jblock, jarg, jband, jmu0

  ! Range of spectral indices to process
  integer :: istartspec,iendspec

  ! Range of spectral indices to read from file
  integer :: istartspecin, iendspecin

  ! Range of wavenumbers to read from file (cm-1)
  real(jprb) :: wn_start, wn_end

  ! Range of columns to process
  integer :: istartcol, iendcol

  ! Number and index to output column
  integer ::ncolout, icolout, iendcolout
  
  ! Pressure and temperature at half levels (nlev+1,ncol)
  real(jprb), allocatable, dimension(:,:) :: pressure_hl, temperature_hl
  ! Surface albedo and incoming solar spectral irradiance (W m-2)
  real(jprb), allocatable, dimension(:) :: surf_albedo, ssi, ssi_tmp
  ! Wavenumber in cm-1 (nspec)
  real(jprb), allocatable, dimension(:)   :: wavenumber_cm1
  ! Optical depth of 1 gas (nspecin,nlev) and all gases (nspec,nlev)
  real(jprb), allocatable, dimension(:,:) :: od_1gas
  real(jprb), allocatable, dimension(:,:) :: od, ssa, asymmetry

  ! Other optical properties of single species
  real(jprb) :: ssa_scalar, asymmetry_scalar
  real(jprb), allocatable, dimension(:,:) :: ssa_1gas, asymmetry_1gas

  ! Fluxes (W m-2)
  real(jprb), allocatable, dimension(:,:) :: flux_dn, flux_up, flux_dn_direct

  ! Broadband fluxes (W m-2)
  real(jprb), allocatable, dimension(:) :: bb_flux_dn, bb_flux_up, bb_flux_dn_direct

  ! Mole fraction (mol mol-1), dimensioned (ngas,nlev)
  real(jprb), allocatable, dimension(:,:) :: mole_fraction_fl

  ! Mole fraction (mol mol-1) of a single gas (nlev)
  real(jprb), allocatable, dimension(:) :: mole_fraction_fl_1gas

  ! Mole fraction (mol mol-1) of a single gas on half levels (nlev+1)
  real(jprb), allocatable, dimension(:) :: mole_fraction_hl_1gas

  ! Concentration scaling as a function of pressure
  real(jprb), allocatable, dimension(:) :: mole_fraction_scaling

  ! Band boundaries
  real(jprb), allocatable, dimension(:) :: band_wavenumber1, band_wavenumber2

  ! Number of band to which each wavenumber belongs
  integer, allocatable, dimension(:) :: band_number

  ! Fluxes in bands (W m-2), dimensioned (nband, nlev+1)
  real(jprb), allocatable, dimension(:,:) :: band_flux_dn, band_flux_up, band_flux_dn_direct

  ! Number and index of command-line arguments, and number of bands
  integer :: narg, iarg, nbands

  integer :: istatus, itoken

  ! Optionally split layers between isplit1 and isplit2 (inclusive)
  ! into nsplit sub-layers
  integer :: nsplit, isplit1, isplit2

  ! File name
  character(len=512) :: argument_str, file_name, out_file_name

  character(len=512) :: molecule_str, molecule_list, scenario_str
  character(len=32)  :: spectral_dim_name, flux_dim3_name, flux_dim2_name
  character(len=32)  :: spectral_level_dim_name

  ! Spectral fluxes can be written on a subset of the full
  ! half-levels, in which case we use
  ! config%i_spectral_level_index(1:nspeclev) as indices to the full
  ! spectral flux profile
  integer :: nspeclev

  integer :: nfluxdims

  ! Do we output anything a function of wavenumber or g-point?
  logical :: do_save_spectrum

  ! Do we simply merge optical depths without any radiative transfer?
  logical :: do_merge_only

  ! Is there a scenario string?
  logical :: is_scenario

  ! Is the input in the form of optical depths in the g-points of a
  ! CKD model?
  logical :: is_ckd

  ! Are we only processing part of the spectrum?
  logical :: do_spectral_subsample

  ! Do we have any scatterers to consider
  logical :: have_scatterers

  ! Initialize configuration variables
  do_merge_only = .false.
  is_scenario   = .false.
  is_ckd        = .false.
  istartcol     = 1
  iendcol       = 0
  nbands        = 0
  istartspecin  = 1
  iendspecin    = 0
  do_spectral_subsample = .false.
  have_scatterers = .false.
  wn_start      = -1.0
  wn_end        = -1.0

  ! Initialize concentration scaling variables
  conc_operator = I_AS_IS   ! No scaling by default
  conc_factor   = -1.0_jprb ! No scaling factor by default
  ref_conc      = -1.0_jprb 
  target_conc   = -1.0_jprb

  ! Initialize gas and argument counters
  ngas = 0
  iarg = 1

  ! No splitting of layers
  nsplit = 1

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
      config%nmu0 = 1

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

    else if (argument_str == '--spectral-range') then

      if (iarg >= narg-1) then
        write(*,'(a)') '"--spectral-range" must be followed by two arguments'
        call print_usage_and_exit()
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) istartspecin

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) iendspecin

      if (istartspecin < 1 .or. iendspecin < istartspecin) then
        write(*,'(a,i0,a,i0,a)') 'Error: spectral range of ', istartspecin, ' to ', iendspecin, ' is invalid'
        error stop
      end if

      if (config%iverbose >= 3) then
        write(*,'(a,i0,a,i0)') 'Reading spectral points ', istartspecin, ' to ', iendspecin
      end if

    else if (argument_str == '--wavenumber-range') then

      if (iarg >= narg-1) then
        write(*,'(a)') '"--wavenumber-range" must be followed by two arguments'
        call print_usage_and_exit()
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) wn_start

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) wn_end

      if (wn_start < 0.0_jprb .or. wn_end < wn_start) then
        write(*,'(a,f0.1,a,f0.1,a)') 'Error: wavenumber range of ', wn_start, ' to ', wn_end, ' is invalid'
        error stop
      end if

      if (config%iverbose >= 3) then
        write(*,'(a,f0.1,a,f0.1,a)') 'Reading wavenumbers in the range [', wn_start, ', ', wn_end, ')'
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

      ngas = ngas+1
      if (argument_str == '--scale') then
        conc_operator(ngas) = I_SCALE
      else if (argument_str == '--const') then
        conc_operator(ngas) = I_CONST
      else if (argument_str == '--conc') then
        conc_operator(ngas) = I_CONC
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) conc_factor(ngas)

      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of input file as string of length < 512'
        call print_usage_and_exit()
      end if

      call in_file(ngas)%open(trim(file_name), iverbose=config%iverbose)

    else if (argument_str == '--ckd') then
      ! Read a Correlated K-Distribution file containing optical depth
      ! in a number of g-points
      if (ngas /= 0) then
        write(*,'(a)') 'Error: cannot read CKD optical depth file if already using spectral optical depth file(s)'
        error stop
      else if (iarg == narg) then
        write(*,'(a)') 'Error: "--ckd" must be followed by a file name'
        error stop
      end if

      is_ckd = .true.

      ngas = ngas+1
      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of CKD input file as string of length < 512'
        call print_usage_and_exit()
      end if
      call in_file(ngas)%open(trim(file_name), iverbose=config%iverbose)

    else if (argument_str == '--layer-split') then
      if (iarg >= narg-1) then
        write(*,'(a)') 'Error: "--layer-split" must be followed by three arguments'
        error stop
      end if
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) isplit1
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) isplit2
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) nsplit

    else if (argument_str == '--ssi') then

      if (iarg == narg) then
        write(*,'(a)') 'Error: "--ssi" must be followed by a file name'
        error stop
      end if

      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of SSI input file as string of length < 512'
        call print_usage_and_exit()
      end if
      call ssi_file%open(trim(file_name), iverbose=config%iverbose)

      call ssi_file%get("solar_spectral_irradiance", ssi);

      call ssi_file%close()

    else if (argument_str == '--cos-sza') then

      if (iarg == narg) then
        write(*,'(a)') 'Error: "--cos-sza" must be followed by a number'
        error stop
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) config%cos_solar_zenith_angle(1)
      config%nmu0 = 1
      
    else if (argument_str == '--sza') then

      if (iarg == narg) then
        write(*,'(a)') 'Error: "--sza" must be followed by the solar zenith angle in degrees'
        error stop
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) config%cos_solar_zenith_angle(1)
      config%cos_solar_zenith_angle(1) &
           &  = cos(config%cos_solar_zenith_angle(1) &
           &        * 3.14159265358979323846_jprb / 180.0_jprb )
      config%nmu0 = 1
      
    else if (argument_str(1:1) == '-') then
      
      write(*,'(a,a,a)') 'Argument "', trim(argument_str), '" not understood'
      call print_usage_and_exit()

    else
      ! Open file containing the spectral optical depth of a single
      ! gas
      if (is_ckd) then
        write(*,'(a)') 'Error: cannot read spectral optical depth files if already using CKD optical depth file'
        error stop
      end if

      ngas = ngas+1

      file_name = trim(argument_str)

      call in_file(ngas)%open(trim(file_name), iverbose=config%iverbose)

      conc_operator(ngas) = I_AS_IS
      conc_factor(ngas)   = 1.0_jprb

    end if

    iarg = iarg + 1

  end do

  if (ngas == 0 .or. narg == 0) then
    call print_usage_and_exit()
  end if

  if (.not. out_file%is_open()) then
    error stop 'Output file not specified'
  end if

  if (.not. do_merge_only) then
    if (config%nmu0 < 1) then
      error stop 'Cosine of solar zenith angle not specified in namelist file or on command line'
    else if (any(config%cos_solar_zenith_angle(1:config%nmu0) < -1.0_jprb &
         &       .or. config%cos_solar_zenith_angle(1:config%nmu0) > 1.0_jprb)) then
      error stop 'Cosine of the solar zenith angle not specified or outside the valid range of -1 to 1'
    end if
  end if

  where (config%cos_solar_zenith_angle < 0.0_jprb)
    config%cos_solar_zenith_angle = 0.0_jprb
  end where

  if (.not. (is_ckd .or. do_merge_only) .and. .not. allocated(ssi)) then
    error stop 'Solar spectral irradiance file not specified on the command line';
  end if

  if (is_ckd .and. do_merge_only) then
    error stop 'Error: "--ckd" and "--merge-only" arguments are incompatible'
  end if

  if (is_ckd) then
    spectral_dim_name = 'gpoint_sw'
  else
    spectral_dim_name = 'wavenumber'    
  end if

  if (config%nspectralstride > 1 .or. istartspecin > 1 .or. iendspecin > 0 &
       &  .or. wn_start >= 0.0_jprb .or. wn_end > 0.0_jprb) then
    do_spectral_subsample = .true.
  end if

  if (.not. is_ckd) then
    ! Report what gases will be merged

    do jgas = 1,ngas
      if (in_file(jgas)%global_attribute_exists("constituent_id")) then
        call in_file(jgas)%get_global_attribute("constituent_id", molecule_str)
      else
        write(*,'(a,i0)') 'Error: constituent_id not present for gas ', jgas
        error stop
      end if
      if (jgas == 1) then
        molecule_list = molecule_str
      else
        molecule_list = trim(molecule_list) // ' ' // trim(molecule_str)
      end if
      
      ! Load reference concentration, if present
      if (in_file(jgas)%exists("reference_surface_mole_fraction")) then
        call in_file(jgas)%get("reference_surface_mole_fraction", ref_conc(jgas))
        target_conc(jgas) = ref_conc(jgas)
      else if (conc_operator(jgas) == I_CONC .or. conc_operator(jgas) == I_CONST) then
        write(*,'(a,a,a)') 'Error: cannot set concentration of ', to_upper(trim(molecule_str)), &
             &         '  as no reference surface concentration is present'
        error stop
      end if
      
      if (config%iverbose >= 3) then
        write(*,'(a,i0,a,a)') '  Gas ', jgas, ': ', to_upper(trim(molecule_str))
        if (conc_operator(jgas) == I_AS_IS) then
          write(*,'(a)') '    No scaling to be applied'
        else if (conc_operator(jgas) == I_SCALE) then
          write(*,'(a,e11.5)') '    Scaling by a factor of ', conc_factor(jgas)
        else if (conc_operator(jgas) == I_CONC) then
          write(*,'(a,e11.5,a,e11.5,a)') '    Setting surface mole fraction to ', conc_factor(jgas), ' mol mol-1 ' &
               &             //'while reference value is ', ref_conc(jgas), ' mol mol-1:'
          write(*,'(a,f8.3)') '    Scaling by a factor of ', conc_factor(jgas) / ref_conc(jgas)
        else if (conc_operator(jgas) == I_CONST) then
          write(*,'(a,e11.5,a)') '    Setting mole fraction to a constant value of ', conc_factor(jgas), ' mol mol-1'
        end if
      end if
      
      ! If a concentration has been specified then change it to a
      ! scaling
      if (conc_operator(jgas) == I_CONC) then
        conc_factor(jgas) = conc_factor(jgas) / ref_conc(jgas)
        conc_operator(jgas) = I_SCALE
      end if
      if (conc_operator(jgas) == I_SCALE) then
        target_conc(jgas) = ref_conc(jgas) * conc_factor(jgas)
      else if (conc_operator(jgas) == I_CONST) then
        target_conc(jgas) = conc_factor(jgas)
      end if
      
    end do
    
  end if

  ! Read basic properties from the first input file

  call in_file(1)%get(trim(config%pressure_name), pressure_hl)
  if (config%pressure_scaling /= 1.0_jprb) then
    pressure_hl = pressure_hl * config%pressure_scaling
  end if

  ncol     = size(pressure_hl,2)
  nlevorig = size(pressure_hl,1)-1

  if (nsplit > 1) then
    call split_array_hl(isplit1, isplit2, nsplit, pressure_hl)
  end if

  if (nsplit > 1) then
    nlev = nlevorig + (nsplit-1)*(isplit2-isplit1+1)
  else
    nlev = nlevorig
  end if

  ! Don't output bands if we are only merging spectra
  if (do_merge_only) then
    nbands = 0
    config%use_mu0_dimension = .false.
  end if

  if (.not. is_ckd) then
    ! Spectral files
    call in_file(1)%get(trim(config%wavenumber_name), wavenumber_cm1)
    if (.not. do_spectral_subsample) then
      nspec = size(wavenumber_cm1)
      iendspecin = nspec
    else
      ! Are we selecting by wavenumber or index?  Wavenumber takes precedence
      if (wn_start >= 0.0_jprb .or. wn_end > 0.0_jprb) then
        istartspecin = minloc(wavenumber_cm1, 1, mask=(wavenumber_cm1 >= wn_start))
        iendspecin   = maxloc(wavenumber_cm1, 1, mask=(wavenumber_cm1 < wn_end))
      end if

      if (iendspecin < 1) then
        iendspecin = size(wavenumber_cm1)
      end if
      nspec = (iendspecin-istartspecin+config%nspectralstride) / config%nspectralstride
      iendspecin = istartspecin+(nspec-1)*config%nspectralstride
      write(*,'(a,i0,a,i0,a,i0,a,i0)') 'Subsampling ', nspec, ' wavenumbers ', istartspecin, '-', &
           &  iendspecin, ' with stride ', config%nspectralstride
      do jspec = 1,nspec
        wavenumber_cm1(jspec) = wavenumber_cm1(istartspecin+(jspec-1)*config%nspectralstride)
      end do

      ! Average spectral solar irradiance
      do jspec = 1,nspec
        ssi(jspec) = sum(ssi(istartspecin+(jspec-1)*config%nspectralstride:istartspecin-1+jspec*config%nspectralstride))
      end do

    end if
      
    call in_file(1)%get(trim(config%temperature_name), temperature_hl)
    if (nsplit > 1) then
      call split_array_hl(isplit1, isplit2, nsplit, temperature_hl)
    end if

    ! Find locations of the requested bands
    if (nbands > 1) then
      allocate(band_number(nspec))
      band_number = 0
      allocate(band_flux_dn(nbands,nlev+1))
      allocate(band_flux_dn_direct(nbands,nlev+1))
      allocate(band_flux_up(nbands,nlev+1))
      do jband = 1,nbands
        where (wavenumber_cm1 > band_wavenumber1(jband) .and. wavenumber_cm1 <= band_wavenumber2(jband))
          band_number = jband
        end where
      end do
    end if

  else
    ! CKD file
    call in_file(1)%get(trim(config%optical_depth_name), od, 1)
    nspec = size(od,1)
    deallocate(od);
    !ninspec = nspec
  end if

  allocate(surf_albedo(nspec))
  surf_albedo = config%surf_albedo ! 0.06 is the default

  if (.not. do_merge_only) then
    allocate(flux_dn(nspec,nlev+1))
    allocate(flux_dn_direct(nspec,nlev+1))
    allocate(flux_up(nspec,nlev+1))

    allocate(bb_flux_dn(nlev+1))
    allocate(bb_flux_dn_direct(nlev+1))
    allocate(bb_flux_up(nlev+1))
  end if

  allocate(mole_fraction_scaling(nlevorig))

  if (iendcol <= 0) then
    iendcol = ncol
  else
    iendcol = iendcol
  end if
  ncolout = iendcol + 1 - istartcol

  if (do_merge_only) then
    config%do_write_optical_depth = .true.
    config%do_write_spectral_fluxes = .false.
  end if

  if (config%do_write_optical_depth &
       &  .or. config%do_write_spectral_fluxes &
       &  .or. config%do_write_spectral_boundary_fluxes) then
    do_save_spectrum = .true.
  else
    do_save_spectrum = .false.
  end if

  ! Define dimensions of output file
 
  if (config%use_mu0_dimension) then
    flux_dim3_name = "column"
    flux_dim2_name = "mu0"
    nfluxdims = 3
  else
    flux_dim3_name = " "
    flux_dim2_name = "column"
    nfluxdims = 2
  end if

  if (config%do_write_optical_depth &
       &  .or. config%do_write_spectral_fluxes) then
    ! File will be large: don't use an unlimited dimension, enabling
    ! compression
    if (config%use_mu0_dimension) then
      call out_file%define_dimension("column", ncolout)
      call out_file%define_dimension("mu0", config%nmu0)
    else
      call out_file%define_dimension("column", ncolout*config%nmu0)
    end if
  else
    ! File will be small: use an unlimited dimension so that files can
    ! be easily concatenated along this dimension
    call out_file%define_dimension("column", 0)
    if (config%use_mu0_dimension) then
      call out_file%define_dimension("mu0", config%nmu0)
    end if
  end if

  call out_file%define_dimension("gas", ngas)
  call out_file%define_dimension("level", nlev)
  call out_file%define_dimension("half_level", nlev+1)
  if (nbands > 0) then
    call out_file%define_dimension("band_sw", nbands)
  end if
  if (do_save_spectrum) then
    call out_file%define_dimension(trim(spectral_dim_name), nspec)
  end if

  ! Define coordinate variables
  call out_file%define_variable("pressure_hl", &
       &   dim2_name="column", dim1_name="half_level", &
       &   units_str="Pa", long_name="Pressure", &
       &   standard_name="air_pressure")
  if (allocated(temperature_hl)) then
    call out_file%define_variable("temperature_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="K", long_name="Temperature", &
         &   standard_name="air_temperature")
  end if
  if (.not. do_merge_only) then
    if (config%use_mu0_dimension) then
      call out_file%define_variable("mu0", &
           &  dim1_name="mu0", &
           &  long_name="Cosine of solar zenith angle")
    else
      call out_file%define_variable("mu0", &
           &  dim1_name="column", &
           &  long_name="Cosine of solar zenith angle")
    end if
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
         &  units_str="1", long_name="Mole fraction", fill_value=-1.0_jprb)
  end if

  if (.not. do_merge_only) then
    if (.not. config%do_write_direct_only) then
      call out_file%define_variable("flux_up_sw", &
           &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims, dim1_name="half_level", &
           &   units_str=SW_UNITS_STR, long_name="Upwelling shortwave flux", &
           &   standard_name="upwelling_shortwave_flux_in_air")
      call out_file%define_variable("flux_dn_sw", &
           &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims, dim1_name="half_level", &
           &   units_str=SW_UNITS_STR, long_name="Downwelling shortwave flux", &
           &   standard_name="downwelling_shortwave_flux_in_air")
    end if
    call out_file%define_variable("flux_dn_direct_sw", &
         &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
         &   ndims=nfluxdims, dim1_name="half_level", &
         &   units_str=SW_UNITS_STR, long_name="Direct downwelling shortwave flux")
  end if

  if (nbands > 0) then
    call out_file%define_variable("band_wavenumber1_sw", &
         &   dim1_name="band_sw", long_name="Lower bound wavenumber for shortwave band", &
         &   units_str="cm-1")
    call out_file%define_variable("band_wavenumber2_sw", &
         &   dim1_name="band_sw", long_name="Upper bound wavenumber for shortwave band", &
         &   units_str="cm-1")
    if (.not. config%do_write_direct_only) then
      call out_file%define_variable("band_flux_up_sw", &
           &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims+1, dim2_name="half_level", dim1_name="band_sw", &
           &   units_str=SW_UNITS_STR, long_name="Upwelling shortwave flux in bands")
      call out_file%define_variable("band_flux_dn_sw", &
           &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims+1, dim2_name="half_level", dim1_name="band_sw", &
           &   units_str=SW_UNITS_STR, long_name="Downwelling shortwave flux in bands")
    end if
    call out_file%define_variable("band_flux_dn_direct_sw", &
         &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
         &   ndims=nfluxdims+1, dim2_name="half_level", dim1_name="band_sw", &
         &   units_str=SW_UNITS_STR, long_name="Direct downwelling shortwave flux in bands")
  end if

  if (do_save_spectrum .and. .not. is_ckd) then
    call out_file%define_variable(trim(spectral_dim_name), &
         &   dim1_name=trim(spectral_dim_name), units_str="cm-1", &
         &   long_name=trim(spectral_dim_name), &
         &   is_double=.true., &
         &   deflate_level=2, shuffle=.true.)
  end if

  if (config%do_write_optical_depth) then
    call out_file%define_variable("optical_depth", &
         &   dim3_name="column", dim2_name="level", dim1_name=trim(spectral_dim_name), &
         &   long_name="Layer optical depth")
  end if

  if (config%do_write_single_scattering_albedo) then
    call out_file%define_variable("single_scattering_albedo", &
         &   dim3_name="column", dim2_name="level", dim1_name=trim(spectral_dim_name), &
         &   long_name="Single scattering albedo")
  end if

  if (config%do_write_asymmetry_factor) then
    call out_file%define_variable("asymmetry_factor", &
         &   dim3_name="column", dim2_name="level", dim1_name=trim(spectral_dim_name), &
         &   long_name="Asymmetry factor")
  end if

  if (config%do_write_spectral_fluxes) then
    if (config%i_spectral_level_index(1) > 0) then
      ! Spectral fluxes are on a different grid from broadband fluxes
      spectral_level_dim_name = "spectral_half_level"
      nspeclev = 1
      do while (config%i_spectral_level_index(nspeclev+1) > 0 .and. nspeclev < 1001)
        nspeclev = nspeclev + 1
      end do

      call out_file%define_dimension(spectral_level_dim_name, nspeclev)
      call out_file%define_variable("spectral_half_level_index", &
           &  data_type_name="int", &
           &  dim1_name=trim(spectral_level_dim_name), &
           &  long_name="Index to half-levels where spectral fluxes are stored")
      call out_file%define_variable("pressure_spectral_hl", &
           &   dim2_name="column", dim1_name=trim(spectral_level_dim_name), &
           &   units_str="Pa", long_name="Pressure where spectral fluxes are stored", &
           &   standard_name="air_pressure")
    else
      nspeclev = 0
      spectral_level_dim_name = "half_level"
    end if
    if (.not. config%do_write_direct_only) then
      call out_file%define_variable("spectral_flux_up_sw", &
           &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims+1, dim2_name=trim(spectral_level_dim_name),  &
           &   dim1_name=trim(spectral_dim_name), &
           &   units_str=SW_UNITS_STR, long_name="Upwelling spectral shortwave flux", &
           &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1,1])
      call out_file%define_variable("spectral_flux_dn_sw", &
           &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims+1, dim2_name=trim(spectral_level_dim_name), &
           &   dim1_name=trim(spectral_dim_name), &
           &   units_str=SW_UNITS_STR, long_name="Downwelling spectral shortwave flux", &
           &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1,1])
    end if
    call out_file%define_variable("spectral_flux_dn_direct_sw", &
         &   dim4_name=trim(flux_dim3_name), dim3_name=trim(flux_dim2_name), &
         &   ndims=nfluxdims+1, dim2_name=trim(spectral_level_dim_name), &
         &   dim1_name=trim(spectral_dim_name), &
         &   units_str=SW_UNITS_STR, long_name="Downwelling direct spectral shortwave flux", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1,1])
  end if

  if (config%do_write_spectral_boundary_fluxes) then
    if (.not. config%do_write_direct_only) then
      call out_file%define_variable("spectral_flux_up_toa_sw", &
           &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims, &
           &   dim1_name=trim(spectral_dim_name), &
           &   units_str=SW_UNITS_STR, long_name="Upwelling spectral shortwave flux at top-of-atmosphere", &
           &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
      call out_file%define_variable("spectral_flux_dn_surf_sw", &
           &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
           &   ndims=nfluxdims, &
           &   dim1_name=trim(spectral_dim_name), &
           &   units_str=SW_UNITS_STR, long_name="Downwelling spectral shortwave flux at surface", &
           &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])
    end if
    call out_file%define_variable("spectral_flux_dn_direct_surf_sw", &
         &   dim3_name=trim(flux_dim3_name), dim2_name=trim(flux_dim2_name), &
         &   ndims=nfluxdims, &
         &   dim1_name=trim(spectral_dim_name), &
         &   units_str=SW_UNITS_STR, long_name="Downwelling direct spectral shortwave flux at surface", &
         &   deflate_level=2, shuffle=.true.,chunksizes=[nspec,1,1])

  end if

  ! Write global variables
  if (is_ckd) then
    call out_file%put_global_attributes(title_str="Shortwave fluxes from a CKD model", &
         &  project_str="CKDMIP", conventions_str="CF-1.7")
  else
    if (do_merge_only) then
      call out_file%put_global_attributes(title_str="Spectral optical depth profiles from multiple gases", &
           &  project_str="CKDMIP", conventions_str="CF-1.7")
      call out_file%put_global_attribute("profile_type", "merge")
    else
      call out_file%put_global_attributes(title_str="Line-by-line shortwave fluxes", &
           &  project_str="CKDMIP", conventions_str="CF-1.7")
    end if
    call out_file%put_global_attribute("constituent_id", molecule_list)
  end if

  if (is_scenario) then
    call out_file%put_global_attribute("scenario", trim(scenario_str))
  end if
  call out_file%put_global_attribute("software_version", CKDMIP_VERSION)

  ! Write coordinate variables
  call out_file%put("pressure_hl", pressure_hl(:,istartcol:iendcol))
  if (allocated(temperature_hl)) then
    call out_file%put("temperature_hl", temperature_hl(:,istartcol:iendcol))
  end if

  if (nspeclev > 0) then
    call out_file%put("pressure_spectral_hl", &
         &  pressure_hl(config%i_spectral_level_index(1:nspeclev),istartcol:iendcol))
  end if

  if (.not. do_merge_only) then
    if (config%use_mu0_dimension) then
      call out_file%put("mu0", config%cos_solar_zenith_angle(1:config%nmu0))
    end if
  end if

  if (nbands > 0) then
    call out_file%put("band_wavenumber1_sw", band_wavenumber1(1:nbands))
    call out_file%put("band_wavenumber2_sw", band_wavenumber2(1:nbands))
  end if

  if (do_save_spectrum .and. .not. is_ckd) then
    call out_file%put(trim(spectral_dim_name), wavenumber_cm1(1:nspec))
  end if

  if (.not. is_ckd) then
    call out_file%put("reference_surface_mole_fraction", target_conc(1:ngas))
  end if

  if (nspeclev > 0) then
    call out_file%put("spectral_half_level_index", &
         &  config%i_spectral_level_index(1:nspeclev))
  end if

  do jcol = istartcol,iendcol

    if (config%iverbose >= 3) then
      write(*,'(a,i0)') 'Column ', jcol
    end if

    allocate(od(nspec,nlev))
    allocate(ssa(nspec,nlev))
    allocate(asymmetry(nspec,nlev))

    ssa       = 0.0_jprb ! Absorbing by default
    asymmetry = 0.0_jprb ! Isotropic scattering by default

    do jmu0 = 1,config%nmu0

      if (is_ckd) then

        if (jmu0 == 1) then
          ! Read from a single CKD file
          call in_file(1)%get(trim(config%optical_depth_name), od, jcol)

          if (in_file(1)%exists(trim(config%single_scattering_albedo_name))) then
            ! File has provided SSA and g
            call in_file(1)%get(trim(config%single_scattering_albedo_name), ssa, jcol)
            if (in_file(1)%exists(trim(config%asymmetry_factor_name))) then
              call in_file(1)%get(trim(config%asymmetry_factor_name), asymmetry, jcol)
            end if
          else if (in_file(1)%exists(trim(config%rayleigh_optical_depth_name))) then
            ! File has provided layer-wise Rayleigh optical depth: add
            ! to total optical depth and calculate SSA (temporarily
            ! using "asymmetry" to store Rayleigh optical depth)
            call in_file(1)%get(trim(config%rayleigh_optical_depth_name), asymmetry, jcol)
            od  = od + asymmetry
            ssa = asymmetry / max(od, LAYER_OD_TOL)
            asymmetry = 0.0_jprb
          else
            error stop 'No Rayleigh scattering variables in CKD optical depth file'
          end if

          call in_file(1)%get(trim(config%incoming_flux_name), ssi, jcol)

        end if

        call calc_shortwave_fluxes(nlev,1,nspec, &
             &  config%cos_solar_zenith_angle(jmu0),ssi,surf_albedo, &
             &  od, ssa, asymmetry, flux_dn_direct, flux_dn, flux_up)

      else

        if (jmu0 == 1) then

          allocate(mole_fraction_fl(nlevorig,ngas))

          ! Read and sum optical depths from multiple gas files
          od = 0.0_jprb

          do jgas = 1,ngas
            if (in_file(jgas)%exists("mole_fraction_fl")) then
              ! We are considering gas absorption: read the mole fraction
              call in_file(jgas)%get("mole_fraction_fl", mole_fraction_fl_1gas, jcol)
              if (size(mole_fraction_fl_1gas) /= nlevorig) then
                error stop 'Error: incorrect size for mole_fraction_fl'
              end if
              if (nsplit > 1) then
                if (in_file(jgas)%exists("mole_fraction_hl")) then
                  call in_file(jgas)%get("mole_fraction_hl", mole_fraction_hl_1gas, jcol)
                else
                  if (.not. allocated(mole_fraction_hl_1gas)) then
                    allocate(mole_fraction_hl_1gas(nlevorig+1))
                  end if
                  mole_fraction_hl_1gas = 1.0_jprb
                end if
              end if
              mole_fraction_scaling = 1.0_jprb
              if (conc_operator(jgas) == I_SCALE) then
                mole_fraction_scaling = conc_factor(jgas)
                if (config%iverbose >= 3) then
                  write(*,'(a,i0,a,f8.3)') '  Reading gas ', jgas, ': scaling by ', conc_factor(jgas)
                end if
              else if (conc_operator(jgas) == I_CONST) then
                mole_fraction_scaling = conc_factor(jgas) / mole_fraction_fl_1gas
                if (config%iverbose >= 3) then
                  write(*,'(a,i0,a)') '  Reading gas ', jgas, ': scaling by variable amount with pressure'
                end if
              else
                if (config%iverbose >= 3) then
                  write(*,'(a,i0)') '  Reading gas ', jgas
                end if
              end if
              mole_fraction_fl_1gas = mole_fraction_fl_1gas * mole_fraction_scaling
              mole_fraction_fl(:,jgas) = mole_fraction_fl_1gas

            else
              ! We are considering Rayleigh scattering, cloud or aerosol:
              ! set mole fraction to missing
              mole_fraction_fl(:,jgas) = -1.0_jprb
              if (conc_operator(jgas) == I_SCALE) then
                mole_fraction_scaling = conc_factor(jgas)
                if (config%iverbose >= 3) then
                  write(*,'(a,i0,a,f8.3)') '  Reading species ', jgas, ': scaling by ', conc_factor(jgas)
                end if
              else
                write(*,'(a,i0)') '  Reading species ', jgas
                mole_fraction_scaling = 1.0_jprb
              end if
            end if
            
            call in_file(jgas)%get(trim(config%optical_depth_name), od_1gas, jcol)
            if (size(od_1gas,2) /= nlevorig) then
              error stop 'Error: incorrect size for optical depth'
            end if
            
            do jlev = 1,nlevorig
              if (mole_fraction_scaling(jlev) /= 1.0_jprb) then
                od_1gas(:,jlev) = od_1gas(:,jlev) * mole_fraction_scaling(jlev)
              end if
            end do
            
            ! Optionally split some layers
            if (nsplit > 1) then
              call split_array_od(isplit1, isplit2, nsplit, od_1gas, mole_fraction_hl_1gas)
            end if

            ! Check if we have single scattering albedo and asymmetry factor
            if (in_file(jgas)%exists(trim(config%single_scattering_albedo_name))) then
              if (.not. in_file(jgas)%exists(trim(config%asymmetry_factor_name))) then
                error stop 'Error: single_scattering_albedo must be accompanied by asymmetry_factor'
              end if
              have_scatterers = .true.
              if (in_file(jgas)%get_rank(trim(config%single_scattering_albedo_name)) == 0 &
                   &  .and. in_file(jgas)%get_rank(trim(config%asymmetry_factor_name)) == 0) then
                ! SSA and g are scalars
                call in_file(jgas)%get(trim(config%single_scattering_albedo_name), ssa_scalar)
                call in_file(jgas)%get(trim(config%asymmetry_factor_name), asymmetry_scalar)
                
                if (do_spectral_subsample) then
                  asymmetry = (asymmetry*ssa*od + (asymmetry_scalar*ssa_scalar) &
                       &                         * od_1gas(istartspecin:iendspecin:config%nspectralstride,:)) &
                       &    / max(ssa*od + ssa_scalar &
                       &                   *od_1gas(istartspecin:iendspecin:config%nspectralstride,:), &
                       &          LAYER_OD_TOL)
                  ssa = (ssa*od + ssa_scalar &
                       &          *od_1gas(istartspecin:iendspecin:config%nspectralstride,:)) &
                       &    / max(od + od_1gas(istartspecin:iendspecin:config%nspectralstride,:), &
                       &          LAYER_OD_TOL)
                else
                  asymmetry = (asymmetry*ssa*od + (asymmetry_scalar*ssa_scalar)*od_1gas) &
                       &    / max(ssa*od + ssa_scalar*od_1gas, LAYER_OD_TOL)
                  ssa = (ssa*od + ssa_scalar*od_1gas) / max(od + od_1gas, LAYER_OD_TOL)
                end if
              else
                ! SSA and g should be matrices
                call in_file(jgas)%get(trim(config%single_scattering_albedo_name), ssa_1gas, jcol)
                call in_file(jgas)%get(trim(config%asymmetry_factor_name), asymmetry_1gas, jcol)

                if (nsplit > 1) then
                  call split_array_dup(isplit1, isplit2, nsplit, ssa_1gas)
                  call split_array_dup(isplit1, isplit2, nsplit, asymmetry_1gas)
                end if

                if (do_spectral_subsample) then
                  asymmetry = (asymmetry*ssa*od + asymmetry_1gas(istartspecin:iendspecin:config%nspectralstride,:) &
                       &                         *ssa_1gas(istartspecin:iendspecin:config%nspectralstride,:) &
                       &                         * od_1gas(istartspecin:iendspecin:config%nspectralstride,:)) &
                       &    / max(ssa*od + ssa_1gas(istartspecin:iendspecin:config%nspectralstride,:) &
                       &                   *od_1gas(istartspecin:iendspecin:config%nspectralstride,:), &
                       &          LAYER_OD_TOL)
                  ssa = (ssa*od + ssa_1gas(istartspecin:iendspecin:config%nspectralstride,:) &
                       &          *od_1gas(istartspecin:iendspecin:config%nspectralstride,:)) &
                       &    / max(od + od_1gas(istartspecin:iendspecin:config%nspectralstride,:), &
                       &          LAYER_OD_TOL)
                else
                  asymmetry = (asymmetry*ssa*od + asymmetry_1gas*ssa_1gas*od_1gas) &
                       &    / max(ssa*od + ssa_1gas*od_1gas, LAYER_OD_TOL)
                  ssa = (ssa*od + ssa_1gas*od_1gas) / max(od + od_1gas, LAYER_OD_TOL)
                end if
                deallocate(ssa_1gas)
                deallocate(asymmetry_1gas)
              end if
              
            else if (have_scatterers) then
              ! Still need to modify ssa due to previous scatterers
              if (do_spectral_subsample) then
                ssa = ssa*od / max(od + od_1gas(istartspecin:iendspecin:config%nspectralstride,:), LAYER_OD_TOL)
              else
                ssa = ssa*od / max(od + od_1gas, LAYER_OD_TOL)
              end if
            end if
            
            ! Add optical depth to the total
            if (do_spectral_subsample) then
              od = od + od_1gas(istartspecin:iendspecin:config%nspectralstride,:)
            else
              od = od + od_1gas
            end if
            
          end do

          deallocate(od_1gas)
          deallocate(mole_fraction_fl_1gas)

        end if ! Is first solar zenith angle

        nblock = (nspec - 1 + config%nblocksize) / config%nblocksize
        
        if (.not. do_merge_only) then
          
          if (config%iverbose >= 3) then
            write(*,'(a,f8.3)') '  Radiative transfer for mu0 = ', &
                 &              config%cos_solar_zenith_angle(jmu0)
            
            write(*,'(a,i0,a,i0)') '    ', nblock, ' blocks of length ', config%nblocksize
          end if
          
          !$OMP PARALLEL DO PRIVATE(istartspec,iendspec) SCHEDULE(RUNTIME)
          do jblock = 1,nblock
            istartspec = (jblock-1) * config%nblocksize + 1
            iendspec = min(istartspec + config%nblocksize - 1, nspec)
            
            ! print *, 'Block ', jblock, ' of ', nblock
            
            call calc_shortwave_fluxes(nlev,istartspec,iendspec, &
                 &  config%cos_solar_zenith_angle(jmu0),ssi,surf_albedo, &
                 &  od, ssa, asymmetry, &
                 &  flux_dn_direct, flux_dn, flux_up)
          end do
          !$OMP END PARALLEL DO
          
        end if

      end if ! CKD vs LBL

      if (config%use_mu0_dimension) then
        icolout = jcol + 1 - istartcol
      else
        icolout = jmu0 + (jcol - istartcol) * config%nmu0
      end if

      if (.not. is_ckd .and. jmu0 == 1) then
        ! Write a slice of mole fraction to this column
        if (nsplit > 1) then
          call split_array(isplit1, isplit2, nsplit, mole_fraction_fl)
        end if
        ! Write a slice of mole fraction to this column
        call out_file%put("mole_fraction_fl", mole_fraction_fl, icolout)
        deallocate(mole_fraction_fl)
      end if

      if (.not. do_merge_only) then
        
        bb_flux_dn_direct = sum(flux_dn_direct, 1)
        bb_flux_dn = sum(flux_dn, 1)
        bb_flux_up = sum(flux_up, 1)
        
        ! Write slice of broadband flux corresponding to this column
        if (config%use_mu0_dimension) then
          if (.not. config%do_write_direct_only) then
            call out_file%put("flux_up_sw", bb_flux_up, jmu0, icolout)
            call out_file%put("flux_dn_sw", bb_flux_dn, jmu0, icolout)
          end if
          call out_file%put("flux_dn_direct_sw", bb_flux_dn_direct, jmu0, icolout)
        else
          call out_file%put("mu0", config%cos_solar_zenith_angle(jmu0), icolout)
          if (.not. config%do_write_direct_only) then
            call out_file%put("flux_up_sw", bb_flux_up, icolout)
            call out_file%put("flux_dn_sw", bb_flux_dn, icolout)
          end if
          call out_file%put("flux_dn_direct_sw", bb_flux_dn_direct, icolout)
        end if
        
        ! Write band fluxes
        if (nbands > 0) then
          do jband = 1,nbands
            do jlev = 1,nlev+1
              band_flux_dn(jband,jlev) = sum(flux_dn(:,jlev), mask=(band_number==jband))
              band_flux_up(jband,jlev) = sum(flux_up(:,jlev), mask=(band_number==jband))
              band_flux_dn_direct(jband,jlev) = sum(flux_dn_direct(:,jlev), mask=(band_number==jband))
            end do
          end do
          if (config%use_mu0_dimension) then
            if (.not. config%do_write_direct_only) then
              call out_file%put("band_flux_up_sw", band_flux_up, &
                   &  jmu0, icolout)
              call out_file%put("band_flux_dn_sw", band_flux_dn, &
                   &  jmu0, icolout)
            end if
            call out_file%put("band_flux_dn_direct_sw", band_flux_dn_direct, &
                 &  jmu0, icolout)
          else
            if (.not. config%do_write_direct_only) then
              call out_file%put("band_flux_up_sw", band_flux_up, icolout)
              call out_file%put("band_flux_dn_sw", band_flux_dn, icolout)
            end if
            call out_file%put("band_flux_dn_direct_sw", band_flux_dn_direct, icolout)
          end if
        end if
        
      end if

      if (jmu0 == 1) then
        ! Write optical depth and other scattering properties of the
        ! mixture
        if (config%do_write_optical_depth) then
          call out_file%put("optical_depth", od, icolout)
        end if
        if (config%do_write_single_scattering_albedo) then
          call out_file%put("single_scattering_albedo", ssa, icolout)
        end if
        if (config%do_write_asymmetry_factor) then
          call out_file%put("asymmetry_factor", asymmetry, icolout)
        end if
      end if
        
      ! Write slice of spectral flux corresponding to this column
      if (config%do_write_spectral_fluxes) then
        if (nspeclev == 0) then
          if (config%use_mu0_dimension) then
            if (.not. config%do_write_direct_only) then
              call out_file%put("spectral_flux_up_sw", flux_up, jmu0, icolout)
              call out_file%put("spectral_flux_dn_sw", flux_dn, jmu0, icolout)
            end if
            call out_file%put("spectral_flux_dn_direct_sw", flux_dn_direct, &
                 &  jmu0, icolout)
          else
            if (.not. config%do_write_direct_only) then
              call out_file%put("spectral_flux_up_sw", flux_up, icolout)
              call out_file%put("spectral_flux_dn_sw", flux_dn, icolout)
            end if
            call out_file%put("spectral_flux_dn_direct_sw", flux_dn_direct, icolout)
          end if
        else
          if (config%use_mu0_dimension) then
            if (.not. config%do_write_direct_only) then
              call out_file%put("spectral_flux_up_sw", &
                   &  flux_up(:,config%i_spectral_level_index(1:nspeclev)), jmu0, icolout)
              call out_file%put("spectral_flux_dn_sw", &
                   &  flux_dn(:,config%i_spectral_level_index(1:nspeclev)), jmu0, icolout)
            end if
            call out_file%put("spectral_flux_dn_direct_sw", &
                 &  flux_dn_direct(:,config%i_spectral_level_index(1:nspeclev)), &
                 &  jmu0, icolout)
          else
            if (.not. config%do_write_direct_only) then
              call out_file%put("spectral_flux_up_sw", &
                   &  flux_up(:,config%i_spectral_level_index(1:nspeclev)), icolout)
              call out_file%put("spectral_flux_dn_sw", &
                   &  flux_dn(:,config%i_spectral_level_index(1:nspeclev)), icolout)
            end if
            call out_file%put("spectral_flux_dn_direct_sw", &
                 &  flux_dn_direct(:,config%i_spectral_level_index(1:nspeclev)), icolout)
          end if
        end if
      end if

      if (config%do_write_spectral_boundary_fluxes) then
        if (config%use_mu0_dimension) then
          if (.not. config%do_write_direct_only) then
            call out_file%put("spectral_flux_up_toa_sw", &
                 &  flux_up(:,1), jmu0, icolout)
            call out_file%put("spectral_flux_dn_surf_sw", &
                 &  flux_dn(:,nlev), jmu0, icolout)
          end if
          call out_file%put("spectral_flux_dn_direct_surf_sw", &
               &  flux_dn_direct(:,nlev), jmu0, icolout)
        else
          if (.not. config%do_write_direct_only) then
            call out_file%put("spectral_flux_up_toa_sw", &
                 &  flux_up(:,1), icolout)
            call out_file%put("spectral_flux_dn_surf_sw", &
                 &  flux_dn(:,nlev), icolout)
          end if
          call out_file%put("spectral_flux_dn_direct_surf_sw", &
               &  flux_dn_direct(:,nlev), icolout)
        end if
      end if
      
    end do

    deallocate(od)
    deallocate(ssa)
    deallocate(asymmetry)

  end do
    
  do jgas = 1,ngas
    call in_file(jgas)%close()
  end do
  
  call out_file%close()

contains

  ! Split layers isplit1-isplit2 each into nsplit separate layers,
  ! duplicating the data therein
  subroutine split_array(isplit1, isplit2, nsplit, arr)

    integer, intent(in) :: isplit1, isplit2, nsplit
    real(jprb), intent(inout), allocatable :: arr(:,:)
    real(jprb), allocatable :: tmp(:,:)
    
    integer :: nspec
    integer :: nlev

    integer :: ilev, jlev, nlevnew

    nspec = size(arr,2)
    nlev  = size(arr,1)

    if (nsplit > 1) then
      write(*,'(a)') '  Splitting layers'
      tmp = arr
      nlevnew = nlev + (nsplit-1)*(isplit2-isplit1+1)
      deallocate(arr)
      allocate(arr(nlevnew,nspec))
      arr(1:isplit1-1,:) = tmp(1:isplit1-1,:)
      ilev = isplit1
      do jlev = isplit1,isplit2
        arr(ilev:ilev+nsplit-1,:) = spread(tmp(jlev,:),1,nsplit)
        ilev = ilev+nsplit
      end do
      arr(ilev:nlevnew,:) = tmp(isplit2+1:nlev,:)
    end if

  end subroutine split_array
  
  ! ! Split layers isplit1-isplit2 each into nsplit separate layers,
  ! ! partitioning the data equally, transposed from the convention
  ! ! above
  ! subroutine split_array_od(isplit1, isplit2, nsplit, arr)

  !   integer, intent(in) :: isplit1, isplit2, nsplit
  !   real(jprb), intent(inout), allocatable :: arr(:,:)
  !   real(jprb), allocatable :: tmp(:,:)
    
  !   integer :: nspec
  !   integer :: nlev

  !   integer :: ilev, jlev, nlevnew

  !   nspec = size(arr,1)
  !   nlev  = size(arr,2)

  !   if (nsplit > 1) then
  !     write(*,'(a)') '  Splitting layers'
  !     tmp = arr
  !     nlevnew = nlev + (nsplit-1)*(isplit2-isplit1+1)
  !     deallocate(arr)
  !     allocate(arr(nspec,nlevnew))
  !     arr(:,1:isplit1-1) = tmp(:,1:isplit1-1)
  !     ilev = isplit1
  !     do jlev = isplit1,isplit2
  !       arr(:,ilev:ilev+nsplit-1) = (1.0_jprb/nsplit)*spread(tmp(:,jlev),2,nsplit)
  !       ilev = ilev+nsplit
  !     end do
  !     arr(:,ilev:nlevnew) = tmp(:,isplit2+1:nlev)
  !   end if

  ! end subroutine split_array_od
  
  ! Split layers isplit1-isplit2 each into nsplit separate layers,
  ! partitioning the data equally, transposed from the convention
  ! above
  subroutine split_array_od(isplit1, isplit2, nsplit, arr, mole_frac_hl)

    integer, intent(in) :: isplit1, isplit2, nsplit
    real(jprb), intent(inout), allocatable :: arr(:,:)
    real(jprb), intent(in) :: mole_frac_hl(:)
    real(jprb), allocatable :: tmp(:,:)
    real(jprb) :: scaling, frac_dev

    integer :: nspec
    integer :: nlev

    integer :: ilev, jlev, nlevnew, jsublev

    nspec = size(arr,1)
    nlev  = size(arr,2)

    if (nsplit > 1) then
      write(*,'(a)') '  Splitting layers'
      tmp = arr
      nlevnew = nlev + (nsplit-1)*(isplit2-isplit1+1)
      deallocate(arr)
      allocate(arr(nspec,nlevnew))
      arr(:,1:isplit1-1) = tmp(:,1:isplit1-1)
      ilev = isplit1
      do jlev = isplit1,isplit2
        ! Fractional deviation of half-level mole fractions from mean
        frac_dev = (mole_frac_hl(jlev+1)-mole_frac_hl(jlev)) &
             &   / real(mole_frac_hl(jlev)+mole_frac_hl(jlev+1),jprb)
        ! Scale to get the same but for the furthest sub-layer
        frac_dev = frac_dev * (nsplit-1)/real(nsplit,jprb)
        !arr(:,ilev:ilev+nsplit-1) = (1.0_jprb/nsplit)*spread(tmp(:,jlev),2,nsplit)
        do jsublev = 1,nsplit
          scaling = frac_dev * (2.0_jprb*(jsublev-1)/real(nsplit-1,jprb) - 1.0_jprb)
          arr(:,ilev) = (1.0_jprb/nsplit)*(scaling+1.0_jprb)*tmp(:,jlev)
          ilev = ilev+1
        end do
      end do
      arr(:,ilev:nlevnew) = tmp(:,isplit2+1:nlev)
    end if

  end subroutine split_array_od
  
  ! Split layers isplit1-isplit2 each into nsplit separate layers,
  ! duplicating the data, transposed from the convention above
  subroutine split_array_dup(isplit1, isplit2, nsplit, arr)

    integer, intent(in) :: isplit1, isplit2, nsplit
    real(jprb), intent(inout), allocatable :: arr(:,:)
    real(jprb), allocatable :: tmp(:,:)
    
    integer :: nspec
    integer :: nlev

    integer :: ilev, jlev, nlevnew

    nspec = size(arr,1)
    nlev  = size(arr,2)

    if (nsplit > 1) then
      write(*,'(a)') '  Splitting layers'
      tmp = arr
      nlevnew = nlev + (nsplit-1)*(isplit2-isplit1+1)
      deallocate(arr)
      allocate(arr(nspec,nlevnew))
      arr(:,1:isplit1-1) = tmp(:,1:isplit1-1)
      ilev = isplit1
      do jlev = isplit1,isplit2
        arr(:,ilev:ilev+nsplit-1) = spread(tmp(:,jlev),2,nsplit)
        ilev = ilev+nsplit
      end do
      arr(:,ilev:nlevnew) = tmp(:,isplit2+1:nlev)
    end if

  end subroutine split_array_dup
  

  ! Split layers isplit1-isplit2 each into nsplit separate layers,
  ! where the array is taken to be on half levels and the data are
  ! interpolated across the additional half levels
  subroutine split_array_hl(isplit1, isplit2, nsplit, arr)

    integer, intent(in) :: isplit1, isplit2, nsplit
    real(jprb), intent(inout), allocatable :: arr(:,:)
    real(jprb), allocatable :: tmp(:,:)
    
    integer :: ncol
    integer :: nlev

    integer :: ilev, jlev, jsublev, nlevnew
    
    nlev = size(arr,1)-1
    ncol = size(arr,2)

    if (nsplit > 1) then
      write(*,*) ' Splitting layers (half-level version), array shape=', shape(arr)
      tmp = arr
      nlevnew = nlev + (nsplit-1)*(isplit2-isplit1+1)
      deallocate(arr)
      allocate(arr(nlevnew+1,ncol))
      arr(1:isplit1,:) = tmp(1:isplit1,:)
      ilev = isplit1
      do jlev = isplit1,isplit2
        do jsublev = 1,nsplit
          ! Linear interpolation
          arr(ilev+jsublev,:) = ((nsplit-real(jsublev,jprb))/nsplit)*tmp(jlev,:) &
               &              +         (real(jsublev,jprb) /nsplit)*tmp(jlev+1,:)
        end do
        ilev = ilev+nsplit
      end do
      arr(ilev:nlevnew+1,:) = tmp(isplit2+1:nlev+1,:)
    end if

  end subroutine split_array_hl
  
  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: ckdmip_sw [arguments]'
    write(*,'(a)') '  where the possible arguments are:'
    write(*,'(a)') '    --cos-sza X         : Set the cosine of the solar zenith angle to X'
    write(*,'(a)') '    --sza X             : Set the solar zenith angle to X degrees'
    write(*,'(a)') '              input.h5  : Add unscaled spectral optical depth from this file'
    write(*,'(a)') '    --scale X input.h5  : Add spectral optical depths scaled by X'
    write(*,'(a)') '    --conc  X input.h5  : Add optical depths scaled so surface mole fraction is X'
    write(*,'(a)') '    --const X input.h5  : Add optical depths with constant mole fraction of X'
    write(*,'(a)') '    --ssi     input.h5  : Read solar spectral irradiance from file'
    write(*,'(a)') '    --ckd     input.h5  : Read optical depth and incoming solar from CKD model'
    write(*,'(a)') '           Note: this option is not compatible with spectral optical depth inputs'
    write(*,'(a)') ' -o|--output  output.h5 : Output file'
    write(*,'(a)') ' -c|--config  conf.nam  : Configure using namelist file'
    write(*,'(a)') ' -m|--merge-only        : Only merge optical depths'
    write(*,'(a)') '    --scenario str      : Add a "scenario" global attribute to the output'
    write(*,'(a)') '    --column-range M N  : Only process columns M to N'
    write(*,'(a)') '    --spectral-range M N: Only process spectral points M to N'
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
  
end program ckdmip_sw
