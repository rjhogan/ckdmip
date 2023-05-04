! ckdmip_tool.F90 - general utility for CKDMIP
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

program ckdmip_tool

  use parkind1,        only : jprb, jpib
  use easy_netcdf,     only : netcdf_file
  use interpolation,   only : interpolate
  
  implicit none

#include "ckdmip_version.h"

  integer, parameter :: TASK_NONE     = 0
  integer, parameter :: TASK_RAYLEIGH = 1
  integer, parameter :: TASK_CLOUD    = 2
  integer, parameter :: TASK_SSI      = 3
  integer, parameter :: TASK_MASSEXT  = 4

  integer, parameter :: iverbose = 3

  real(jprb), parameter :: ACCEL_DUE_TO_GRAVITY  = 9.80665_jprb  ! m s-2
  real(jprb), parameter :: DRY_AIR_MOLAR_MASS    = 0.028970_jprb ! kg mol-1
  real(jprb), parameter :: AVOGADRO_CONSTANT     = 6.02214076e23_jprb ! mol-1
  real(jprb), parameter :: H2O_MOLAR_MASS        = 0.0180152833_jprb ! kg mol-1

  integer :: task

  ! Number of levels, spectral intervals, gases and profiles
  integer :: nlev, nspecin, ncol, ncolout, nspec, ninspec, narg, iarg, ire, nwnorig, nre

  ! Column, level loop index
  integer :: jcol, jlev

  integer :: istartcol, iendcol
  integer :: istartlev, iendlev

  real(jprb) :: total_water_path
  real(jprb) :: effective_radius

  ! Wavenumber in cm-1 at original and final resolution
  real(jprb), allocatable, dimension(:)   :: wavenumber_orig, wavenumber

  ! Wavelength in microns, for Rayleigh scattering
  real(jprb), allocatable, dimension(:)   :: wavelength_um

  ! Rayleigh scattering properties; arayleigh is the extinction per
  ! unit mass of air (m2 kg-1)
  real(jprb), allocatable, dimension(:)   :: xrayleigh, arayleigh

  ! Mass of layer, kg m-2
  real(jprb) :: layer_mass
  real(jprb), allocatable, dimension(:) :: layer_mass_full
  real(jprb), allocatable, dimension(:) :: mole_fraction_fl

  ! Wavenumber spacing
  real(jprb), allocatable, dimension(:)   :: dwav

  ! Pressure at half-levels (Pa)
  real(jprb), allocatable, dimension(:,:)   :: pressure_hl

  ! Optical depth (nspec,nlev)
  real(jprb), allocatable, dimension(:,:)   :: optical_depth

  ! Single scattering albedo and asymmetry factor on full wavenumber
  ! and height grid (nspec,nlev)
  real(jprb), allocatable, dimension(:,:)   :: ssa_full, asymmetry_full

  ! Solar spectral irradiance at original resolution, and final
  ! resolution (W m-2 cm)
  real(jprb), allocatable, dimension(:) :: ssi_orig, ssi

  ! Cloud mass-extinction coefficient, single-scattering albedo and
  ! asymmetry factor at the original spectral resolution
  real(jprb), allocatable, dimension(:,:) :: mec_orig, ssa_orig, asymmetry_orig

  ! Effective radius look-up table
  real(jprb), allocatable, dimension(:)   :: effective_radius_orig

  ! Cloud properties interpolated to the requested effective radius
  real(jprb), allocatable, dimension(:) :: mec_re, ssa_re, asymmetry_re

  ! Cloud properties interpolated to the full wavenumber grid
  real(jprb), allocatable, dimension(:) :: mec, ssa, asymmetry

  ! Chunk sizes in compression of optical depth
  integer :: chunksizes(3)

  ! Total solar irradiance (W m-2)
  real(jprb) :: tsi

  ! Scaling to obtain requested water path; interpolation weights
  real(jprb) :: scaling, weight1, weight2

  ! Input and output files
  type(netcdf_file) :: in_file, out_file

  character(len=512) :: argument_str, file_name, out_file_name
  character(len=512) :: molecule_str, scenario_str, title_str, references_str

  logical :: is_scenario, is_out_file, do_delta_eddington, do_cloud_scattering

  integer :: istatus

  task = TASK_NONE
  is_scenario = .false.
  istartcol   = 1
  iendcol     = 0
  is_out_file = .false.
  ncol        = 0
  ncolout     = 0
  nlev        = 0
  iarg        = 1
  do_cloud_scattering = .true.
  do_delta_eddington = .false.

  ! Loop over arguments
  narg = command_argument_count()

  do while (iarg <= narg)
    
    call get_command_argument(iarg, argument_str, status=istatus)
    if (istatus /= 0) then
      write(*,'(a)') 'Failed to read argument as string of length < 512'
      call print_usage_and_exit()
    end if
    
    if (argument_str == '-o' .or. argument_str == '--output') then

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

      is_out_file = .true.

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
        stop
      end if

      if (iverbose >= 3) then
        write(*,'(a,i0,a,i0)') 'Processing columns ', istartcol, ' to ', iendcol
      end if

    else if (argument_str == '--rayleigh') then

      if (task /= TASK_NONE) then
        stop 'More than one task requested'
      end if
      task = TASK_RAYLEIGH

    else if (argument_str == '--grid') then

      if (iarg == narg) then
        write(*,'(a)') '"--grid" must be followed by the name of a file'
        call print_usage_and_exit()
      end if

      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of grid file as string of length < 512'
        call print_usage_and_exit()
      end if

      call in_file%open(trim(file_name), iverbose=iverbose)
      call in_file%get("wavenumber", wavenumber)
      call in_file%get("pressure_hl", pressure_hl)

      call in_file%close()

    else if (argument_str == '--ssi') then
      
      if (iarg >= narg-1) then
        write(*,'(a)') '"--ssi" must be followed by two arguments'
        call print_usage_and_exit()
      else if (task /= TASK_NONE) then
        stop 'More than one task requested'
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) tsi
      
      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of SSI file as string of length < 512'
        call print_usage_and_exit()
      end if

      call in_file%open(trim(file_name), iverbose=iverbose)
      call in_file%get("wavenumber", wavenumber_orig)
      call in_file%get("mean_solar_spectral_irradiance", ssi_orig)
      call in_file%get_global_attribute("title", title_str)
      call in_file%close()

      task = TASK_SSI

    else if (argument_str == '--cloud' .or. argument_str == '--delta-cloud' &
         & .or. argument_str == '--absorption-cloud') then

      do_cloud_scattering = .true.
      do_delta_eddington = .false.
      if (argument_str == '--delta-cloud') then
        do_delta_eddington = .true.
      else if (argument_str == '--absorption-cloud') then
        do_cloud_scattering = .false.
      end if

      if (iarg >= narg-4) then
        write(*,'(a)') '"--cloud" must be followed by four arguments'
        call print_usage_and_exit()
      else if (task /= TASK_NONE) then
        stop 'More than one task requested'
      end if

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) total_water_path
      
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) effective_radius
      ! Convert from um to m
      effective_radius = effective_radius * 1.0e-6
      
      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) istartlev

      iarg = iarg+1
      call get_command_argument(iarg, argument_str)
      read(argument_str, *) iendlev
      
      iarg = iarg+1
      call get_command_argument(iarg, file_name, status=istatus)
      if (istatus /= 0) then
        write(*,'(a)') 'Failed to read name of cloud optics file as string of length < 512'
        call print_usage_and_exit()
      end if

      ! Load cloud optics file
      call in_file%open(trim(file_name), iverbose=iverbose)
      call in_file%get("wavenumber", wavenumber_orig)
      call in_file%get("effective_radius", effective_radius_orig)
      call in_file%get("mass_extinction_coefficient", mec_orig)
      call in_file%get("single_scattering_albedo", ssa_orig)
      call in_file%get("asymmetry_factor", asymmetry_orig)
      call in_file%close()

      ! Interpolate to requested effective radius
      nwnorig = size(wavenumber_orig)
      nre = size(effective_radius_orig)
      if (effective_radius < effective_radius_orig(1) &
           &  .or. effective_radius > effective_radius_orig(nre)) then
        write(*,'(a,e11.5,a,e11.5,a,e11.5,a)') 'Error: effective radius of ', effective_radius, &
             &  ' m outside valid range of ', effective_radius_orig(1), '-', &
             &  effective_radius_orig(nre), ' m'
        stop
      end if

      ire = 1
      do while (ire < nre-1 .and. effective_radius > effective_radius_orig(ire+1))
        ire = ire+1
      end do

      weight1 = (effective_radius_orig(ire+1)-effective_radius) &
           &  / (effective_radius_orig(ire+1)-effective_radius_orig(ire))
      weight2 = (effective_radius-effective_radius_orig(ire)) &
           &  / (effective_radius_orig(ire+1)-effective_radius_orig(ire))

      !print *, weight1, weight2

      if (do_cloud_scattering) then

        allocate(mec_re(nwnorig))
        allocate(ssa_re(nwnorig))
        allocate(asymmetry_re(nwnorig))
        mec_re = weight1 * mec_orig(:,ire) + weight2 * mec_orig(:,ire+1)
        ssa_re = weight1 * ssa_orig(:,ire) + weight2 * ssa_orig(:,ire+1)
        asymmetry_re = weight1 * asymmetry_orig(:,ire) + weight2 * asymmetry_orig(:,ire+1)

      else

        write(*,'(a)') 'Computing cloud absorption optical depth'

        ! When scattering is ignored we store the absorption optical depth
        allocate(mec_re(nwnorig))
        mec_re = weight1 * mec_orig(:,ire)*(1.0_jprb-ssa_orig(:,ire)) &
             & + weight2 * mec_orig(:,ire+1)*(1.0_jprb-ssa_orig(:,ire+1))

      end if

      deallocate(mec_orig)
      deallocate(ssa_orig)
      deallocate(asymmetry_orig)

      if (do_delta_eddington) then
        write(*,'(a)') 'Performing in-place delta-Eddington scaling of cloud optical properties'
        call delta_eddington(mec_re, ssa_re, asymmetry_re)
      else
        write(*,'(a)') 'Not performing delta-Eddington scaling of cloud optical properties'
      end if

      task = TASK_CLOUD

    else
      
      write(*,'(a,a,a)') 'Argument "', trim(argument_str), '" not understood'
      call print_usage_and_exit()

    end if

    iarg = iarg + 1

  end do

  if (task == TASK_NONE) then
    write(*,'(a)') 'No task specified'
    call print_usage_and_exit()
  else if (.not. is_out_file) then
    write(*,'(a)') 'Output file not specified'
    call print_usage_and_exit()
  else if (.not. allocated(wavenumber) .and. task /= TASK_MASSEXT) then
    write(*,'(a)') 'Task requested requires a --grid argument'
    call print_usage_and_exit()
  end if

  nspec = size(wavenumber)

  if (allocated(pressure_hl) .and. task /= TASK_SSI) then
    ncol = size(pressure_hl,2)
    nlev = size(pressure_hl,1)-1
    
    if (iendcol <= 0) then
      iendcol = ncol
    else
      iendcol = iendcol
    end if
    ncolout = iendcol + 1 - istartcol
  end if
  ! Open output file
  call out_file%create(trim(out_file_name), is_hdf5_file=.true., iverbose=iverbose)

  if (ncolout > 0) then
    ! File will be large: don't use an unlimited dimension, enabling
    ! compression
    call out_file%define_dimension("column", ncolout)
  end if
  if (nlev > 0) then
    call out_file%define_dimension("level", nlev)
    call out_file%define_dimension("half_level", nlev+1)
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="Pa", long_name="Pressure", &
         &   standard_name="air_pressure")
   end if

   call out_file%define_dimension("wavenumber", nspec)
   call out_file%define_variable("wavenumber", &
        &   dim1_name="wavenumber", units_str="cm-1", long_name="Wavenumber", &
        &   deflate_level=2, shuffle=.true., is_double=.true.)

   if (task == TASK_RAYLEIGH) then
     allocate(wavelength_um(nspec))
     wavelength_um = 1.0e4_jprb / wavenumber
     allocate(arayleigh(nspec))
!#define USE_NICOLET_RAYLEIGH 1
#ifdef USE_NICOLET_RAYLEIGH
     ! Marcel Nicolet, On the molecular scattering in the terrestrial
     ! atmosphere : An empirical formula for its calculation in the
     ! homosphere, Planetary and Space Science, Volume 32, Issue 11,
     ! 1984, Pages 1467-1468.
     allocate(xrayleigh(nspec))
     where (wavelength_um < 0.550_jprb)
       xrayleigh = 0.389_jprb * wavelength_um + 0.09426_jprb / wavelength_um - 0.3228_jprb
     elsewhere
       xrayleigh = 0.04_jprb
     end where
     arayleigh = 4.02e-32_jprb * wavelength_um ** -(4.0_jprb + xrayleigh)
     deallocate(xrayleigh)
     references_str = "Nicolet, M., 1984: On the molecular scattering in the terrestrial atmosphere: " &
          &  //"An empirical formula for its calculation in the homosphere. Planetary and Space Sci., " &
          &  //"32, 1467-1468."
#else
     ! Bucholtz (1995) model, probably superior than Nicolet, giving
     ! molecular cross-section in m2
     where (wavelength_um < 0.5_jprb)
       arayleigh = 3.01577e-32_jprb * wavelength_um ** -( 3.55212_jprb &
            &                                            +1.35579_jprb*wavelength_um &
            &                                            +0.11563_jprb/wavelength_um)
     elsewhere
       arayleigh = 4.01061e-32_jprb * wavelength_um ** -( 3.99668_jprb &
            &                                            +0.00110298_jprb*wavelength_um &
            &                                            +0.0271393_jprb/wavelength_um)
     end where
     references_str = "Bucholtz, A., 1995: Rayleigh-scattering calculations for the terrestrial " &
          &  //"atmosphere. Applied Optics, 34, 2765-2773."
#endif

     ! Convert to cross-section per kg of air (m2 kg-1)
     arayleigh = arayleigh * (AVOGADRO_CONSTANT / DRY_AIR_MOLAR_MASS)

     deallocate(wavelength_um)

     call out_file%define_variable("single_scattering_albedo", &
          &  long_name="Single scattering albedo due to Rayleigh scattering")
     call out_file%define_variable("asymmetry_factor", &
          &  long_name="Asymmetry factor due to Rayleigh scattering")
     
     ! Choose chunk sizes for compression, allowing for the limit on
     ! the chunk size being 2**32 (which is apparently less)
     chunksizes = [nspec,nlev,1]
     if (int(nspec,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
       chunksizes(2) = 1
     end if

     call out_file%define_variable("column_optical_depth", &
          &  dim2_name="column", dim1_name="wavenumber", &
          &  long_name="Column optical depth due to Rayleigh scattering", &
          &  deflate_level=2, shuffle=.true., chunksizes=[nspec,1])
     call out_file%define_variable("optical_depth", &
          &  dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
          &  long_name="Layer optical depth due to Rayleigh scattering", &
          &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)

     call out_file%put_global_attributes(title_str="Spectral Rayleigh scattering profiles", &
          &  project_str="CKDMIP", conventions_str="CF-1.7", references_str=references_str)
     call out_file%put_global_attribute("constituent_id", "rayleigh")
     call out_file%put_global_attribute("profile_type", "present")
     call out_file%put_global_attribute("software_version", CKDMIP_VERSION)

     call out_file%put("pressure_hl", pressure_hl)
     call out_file%put("wavenumber", wavenumber)
     call out_file%put("single_scattering_albedo", 1.0_jprb)
     call out_file%put("asymmetry_factor", 0.0_jprb)

     allocate(optical_depth(nspec,nlev))
     do jcol = istartcol,iendcol
       do jlev = 1,nlev
         layer_mass = (pressure_hl(jlev+1,jcol) - pressure_hl(jlev,jcol)) / ACCEL_DUE_TO_GRAVITY
         optical_depth(:,jlev) = layer_mass * arayleigh
       end do
       call out_file%put("column_optical_depth", sum(optical_depth,2), jcol-istartcol+1)
       call out_file%put("optical_depth", optical_depth, jcol-istartcol+1)
     end do

   else if (task == TASK_CLOUD) then

     allocate(mec(nspec))
     call interpolate(wavenumber_orig, mec_re, wavenumber, mec)
     deallocate(mec_re)

     if (do_cloud_scattering) then
       allocate(ssa(nspec))
       allocate(asymmetry(nspec))
       call interpolate(wavenumber_orig, ssa_re, wavenumber, ssa)
       call interpolate(wavenumber_orig, asymmetry_re, wavenumber, asymmetry)
       deallocate(ssa_re)
       deallocate(asymmetry_re)
     end if
     
     ! Choose chunk sizes for compression, allowing for the limit on
     ! the chunk size being 2**32 (which is apparently less)
     chunksizes = [nspec,nlev,1]
     if (int(nspec,jpib)*int(nlev,jpib) >= int(2,jpib)**30) then
       chunksizes(2) = 1
     end if

     call out_file%define_variable("mole_fraction_fl", &
          &  dim2_name="column", dim1_name="level", &
          &  units_str="1", long_name="Mole fraction", fill_value=-1.0_jprb)

     if (do_cloud_scattering) then
       call out_file%define_variable("optical_depth", &
            &  dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
            &  long_name="Layer optical depth", &
            &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)
       call out_file%define_variable("single_scattering_albedo", &
            &  dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
            &  long_name="Single scattering albedo", &
            &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)
       call out_file%define_variable("asymmetry_factor", &
            &  dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
            &  long_name="Asymmetry factor", &
            &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)
       call out_file%put_global_attributes(title_str="Spectral cloud scattering profiles", &
            &  project_str="CKDMIP", conventions_str="CF-1.7")
     else
       call out_file%define_variable("optical_depth", &
            &  dim3_name="column", dim2_name="level", dim1_name="wavenumber", &
            &  long_name="Layer absorption optical depth", &
            &  deflate_level=2, shuffle=.true., chunksizes=chunksizes)
       call out_file%put_global_attributes(title_str="Spectral cloud absorption profiles", &
            &  project_str="CKDMIP", conventions_str="CF-1.7")
     end if

     call out_file%put_global_attribute("constituent_id", "cloud")
     call out_file%put_global_attribute("profile_type", "present")
     call out_file%put_global_attribute("software_version", CKDMIP_VERSION)

     call out_file%put("pressure_hl", pressure_hl(:,istartcol:iendcol))
     call out_file%put("wavenumber", wavenumber)

     allocate(optical_depth(nspec,nlev))
     optical_depth = 0.0_jprb

     if (do_cloud_scattering) then
       allocate(ssa_full(nspec,nlev))
       allocate(asymmetry_full(nspec,nlev))
       ssa_full = 0.0_jprb
       asymmetry_full = 0.0_jprb
       do jlev = istartlev,iendlev
         ssa_full(:,jlev) = ssa
         asymmetry_full(:,jlev) = asymmetry
       end do
     end if

     allocate(layer_mass_full(nlev))
     allocate(mole_fraction_fl(nlev))

     mole_fraction_fl = 0.0_jprb

     do jcol = istartcol,iendcol
       layer_mass_full = (pressure_hl(2:nlev+1,jcol) - pressure_hl(1:nlev,jcol)) / ACCEL_DUE_TO_GRAVITY
       scaling = total_water_path / sum(layer_mass_full(istartlev:iendlev))
       do jlev = istartlev,iendlev
         optical_depth(:,jlev) = scaling * layer_mass_full(jlev) * mec
         mole_fraction_fl(jlev) = scaling / H2O_MOLAR_MASS
       end do
       call out_file%put("mole_fraction_fl", mole_fraction_fl, jcol-istartcol+1)
       call out_file%put("optical_depth", optical_depth, jcol-istartcol+1)
       if (do_cloud_scattering) then
         call out_file%put("single_scattering_albedo", ssa_full, jcol-istartcol+1)
         call out_file%put("asymmetry_factor", asymmetry_full, jcol-istartcol+1)
       end if
     end do

   else if (task == TASK_SSI) then

     allocate(ssi(nspec))
     call interpolate(wavenumber_orig, ssi_orig, wavenumber, ssi)
     allocate(dwav(nspec))
     dwav(1:nspec-1) = abs(wavenumber(2:nspec) - wavenumber(1:nspec-1))
     dwav(nspec) = dwav(nspec-1)
     ssi = ssi * dwav
     ssi = ssi * (tsi / sum(ssi))

     call out_file%define_variable("total_solar_irradiance", &
          &  long_name="Total solar irradiance", units_str="W m-2")
     call out_file%define_variable("solar_spectral_irradiance", &
          &  dim1_name="wavenumber", long_name="Solar spectral irradiance", &
          &  units_str="W m-2", deflate_level=2, shuffle=.true., &
          &  comment_str= &
          &  "This variable is the solar irradiance in each wavenumber interval. To obtain solar irradiance"//new_line('a') &
          &  //"per unit wavenumber at a point in the spectrum, divide by the difference in wavenumber between"//new_line('a') &
          &  //"that point and the next one.")
     call out_file%put_global_attributes(title_str=title_str, project_str="CKDMIP")
     call out_file%put("wavenumber", wavenumber)
     call out_file%put("total_solar_irradiance", tsi)
     call out_file%put("solar_spectral_irradiance", ssi)

   else if (task == TASK_MASSEXT) then

   end if

  call out_file%close()

contains

  subroutine print_usage_and_exit()

    write(*,'(a)') 'Usage: ckdmip_tool [arguments]'
    write(*,'(a)') '  where the possible arguments are:'
    write(*,'(a)') '    --grid input.h5    : Read pressure and wavenumber from file'
    write(*,'(a)') '    --rayleigh         : Create a Rayleigh scattering spectral file'
    write(*,'(a)') '    --cloud WP RE M N input.nc'
    write(*,'(a)') '                       : Create a cloud spectral file: water path WP (kg m-2),'
    write(*,'(a)') '                         effective radius RE (um), between layers M and N'
    write(*,'(a)') '                         inclusive, reading cloud properties from input file'
    write(*,'(a)') '                         input.nc'
    write(*,'(a)') '    --delta-cloud WP RE M N input.nc'
    write(*,'(a)') '                       : As --cloud but also apply delta-Eddington scaling'
    write(*,'(a)') '    --absorption-cloud WP RE M N input.nc'
    write(*,'(a)') '                       : as --cloud but for a non-scattering cloud'
    write(*,'(a)') '    --ssi TSI input.h5 : Create a solar spectral irradiance file with a total'
    write(*,'(a)') '                         solar irradiance of TSI (W m-2)'
    write(*,'(a)') '    --od-to-mass-ext input.h5'
    write(*,'(a)') '                       : Convert gas optical depths to mass extinction'
    write(*,'(a)') '                         coefficients'
    write(*,'(a)') ' -o|--output output.h5 : Output file'
    write(*,'(a)') '    --scenario str     : Add a "scenario" global attribute to the output'
    write(*,'(a)') '    --column-range M N : Only process columns M to N'
    stop

  end subroutine print_usage_and_exit

  
  !---------------------------------------------------------------------
  ! Perform in-place delta-Eddington scaling of the phase function
  elemental subroutine delta_eddington(od, ssa, g)

    use parkind1, only : jprb
    
    ! Total optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), intent(inout) :: od, ssa, g
    
    ! Fraction of the phase function deemed to be in the forward lobe
    ! and therefore treated as if it is not scattered at all
    real(jprb) :: f
    
    f   = g*g
    od  = od * (1.0_jprb - ssa*f)
    ssa = ssa * (1.0_jprb - f) / (1.0_jprb - ssa*f)
    g   = g / (1.0_jprb + g)
  
  end subroutine delta_eddington

end program ckdmip_tool
