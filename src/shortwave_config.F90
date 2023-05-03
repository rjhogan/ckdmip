module shortwave_config

  use parkind1, only : jprb

  implicit none

  integer, parameter :: N_MAX_NAME_LENGTH = 61
  integer, parameter :: N_MAX_BANDS       = 1000
  integer, parameter :: N_MAX_MU0         = 1001
  integer, parameter :: N_MAX_HALF_LEVELS = 1001

  type shortwave_config_type
    character(len=N_MAX_NAME_LENGTH) :: optical_depth_name = 'optical_depth'
    character(len=N_MAX_NAME_LENGTH) :: rayleigh_optical_depth_name &
         &  = 'rayleigh_optical_depth'
    character(len=N_MAX_NAME_LENGTH) :: single_scattering_albedo_name &
         &  = 'single_scattering_albedo'
    character(len=N_MAX_NAME_LENGTH) :: asymmetry_factor_name &
         &  = 'asymmetry_factor'
    character(len=N_MAX_NAME_LENGTH) :: incoming_flux_name &
         &  = 'incoming_flux'
    character(len=N_MAX_NAME_LENGTH) :: pressure_name = 'pressure_hl'
    character(len=N_MAX_NAME_LENGTH) :: temperature_name = 'temperature_hl'
    character(len=N_MAX_NAME_LENGTH) :: wavenumber_name = 'wavenumber'
    character(len=N_MAX_NAME_LENGTH) :: surf_albedo_name = 'surf_albedo'

    real(jprb), dimension(N_MAX_BANDS) :: band_wavenumber1, band_wavenumber2
    real(jprb), dimension(N_MAX_MU0)   :: cos_solar_zenith_angle

    real(jprb) :: pressure_scaling = 1.0_jprb

    real(jprb) :: surf_albedo = 0.06_jprb

    ! Indices of the half-levels where spectral fluxes are to be
    ! written
    integer, dimension(N_MAX_HALF_LEVELS) :: i_spectral_level_index

    integer :: nblocksize = 64
    integer :: iverbose   = 3
    integer :: nspectralstride = 1
    integer :: nmu0 = 0

    logical :: do_write_optical_depth = .false.
    logical :: do_write_single_scattering_albedo = .false.
    logical :: do_write_asymmetry_factor = .false.
    logical :: do_write_spectral_fluxes = .false.
    logical :: do_write_spectral_boundary_fluxes = .false.
    logical :: do_write_direct_only = .false.

    ! Do we have a dimension in the output file for the cosine of the
    ! solar zenith angle?
    logical :: use_mu0_dimension = .false.

  contains
    procedure :: read => read_config_from_namelist
    procedure :: get_bands

  end type shortwave_config_type

contains

  subroutine read_config_from_namelist(this, file_name)
    class(shortwave_config_type), intent(inout) :: this
    character(*), intent(in)                   :: file_name

    integer :: iosopen ! Status after calling open

    character(len=N_MAX_NAME_LENGTH) :: optical_depth_name
    character(len=N_MAX_NAME_LENGTH) :: rayleigh_optical_depth_name
    character(len=N_MAX_NAME_LENGTH) :: single_scattering_albedo_name
    character(len=N_MAX_NAME_LENGTH) :: asymmetry_factor_name
    character(len=N_MAX_NAME_LENGTH) :: incoming_flux_name
    character(len=N_MAX_NAME_LENGTH) :: pressure_name
    character(len=N_MAX_NAME_LENGTH) :: temperature_name
    character(len=N_MAX_NAME_LENGTH) :: wavenumber_name
    character(len=N_MAX_NAME_LENGTH) :: surf_albedo_name

    real(jprb), dimension(N_MAX_BANDS) :: band_wavenumber1, band_wavenumber2
    real(jprb), dimension(N_MAX_MU0)   :: cos_solar_zenith_angle

    real(jprb) :: pressure_scaling

    real(jprb) :: surf_albedo

    integer, dimension(N_MAX_HALF_LEVELS) :: i_spectral_level_index

    integer :: nblocksize, iverbose, nspectralstride

    logical :: do_write_optical_depth
    logical :: do_write_single_scattering_albedo
    logical :: do_write_asymmetry_factor
    logical :: do_write_spectral_fluxes
    logical :: do_write_spectral_boundary_fluxes
    logical :: do_write_direct_only
    logical :: use_mu0_dimension

    namelist /shortwave_config/ optical_depth_name, pressure_name, &
         &  temperature_name, wavenumber_name, surf_albedo, &
         &  surf_albedo_name, pressure_scaling, nblocksize, iverbose, &
         &  do_write_optical_depth, do_write_single_scattering_albedo, &
         &  do_write_asymmetry_factor, do_write_spectral_fluxes, &
         &  do_write_direct_only, do_write_spectral_boundary_fluxes, &
         &  nspectralstride, band_wavenumber1, band_wavenumber2, &
         &  cos_solar_zenith_angle, use_mu0_dimension, &
         &  rayleigh_optical_depth_name, single_scattering_albedo_name, &
         &  asymmetry_factor_name, incoming_flux_name, i_spectral_level_index

    optical_depth_name = this%optical_depth_name
    rayleigh_optical_depth_name = this%rayleigh_optical_depth_name
    single_scattering_albedo_name = this%single_scattering_albedo_name
    asymmetry_factor_name = this%asymmetry_factor_name
    incoming_flux_name = this%incoming_flux_name
    pressure_name      = this%pressure_name
    temperature_name   = this%temperature_name
    wavenumber_name    = this%wavenumber_name
    surf_albedo_name   = this%surf_albedo_name
    pressure_scaling   = this%pressure_scaling
    surf_albedo        = this%surf_albedo
    nblocksize         = this%nblocksize
    iverbose           = this%iverbose
    nspectralstride    = this%nspectralstride
    do_write_optical_depth = this%do_write_optical_depth
    do_write_single_scattering_albedo = this%do_write_single_scattering_albedo
    do_write_asymmetry_factor = this%do_write_asymmetry_factor
    do_write_spectral_fluxes = this%do_write_spectral_fluxes
    do_write_spectral_boundary_fluxes = this%do_write_spectral_boundary_fluxes
    do_write_direct_only = this%do_write_direct_only
    band_wavenumber1 = -1.0_jprb
    band_wavenumber2 = -1.0_jprb
    cos_solar_zenith_angle = -2.0_jprb
    use_mu0_dimension = .false.
    i_spectral_level_index = 0

    open(unit=10, iostat=iosopen, file=trim(file_name))
    if (iosopen /= 0) then
      stop 'Error: could not open namelist file'
    else
      read(unit=10, nml=shortwave_config)
      close(unit=10)
    end if

    this%optical_depth_name = optical_depth_name
    this%pressure_name      = pressure_name
    this%temperature_name   = temperature_name
    this%wavenumber_name    = wavenumber_name
    this%surf_albedo_name   = surf_albedo_name
    this%pressure_scaling   = pressure_scaling
    this%surf_albedo        = surf_albedo
    this%nblocksize         = nblocksize
    this%iverbose           = iverbose
    this%nspectralstride    = nspectralstride
    this%do_write_optical_depth = do_write_optical_depth
    this%do_write_single_scattering_albedo = do_write_single_scattering_albedo
    this%do_write_asymmetry_factor = do_write_asymmetry_factor
    this%do_write_spectral_fluxes = do_write_spectral_fluxes
    this%do_write_spectral_boundary_fluxes = do_write_spectral_boundary_fluxes
    this%do_write_direct_only = do_write_direct_only
    this%band_wavenumber1 = band_wavenumber1
    this%band_wavenumber2 = band_wavenumber2
    this%use_mu0_dimension = use_mu0_dimension
    this%rayleigh_optical_depth_name = rayleigh_optical_depth_name
    this%single_scattering_albedo_name = single_scattering_albedo_name
    this%asymmetry_factor_name = asymmetry_factor_name
    this%incoming_flux_name = incoming_flux_name
    this%i_spectral_level_index = i_spectral_level_index

    if (this%nmu0 == 0) then
      ! Cosine of solar zenith angle not already set on command line
      this%cos_solar_zenith_angle = cos_solar_zenith_angle

      ! Number of solar zenith angles to use depends on how many have
      ! been set by user in namelist.
      !this%nmu0 = findloc(cos_solar_zenith_angle == 2.0_jprb)-1
      ! Pending the availability of findloc in most Fortran
      ! compilers...
      do while (cos_solar_zenith_angle(this%nmu0+1) > -2.0_jprb &
           &    .and. this%nmu0 < N_MAX_MU0)
        this%nmu0 = this%nmu0 + 1
      end do
    end if

  end subroutine read_config_from_namelist

  ! Get the band bounds and return the number of bands
  integer function get_bands(this, band_wavenumber1, band_wavenumber2)

    class(shortwave_config_type),           intent(in)    :: this
    real(jprb), allocatable, dimension(:), intent(inout) :: band_wavenumber1, band_wavenumber2
    
    get_bands = 0
    do while (this%band_wavenumber1(get_bands+1) >= 0.0_jprb &
         &    .and. this%band_wavenumber2(get_bands+1) >= 0.0_jprb &
         &    .and. get_bands < N_MAX_BANDS-1)
      get_bands = get_bands + 1
    end do

    if (get_bands > 0) then
      if (.not. allocated(band_wavenumber1)) then
        allocate(band_wavenumber1(get_bands))
      end if
      if (.not. allocated(band_wavenumber2)) then
        allocate(band_wavenumber2(get_bands))
      end if

      band_wavenumber1(1:get_bands) = this%band_wavenumber1(1:get_bands)
      band_wavenumber2(1:get_bands) = this%band_wavenumber2(1:get_bands)
    end if

  end function get_bands

end module shortwave_config
