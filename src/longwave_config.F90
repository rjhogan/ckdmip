module longwave_config

  use parkind1, only : jprb

  implicit none

  integer, parameter :: N_MAX_NAME_LENGTH = 61
  integer, parameter :: N_MAX_BANDS       = 1000

  type longwave_config_type
    character(len=N_MAX_NAME_LENGTH) :: optical_depth_name = 'optical_depth'
    character(len=N_MAX_NAME_LENGTH) :: pressure_name = 'pressure_hl'
    character(len=N_MAX_NAME_LENGTH) :: temperature_name = 'temperature_hl'
    character(len=N_MAX_NAME_LENGTH) :: wavenumber_name = 'wavenumber'
    character(len=N_MAX_NAME_LENGTH) :: planck_name = 'planck_hl'
    character(len=N_MAX_NAME_LENGTH) :: surf_emission_name = 'surf_emission'
    character(len=N_MAX_NAME_LENGTH) :: surf_emissivity_name = 'surf_emissivity'
    character(len=N_MAX_NAME_LENGTH) :: surf_temperature_name = 'skin_temperature'

    real(jprb), dimension(N_MAX_BANDS) :: band_wavenumber1, band_wavenumber2

    real(jprb) :: pressure_scaling = 1.0_jprb
    real(jprb) :: surf_temperature = -1.0_jprb

    integer :: nblocksize = 64
    integer :: iverbose   = 3
    integer :: nspectralstride = 1

    ! Number of angles per hemisphere; 0 = classic method with
    ! diffusivity of 1.66, otherwise use Legendre-Gauss quadrature in
    ! each hemisphere.  The maximum is 8.
    integer :: nangle = 0

    logical :: do_read_planck = .false.
    logical :: do_write_planck = .false.
    logical :: do_write_optical_depth = .false.
    logical :: do_write_spectral_fluxes = .false.
    logical :: do_write_spectral_boundary_fluxes = .false.
    logical :: input_planck_per_sterad = .false.

  contains
    procedure :: read => read_config_from_namelist
    procedure :: get_bands

  end type longwave_config_type

contains

  subroutine read_config_from_namelist(this, file_name)
    class(longwave_config_type), intent(inout) :: this
    character(*), intent(in)                   :: file_name

    integer :: iosopen ! Status after calling open

    character(len=N_MAX_NAME_LENGTH) :: optical_depth_name
    character(len=N_MAX_NAME_LENGTH) :: pressure_name
    character(len=N_MAX_NAME_LENGTH) :: temperature_name
    character(len=N_MAX_NAME_LENGTH) :: wavenumber_name
    character(len=N_MAX_NAME_LENGTH) :: planck_name 
    character(len=N_MAX_NAME_LENGTH) :: surf_emission_name
    character(len=N_MAX_NAME_LENGTH) :: surf_emissivity_name
    character(len=N_MAX_NAME_LENGTH) :: surf_temperature_name

    real(jprb), dimension(N_MAX_BANDS) :: band_wavenumber1, band_wavenumber2

    real(jprb) :: pressure_scaling, surf_temperature

    integer :: nblocksize, iverbose, nspectralstride, nangle

    logical :: do_read_planck, do_write_planck
    logical :: do_write_optical_depth
    logical :: do_write_spectral_fluxes, input_planck_per_sterad
    logical :: do_write_spectral_boundary_fluxes

    namelist /longwave_config/ optical_depth_name, pressure_name, &
         &  temperature_name, wavenumber_name, planck_name, &
         &  surf_emission_name, surf_emissivity_name, &
         &  pressure_scaling, nblocksize, surf_temperature, &
         &  iverbose, do_read_planck, do_write_planck, nangle, &
         &  do_write_spectral_boundary_fluxes, &
         &  do_write_optical_depth, do_write_spectral_fluxes, nspectralstride, &
         &  band_wavenumber1, band_wavenumber2, input_planck_per_sterad

    optical_depth_name = this%optical_depth_name
    pressure_name      = this%pressure_name
    temperature_name   = this%temperature_name
    wavenumber_name    = this%wavenumber_name
    planck_name        = this%planck_name
    surf_emission_name = this%surf_emission_name
    surf_emissivity_name = this%surf_emissivity_name
    surf_temperature_name = this%surf_temperature_name
    pressure_scaling   = this%pressure_scaling
    surf_temperature   = this%surf_temperature
    nblocksize         = this%nblocksize
    iverbose           = this%iverbose
    nspectralstride    = this%nspectralstride
    nangle             = this%nangle
    do_read_planck     = this%do_read_planck
    do_write_planck    = this%do_write_planck
    do_write_optical_depth = this%do_write_optical_depth
    do_write_spectral_fluxes = this%do_write_spectral_fluxes
    do_write_spectral_boundary_fluxes = this%do_write_spectral_boundary_fluxes
    band_wavenumber1   = -1.0_jprb
    band_wavenumber2   = -1.0_jprb
    input_planck_per_sterad = this%input_planck_per_sterad

    open(unit=10, iostat=iosopen, file=trim(file_name))
    if (iosopen /= 0) then
      stop 'Error: could not open namelist file'
    else
      read(unit=10, nml=longwave_config)
      close(unit=10)
    end if

    this%optical_depth_name = optical_depth_name
    this%pressure_name      = pressure_name
    this%temperature_name   = temperature_name
    this%wavenumber_name    = wavenumber_name
    this%planck_name        = planck_name
    this%surf_emission_name = surf_emission_name
    this%surf_emissivity_name = surf_emissivity_name
    this%surf_temperature_name = surf_temperature_name
    this%pressure_scaling   = pressure_scaling
    this%surf_temperature   = surf_temperature
    this%nblocksize         = nblocksize
    this%iverbose           = iverbose
    this%nspectralstride    = nspectralstride
    this%nangle             = nangle
    this%do_read_planck     = do_read_planck
    this%do_write_planck    = do_write_planck
    this%do_write_optical_depth = do_write_optical_depth
    this%do_write_spectral_fluxes = do_write_spectral_fluxes
    this%do_write_spectral_boundary_fluxes = do_write_spectral_boundary_fluxes
    this%band_wavenumber1 = band_wavenumber1
    this%band_wavenumber2 = band_wavenumber2
    this%input_planck_per_sterad = input_planck_per_sterad

  end subroutine read_config_from_namelist

  ! Get the band bounds and return the number of bands
  integer function get_bands(this, band_wavenumber1, band_wavenumber2)

    class(longwave_config_type),           intent(in)    :: this
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

end module longwave_config
