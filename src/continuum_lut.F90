! continuum_lut.F90 - Compute continuum absorption from a NetCDF look-up table

module continuum_lut

  use parkind1,      only : jprb

  implicit none

  public
  
  type continuum_data
    ! Number of wavenumbers
    integer :: nwav = 0
    ! Wavenumbers (cm-1)
    real(jprb), allocatable :: wavenumber(:)
    ! Self and foreign continuum absorption coefficients at reference
    ! temperature (cm2 molec-1 cm-1)
    real(jprb), allocatable :: self_abs_coeff_ref(:)
    real(jprb), allocatable :: foreign_abs_coeff_ref(:)
    ! Self continuum temperature absorption
    real(jprb), allocatable :: self_temperature_exponent(:)
    ! Reference pressure (Pa)
    real(jprb) :: pressure_ref
    ! Reference temperature (K)
    real(jprb) :: temperature_ref
    ! Do we multiply by the radiation term?
    logical :: use_radiation_term = .true.

  contains
    procedure :: read => read_lut
    procedure :: calc => calc_continuum
    procedure :: is_available
    
  end type continuum_data
  
contains

  function is_available(self)
    class(continuum_data), intent(in) :: self
    logical :: is_available
    is_available = (self%nwav > 0)
  end function is_available

  
  subroutine read_lut(self, file_name)
    
    use easy_netcdf, only : netcdf_file

    class(continuum_data), intent(inout) :: self
    character(*),                   intent(in)    :: file_name
    type(netcdf_file) :: file

    integer :: n_use_radiation_term
    
    call file%open(file_name)

    call file%get('wavenumbers', self%wavenumber)
    self%nwav = size(self%wavenumber)
    call file%get('self_absco_ref', self%self_abs_coeff_ref)
    call file%get('for_absco_ref',  self%foreign_abs_coeff_ref)
    call file%get('self_texp',      self%self_temperature_exponent)
    call file%get('ref_press',      self%pressure_ref)
    self%pressure_ref = self%pressure_ref * 100.0_jprb
    call file%get('ref_temp',       self%temperature_ref)

    self%use_radiation_term = .true.
    if (file%exists('use_radiation_term')) then
      call file%get('use_radiation_term', n_use_radiation_term)
      if (n_use_radiation_term == 0) then
        self%use_radiation_term = .false.
      end if
    end if
    
    call file%close()
    
  end subroutine read_lut
  
  subroutine calc_continuum(self, nlev, nwav, pressure, temperature, mole_fraction, &
       wavenumber_cm1, continuum)

    use interpolation, only : interpolate

    class(continuum_data), intent(in) :: self
    
    ! Avogadro's number
    real(jprb), parameter :: NAVOGADRO = 6.02214076e23_jprb

    ! Number of levels and wavenumbers
    integer, intent(in) :: nlev, nwav
    
    ! Pressure (Pa), temperature (K) and water vapour mole fraction
    ! (mol/mol)
    real(jprb), intent(in) :: pressure(nlev), temperature(nlev), mole_fraction(nlev)
    
    ! Wavenumber (cm-1)
    real(jprb), intent(in) :: wavenumber_cm1(nwav)

    ! Continuum molar absorption (m2 mol-1)
    real(jprb), intent(out) :: continuum(nwav,nlev)

    ! Continuum molar absorption (m2 mol-1) at spectral resolution of
    ! look-up table
    real(jprb) :: continuum_coarse(self%nwav)

    real(jprb) :: radiation_term(self%nwav)

    real(jprb) :: rho_ratio, scaling
    
    ! Loop index for level
    integer :: jlev

    radiation_term = 1.0_jprb
    
    do jlev = 1,nlev
      if (self%use_radiation_term) then
        call calc_radiation_term(self%nwav, self%wavenumber, temperature(jlev), radiation_term)
      end if
      
      rho_ratio = (pressure(jlev)/self%pressure_ref) &
           &    * (self%temperature_ref/temperature(jlev))
      scaling   = rho_ratio * 0.0001_jprb * NAVOGADRO
      
      ! Self continuum
      continuum_coarse = self%self_abs_coeff_ref * radiation_term &
           &  * (self%temperature_ref/temperature(jlev))**self%self_temperature_exponent &
           &  * (mole_fraction(jlev) * scaling)
      ! Add foreign continuum
      continuum_coarse = continuum_coarse + self%foreign_abs_coeff_ref &
           &  * radiation_term * ((1.0_jprb - mole_fraction(jlev)) * scaling)
      
      ! Interpolate to the input wavenumber
      call interpolate(self%wavenumber, continuum_coarse, wavenumber_cm1, continuum(:,jlev))
      
    end do

  end subroutine calc_continuum

  
  subroutine calc_radiation_term(nwav, wavenumber, temperature, radiation_term)

    ! Radiation constant: Planck * SpeedOfLight / Boltzmann
    real(jprb), parameter :: RADCN2 = 1.4387752_jprb
    
    integer,    intent(in)  :: nwav
    real(jprb), intent(in)  :: wavenumber(nwav) ! cm-1
    real(jprb), intent(in)  :: temperature ! K
    real(jprb), intent(out) :: radiation_term(nwav)

    real(jprb) :: xkt
    real(jprb) :: xviokt(nwav), expvkt(nwav)

    xkt = temperature / RADCN2
    xviokt = wavenumber / xkt

    where (xviokt <= 0.01_jprb)
      radiation_term = 0.5_jprb * xviokt * wavenumber
    elsewhere (xviokt <= 10.0_jprb)
      expvkt = exp(-xviokt)
      radiation_term = wavenumber * (1.0_jprb-expvkt)/(1.0_jprb+expvkt)
    elsewhere
      radiation_term = wavenumber
    end where
    
  end subroutine calc_radiation_term

  
end module continuum_lut
