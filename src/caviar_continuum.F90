! caviar_continuum.F90 - Compute continuum absorption using CAVIAR model

module caviar_continuum

public
  
contains
  
  subroutine calc_caviar_continuum(nlev, nspec, pressure, temperature, mole_fraction, &
       wavenumber_cm1, continuum)

    use parkind1, only : jprb

    ! Number of levels and wavenumbers
    integer, intent(in) :: nlev, nspec

    ! Pressure (Pa), temperature (K) and water vapour mole fraction
    ! (mol/mol)
    real(jprb), intent(in) :: pressure(nlev), temperature(nlev), mole_fraction(nlev)
    
    ! Wavenumber (cm-1)
    real(jprb), intent(in) :: wavenumber_cm1(nspec)

    ! Continuum molar absorption (m2 mol-1)
    real(jprb), intent(out) :: continuum(nspec,nlev)

    ! Insert calculation of water vapour continuum here...
    continuum = 0.0_jprb
    
  end subroutine calc_caviar_continuum
  
end module caviar_continuum
