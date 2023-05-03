module planck_function

contains

  subroutine calc_planck_function(nspec,nlev,istartspec,iendspec, &
       &  temperature, wavenumber_cm1, planck)
    
    use parkind1, only : jprb

    implicit none

    real(jprb), parameter :: PLANCK_CONSTANT = 6.62606896e-34_jprb
    real(jprb), parameter :: SPEED_OF_LIGHT = 2.99792458e8_jprb
    real(jprb), parameter :: BOLTZMANN_CONSTANT =  1.3806504e-23_jprb
    real(jprb), parameter :: INV_CM2_HZ = 100.0_jprb * SPEED_OF_LIGHT
    real(jprb), parameter :: PI = 3.14159265358979323846_jprb

    ! Number of spectral intervals and levels
    integer, intent(in) :: nspec, nlev

    ! Start and end spectral interval
    integer, intent(in) :: istartspec, iendspec

    ! Temperature, K
    real(jprb), intent(in) :: temperature(nlev)

    ! Wavenumber, cm-1
    real(jprb), intent(in) :: wavenumber_cm1(nspec)

    ! Planck function, W m-2
    real(jprb), intent(out) :: planck(nspec,nlev)

    ! Wavenumber interval, cm-1
    real(jprb) :: dwav(istartspec:iendspec)

    real(jprb) :: prefactor, freq_Hz

    integer :: jlev, jspec

    dwav(istartspec:iendspec-1) &
         &  = abs(wavenumber_cm1(istartspec+1:iendspec) &
         &       -wavenumber_cm1(istartspec:iendspec-1))
    dwav(iendspec) = dwav(iendspec-1)

    prefactor = 2.0_jprb * PLANCK_CONSTANT * INV_CM2_HZ * PI &
         &  / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)

    do jlev = 1,nlev
      do jspec = istartspec,iendspec
        freq_Hz = wavenumber_cm1(jspec) * INV_CM2_HZ
        planck(jspec,jlev) = prefactor * dwav(jspec) * (freq_Hz**3) &
             &  / (exp((PLANCK_CONSTANT/BOLTZMANN_CONSTANT) &
             &        *(freq_Hz/temperature(jlev))) - 1.0_jprb)
        
      end do
    end do

  end subroutine calc_planck_function

  subroutine calc_planck_function1(nspec,istartspec,iendspec, &
       &  temperature, wavenumber_cm1, planck)
    
    use parkind1, only : jprb

    implicit none

    real(jprb), parameter :: PLANCK_CONSTANT = 6.62606896e-34_jprb
    real(jprb), parameter :: SPEED_OF_LIGHT = 2.99792458e8_jprb
    real(jprb), parameter :: BOLTZMANN_CONSTANT =  1.3806504e-23_jprb
    real(jprb), parameter :: INV_CM2_HZ = 100.0_jprb * SPEED_OF_LIGHT
    real(jprb), parameter :: PI = 3.14159265358979323846_jprb

    ! Number of spectral intervals
    integer, intent(in) :: nspec

    ! Start and end spectral interval
    integer, intent(in) :: istartspec, iendspec

    ! Temperature, K
    real(jprb), intent(in) :: temperature

    ! Wavenumber, cm-1
    real(jprb), intent(in) :: wavenumber_cm1(nspec)

    ! Planck function, W m-2
    real(jprb), intent(out) :: planck(nspec)

    ! Wavenumber interval, cm-1
    real(jprb) :: dwav(istartspec:iendspec)

    real(jprb) :: prefactor, freq_Hz

    integer :: jspec

    dwav(istartspec:iendspec-1) &
         &  = abs(wavenumber_cm1(istartspec+1:iendspec) &
         &       -wavenumber_cm1(istartspec:iendspec-1))
    dwav(iendspec) = dwav(iendspec-1)

    prefactor = 2.0_jprb * PLANCK_CONSTANT * INV_CM2_HZ * PI &
         &  / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)

    do jspec = istartspec,iendspec
      freq_Hz = wavenumber_cm1(jspec) * INV_CM2_HZ
      planck(jspec) = prefactor * dwav(jspec) * (freq_Hz**3) &
           &  / (exp((PLANCK_CONSTANT/BOLTZMANN_CONSTANT) &
           &        *(freq_Hz/temperature)) - 1.0_jprb)
      
    end do

  end subroutine calc_planck_function1

end module planck_function
