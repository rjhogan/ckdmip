module shortwave_fluxes

contains

  !---------------------------------------------------------------------
  ! Calculate the two-stream coefficients gamma1-gamma3 in the
  ! shortwave
  subroutine calc_two_stream_gammas_sw(ng, mu0, ssa, g, &
       &                               gamma1, gamma2, gamma3)

    use parkind1, only           : jprb

    implicit none

    integer, intent(in) :: ng
    ! Cosine of solar zenith angle, single scattering albedo and
    ! asymmetry factor:
    real(jprb), intent(in)                :: mu0
    real(jprb), intent(in),  dimension(:) :: ssa, g
    real(jprb), intent(out), dimension(:) :: gamma1, gamma2, gamma3

    real(jprb) :: factor

    integer    :: jg

    ! Zdunkowski "PIFM" (Zdunkowski et al., 1980; Contributions to
    ! Atmospheric Physics 53, 147-66)
    do jg = 1, ng
      factor = 0.75_jprb*g(jg)
      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
    end do

  end subroutine calc_two_stream_gammas_sw


  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering.
  subroutine calc_reflectance_transmittance_sw(ng, mu0, od, ssa, &
       &      gamma1, gamma2, gamma3, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir)

    use parkind1, only           : jprb, jprd

    implicit none

    integer, intent(in) :: ng

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od, ssa

    ! The three transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3

    ! The direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

    ! Transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng) :: trans_dir_dir

    real(jprd) :: gamma4, alpha1, alpha2, k_exponent, reftrans_factor
    real(jprd) :: exponential0 ! = exp(-od/mu0)
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprd) :: k_mu0, k_gamma3, k_gamma4
    real(jprd) :: k_2_exponential, od_over_mu0
    integer    :: jg

    do jg = 1, ng
      od_over_mu0 = max(od(jg) / mu0, 0.0_jprd)
      ! In the IFS this appears to be faster without testing the value
      ! of od_over_mu0:
      if (.true.) then
!      if (od_over_mu0 > 1.0e-6_jprd) then
        gamma4 = 1.0_jprd - gamma3(jg)
        alpha1 = gamma1(jg)*gamma4     + gamma2(jg)*gamma3(jg) ! Eq. 16
        alpha2 = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4    ! Eq. 17
        
        ! Note that if the minimum value is reduced (e.g. to 1.0e-24)
        ! then noise starts to appear as a function of solar zenith
        ! angle
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             &       1.0e-12_jprd)) ! Eq 18
        k_mu0 = k_exponent*mu0
        k_gamma3 = k_exponent*gamma3(jg)
        k_gamma4 = k_exponent*gamma4
        ! Check for mu0 <= 0!
        exponential0 = exp(-od_over_mu0)
        trans_dir_dir(jg) = exponential0
        exponential = exp(-k_exponent*od(jg))
        
        exponential2 = exponential*exponential
        k_2_exponential = 2.0_jprd * k_exponent * exponential
        
        if (k_mu0 == 1.0_jprd) then
          k_mu0 = 1.0_jprd - 10.0_jprd*epsilon(1.0_jprd)
        end if
        
        reftrans_factor = 1.0_jprd / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
        ! Meador & Weaver (1980) Eq. 25
        ref_diff(jg) = gamma2(jg) * (1.0_jprd - exponential2) * reftrans_factor
        
        ! Meador & Weaver (1980) Eq. 26
        trans_diff(jg) = k_2_exponential * reftrans_factor
        
        ! Here we need mu0 even though it wasn't in Meador and Weaver
        ! because we are assuming the incoming direct flux is defined
        ! to be the flux into a plane perpendicular to the direction of
        ! the sun, not into a horizontal plane
        reftrans_factor = mu0 * ssa(jg) * reftrans_factor / (1.0_jprd - k_mu0*k_mu0)
        
        ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
        ! exp(-k_exponent*od) in case of very high optical depths
        ref_dir(jg) = reftrans_factor &
             &  * ( (1.0_jprd - k_mu0) * (alpha2 + k_gamma3) &
             &     -(1.0_jprd + k_mu0) * (alpha2 - k_gamma3)*exponential2 &
             &     -k_2_exponential*(gamma3(jg) - alpha2*mu0)*exponential0)
        
        ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
        ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term representing direct
        ! unscattered transmittance.  
        trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4 + alpha1*mu0) &
            & - exponential0 &
            & * ( (1.0_jprd + k_mu0) * (alpha1 + k_gamma4) &
            &    -(1.0_jprd - k_mu0) * (alpha1 - k_gamma4) * exponential2) )

      else
        ! Low optical-depth limit; see equations 19, 20 and 27 from
        ! Meador & Weaver (1980)
        trans_diff(jg)     = 1.0_jprb - gamma1(jg) * od(jg)
        ref_diff(jg)       = gamma2(jg) * od(jg)
        trans_dir_diff(jg) = (1.0_jprb - gamma3(jg)) * ssa(jg) * od(jg)
        ref_dir(jg)        = gamma3(jg) * ssa(jg) * od(jg)
        trans_dir_dir(jg)  = 1.0_jprd - od_over_mu0
      end if
    end do
  
  end subroutine calc_reflectance_transmittance_sw
  
  !---------------------------------------------------------------------
  ! Shortwave adding method treating each column independently
  subroutine adding_ica_sw(ncol, nlev, incoming_toa, &
       &  albedo_surf, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn, flux_dn_direct)

    use parkind1, only           : jprb

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ncol)         :: albedo_surf

    ! Cosine of the solar zenith angle
    real(jprb), intent(in) :: cos_sza

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff

    ! Direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! total downwelling and direct downwelling (dimensioned ncol,nlev+1)
    real(jprb), intent(out), dimension(:,:) :: flux_up, flux_dn, flux_dn_direct
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(:,1) = incoming_toa
    do jlev = 1,nlev
      flux_dn_direct(:,jlev+1) = flux_dn_direct(:,jlev)*trans_dir_dir(:,jlev)
    end do

    albedo(:,nlev+1) = albedo_surf

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(:,nlev+1) = albedo_surf * flux_dn_direct(:,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      do jcol = 1,ncol
        ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev) * transmittance(jcol,jlev) &
             &                                     * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 11:
        source(jcol,jlev) = ref_dir(jcol,jlev)*flux_dn_direct(jcol,jlev) &
             &  + transmittance(jcol,jlev)*(source(jcol,jlev+1) &
             &        + albedo(jcol,jlev+1)*trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) &
             &  * inv_denominator(jcol,jlev)
      end do
    end do

    ! At the start, "flux_dn" is the diffuse downwelling.  At
    ! top-of-atmosphere there is no diffuse downwelling radiation.
    flux_dn(:,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(:,1) = source(:,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      do jcol = 1,ncol
        ! Shonk & Hogan (2008) Eq 14 (after simplification):
        flux_dn(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
        flux_dn_direct(jcol,jlev) = flux_dn_direct(jcol,jlev)*cos_sza
      end do
    end do
    flux_dn_direct(:,nlev+1) = flux_dn_direct(:,nlev+1)*cos_sza

    ! Add the direct to the diffuse downwelling to get the total
    flux_dn = flux_dn + flux_dn_direct

  end subroutine adding_ica_sw

  !---------------------------------------------------------------------
  ! Calculate shortwave flux profiles using the two-stream
  ! approximation
  subroutine calc_shortwave_fluxes(nlev, istartspec, iendspec, &
       &  cos_sza, ssi, albedo, &
       &  od, ssa, asymmetry, &
       &  flux_dn_direct, flux_dn, flux_up)

    use parkind1, only : jprb

    implicit none

    ! Number of height levels and spectral intervals
    integer, intent(in) :: nlev

    ! Range of spectral intervals to process
    integer, intent(in) :: istartspec, iendspec

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: cos_sza

    ! Spectral solar irradiance (W m-2)
    real(jprb), intent(in) :: ssi(:)

    ! Surface albedo
    real(jprb), intent(in) :: albedo(:)

    ! Optical depth, single scattering albedo and asymmetry factor
    real(jprb), intent(in) :: od(:,:)
    real(jprb), intent(in) :: ssa(:,:)
    real(jprb), intent(in) :: asymmetry(:,:)

    ! Output spectral fluxes down and up (W m-2)
    real(jprb), intent(inout), dimension(:,:) :: flux_dn_direct, flux_dn, flux_up

    ! Diffuse reflectance and transmittance for each layer
    real(jprb), dimension(istartspec:iendspec, nlev) :: reflectance, transmittance

    ! Fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams
    real(jprb), dimension(istartspec:iendspec, nlev) :: ref_dir, trans_dir_diff

    ! Transmittance for the direct beam in clear and all skies
    real(jprb), dimension(istartspec:iendspec, nlev) :: trans_dir_dir

    ! Two-stream coefficients
    real(jprb), dimension(istartspec:iendspec) :: gamma1, gamma2, gamma3

    integer :: jlev, nspec_local

    nspec_local = iendspec - istartspec + 1
    
    do jlev = 1,nlev
      call calc_two_stream_gammas_sw(nspec_local, cos_sza, &
           &  ssa(istartspec:iendspec,jlev), asymmetry(istartspec:iendspec,jlev), &
           &  gamma1, gamma2, gamma3)
      call calc_reflectance_transmittance_sw(nspec_local, &
           &  cos_sza, &
           &  od(istartspec:iendspec,jlev), ssa(istartspec:iendspec,jlev), &
           &  gamma1, gamma2, gamma3, &
           &  reflectance(:,jlev), transmittance(:,jlev), &
           &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
           &  trans_dir_dir(:,jlev) )
    end do

    ! Use adding method to compute fluxes
    call adding_ica_sw(nspec_local, nlev, ssi(istartspec:iendspec), &
         &  albedo(istartspec:iendspec), cos_sza, reflectance, transmittance, ref_dir, &
         &  trans_dir_diff, trans_dir_dir, flux_up(istartspec:iendspec,:), &
         &  flux_dn(istartspec:iendspec,:), flux_dn_direct(istartspec:iendspec,:))
    
  end subroutine calc_shortwave_fluxes

end module shortwave_fluxes
