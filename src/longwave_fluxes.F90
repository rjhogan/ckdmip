module longwave_fluxes

  integer, parameter :: MAX_ANGLES = 16

contains

  ! Compute upwelling and downwelling longwave spectral fluxes using a
  ! single zenith angle in each hemisphere of +/-53 degrees (diffusivity 1.66)
  subroutine calc_longwave_fluxes(nlev,nspec,istartspec,iendspec,od,planck_hl, &
       &  surf_emissivity, surf_emission, flux_dn, flux_up)

    use parkind1, only : jprb

    implicit none

    real(jprb), parameter :: DIFFUSIVITY = 1.66_jprb
    !real(jprb), parameter :: DIFFUSIVITY = 2.0_jprb

    ! Number of height levels and spectral intervals
    integer, intent(in) :: nlev, nspec

    ! Range of spectral intervals to process
    integer, intent(in) :: istartspec, iendspec

    ! Optical depth
    real(jprb), intent(in) :: od(nspec,nlev)

    ! Planck function at half-levels, W m-2
    real(jprb), intent(in) :: planck_hl(nspec,nlev+1)

    ! Surface emissivity
    real(jprb), intent(in) :: surf_emissivity(nspec)

    ! Surface emission, W m-2
    real(jprb), intent(in) :: surf_emission(nspec)

    ! Output spectral fluxes down and up (W m-2)
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: flux_dn, flux_up

    real(jprb) :: emissivity(istartspec:iendspec,nlev)
    real(jprb) :: factor(istartspec:iendspec,nlev)

    integer :: jlev

    emissivity = 1.0_jprb - exp(-DIFFUSIVITY*od(istartspec:iendspec,:))
    where (emissivity > 1.0e-5)
      factor = 1.0_jprb - emissivity * (1.0_jprb/DIFFUSIVITY) &
           &  / od(istartspec:iendspec,:)
    elsewhere
      factor = 0.5 * emissivity
    end where

    ! Work down from top of atmosphere
    flux_dn(istartspec:iendspec,1) = 0.0_jprb
    do jlev = 1,nlev
      flux_dn(istartspec:iendspec,jlev+1) = flux_dn(istartspec:iendspec,jlev) &
           &  * (1.0_jprb - emissivity(:,jlev)) &
           &  + planck_hl(istartspec:iendspec,jlev)*(emissivity(:,jlev)-factor(:,jlev)) &
           &  + planck_hl(istartspec:iendspec,jlev+1) * factor(:,jlev)
    end do

    ! Work up from surface
    flux_up(istartspec:iendspec,nlev+1) = surf_emission(istartspec:iendspec) &
         &  + (1.0_jprb - surf_emissivity(istartspec:iendspec)) * flux_dn(istartspec:iendspec,nlev+1)
    do jlev = nlev,1,-1
      flux_up(istartspec:iendspec,jlev) = flux_up(istartspec:iendspec,jlev+1) &
           &  * (1.0_jprb - emissivity(:,jlev)) &
           &  + planck_hl(istartspec:iendspec,jlev+1)*(emissivity(:,jlev)-factor(:,jlev)) &
           &  + planck_hl(istartspec:iendspec,jlev) * factor(:,jlev)
    end do

  end subroutine calc_longwave_fluxes

  ! As the function above, but using "nang" angles in each hemisphere
  ! spaced in cos(angle) according to Legendre-Gauss quadrature.  If
  ! nang=0 then call the function above.  If nang is negative then
  ! adjust the quadrature points to favour slightly smaller zenith
  ! angles.
  subroutine calc_longwave_fluxes_n(nang, nlev,nspec,istartspec,iendspec, &
       &  od,planck_hl, surf_emissivity, surf_emission, flux_dn, flux_up)

    use parkind1, only : jprb

    implicit none

    ! Number of angles, height levels and spectral intervals
    integer, intent(in) :: nang, nlev, nspec

    ! Range of spectral intervals to process
    integer, intent(in) :: istartspec, iendspec

    ! Optical depth
    real(jprb), intent(in) :: od(nspec,nlev)

    ! Planck function at half-levels, W m-2
    real(jprb), intent(in) :: planck_hl(nspec,nlev+1)

    ! Surface emissivity
    real(jprb), intent(in) :: surf_emissivity(nspec)

    ! Surface emission, W m-2
    real(jprb), intent(in) :: surf_emission(nspec)

    ! Output spectral fluxes down and up (W m-2)
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: flux_dn, flux_up

    real(jprb), dimension(MAX_ANGLES) :: mu, weight

    real(jprb) :: flux_before(istartspec:iendspec), flux_after(istartspec:iendspec)

    real(jprb) :: emissivity(istartspec:iendspec)
    real(jprb) :: factor(istartspec:iendspec)

    real(jprb) :: secant, planck_scaling

    integer :: jlev, jang, nang_local

    if (nang == 0) then
      ! Call classic 2-stream (1-angle) model
      call calc_longwave_fluxes(nlev,nspec,istartspec,iendspec,od,planck_hl, &
           &  surf_emissivity, surf_emission, flux_dn, flux_up)
      return
    end if

    ! Replace negative values with positive
    nang_local = min(abs(nang), MAX_ANGLES)
    
    mu = 0.0_jprb
    weight = 0.0_jprb
    if (nang > 0) then
      ! Positive nang leads to classic double-Gauss (Sykes 1951),
      ! which involves Gauss-Legendre quadrature in each hemisphere
      if (nang_local == 1) then
        mu(1:1)     = [0.5000000000]
        weight(1:1) = [1.0000000000]
      else if (nang_local == 2) then
        mu(1:2)     = [0.7886751346, 0.2113248654]
        weight(1:2) = [0.5000000000, 0.5000000000]
      else if (nang_local == 3) then
        mu(1:3)     = [0.8872983346, 0.5000000000, 0.1127016654]
        weight(1:3) = [0.2777777778, 0.4444444444, 0.2777777778]
      else if (nang_local == 4) then
        mu(1:4)     = [0.9305681558, 0.6699905218, 0.3300094782, 0.0694318442]
        weight(1:4) = [0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226]
      else if (nang_local == 5) then
        mu(1:5)     = [0.9530899230, 0.7692346551, 0.5000000000, 0.2307653449, 0.0469100770]
        weight(1:5) = [0.1184634425, 0.2393143352, 0.2844444444, 0.2393143352, 0.1184634425]
      else if (nang_local == 6) then
        mu(1:6)     = [0.9662347571, 0.8306046932, 0.6193095930, 0.3806904070, 0.1693953068, &
             &         0.0337652429]
        weight(1:6) = [0.0856622462, 0.1803807865, 0.2339569673, 0.2339569673, 0.1803807865, &
             &         0.0856622462]
      else if (nang_local == 7) then
        mu(1:7)     = [0.9745539562, 0.8707655928, 0.7029225757, 0.5000000000, 0.2970774243, &
             &         0.1292344072, 0.0254460438]
        weight(1:7) = [0.0647424831, 0.1398526957, 0.1909150253, 0.2089795918, 0.1909150253, &
             &         0.1398526957, 0.0647424831]
      else
        mu(1:8)     = [0.9801449282, 0.8983332387, 0.7627662050, 0.5917173212, 0.4082826788, &
             &         0.2372337950, 0.1016667613, 0.0198550718]
        weight(1:8) = [0.0506142681, 0.1111905172, 0.1568533229, 0.1813418917, 0.1813418917, &
             &         0.1568533229, 0.1111905172, 0.0506142681]
      end if
    else
#define GAUSS_JACOBI_3 1
#ifdef GAUSS_LAGUERRE
      ! Negative nang corresponds to Gauss-Laguerre quadrature in each
      ! hemisphere after a change of variables
      if (nang_local == 1) then
        mu(1:1) = [0.6065306597]
        weight(1:1) = [1.0000000000]
      else if (nang_local == 2) then
        mu(1:2) = [0.7461018061, 0.1813898346]
        weight(1:2) = [0.8535533906, 0.1464466094]
      else if (nang_local == 3) then
        mu(1:3) = [0.8122985952, 0.3175435896, 0.0430681066]
        weight(1:3) = [0.7110930099, 0.2785177336, 0.0103892565]
      else if (nang_local == 4) then
        mu(1:4) = [0.8510589811, 0.4177464746, 0.1034869099, 0.0091177205]
        weight(1:4) = [0.6031541043, 0.3574186924, 0.0388879085, 0.0005392947]
      else if (nang_local == 5) then
        mu(1:5) = [0.8765336712, 0.4932685488, 0.1655945604, 0.0289291656, 0.0017992229]
        weight(1:5) = [0.5217556106, 0.3986668111, 0.0759424497, 0.0036117587, 0.0000233700]
      else if (nang_local == 6) then
        mu(1:6) = [0.8945599997, 0.5518571509, 0.2239420059, &
             &     0.0557113276, 0.0073083795, 0.0003383475]
        weight(1:6) = [0.4589646739, 0.4170008308, 0.1133733821, &
             &         0.0103991975, 0.0002610172, 0.0000008985]
      else if (nang_local == 7) then
        mu(1:7) = [0.9079900684, 0.5984977893, 0.2769444395, 0.0862783534, &
             &     0.0167212198, 0.0017171486, 0.0000614145]
        weight(1:7) = [0.4093189517, 0.4218312779, 0.1471263487, 0.0206335145, &
             &         0.0010740101, 0.0000158655, 0.0000000317]
      else if (nang_local == 8) then
        mu(1:8) = [0.9183838705, 0.6364490646, 0.3244761267, 0.1184398449, &
             &     0.0295121658, 0.0046112422, 0.0003819048, 0.0000108476]
        weight(1:8) = [0.3691885893, 0.4187867808, 0.1757949866, 0.0333434923, &
             &         0.0027945362, 0.0000907651, 0.0000008486, 0.0000000010]
      else
        ! 16 angles for reference calculations
        mu(1:16) = [9.57121721e-01, 7.93463162e-01, 5.65226418e-01, 3.44851351e-01, &
             &      1.79327181e-01, 7.89445709e-02, 2.91538218e-02, 8.92269571e-03, &
             &      2.22697374e-03, 4.43521764e-04, 6.84040566e-05, 7.82683106e-06, &
             &      6.22597950e-07, 3.09250485e-08, 7.81171559e-10, 5.93247174e-12]
        weight(1:16) = [2.06151715e-01, 3.31057855e-01, 2.65795778e-01, 1.36296934e-01, &
             &          4.73289287e-02, 1.12999001e-02, 1.84907094e-03, 2.04271915e-04, &
             &          1.48445869e-05, 6.82831933e-07, 1.88102484e-08, 2.86235024e-10, &
             &          2.12707903e-12, 6.29796700e-15, 5.05047370e-18, 4.16146237e-22]
      end if
#elif defined(GAUSS_JACOBI_2)
      ! Gauss-Jacobi weighted by mu^2
      if (nang_local == 1) then
        mu(1:1) = [0.6297376093]
        weight(1:1) = [1.0000000000]
      else if (nang_local == 2) then
        mu(1:2) = [0.2509907356, 0.7908473988]
        weight(1:2) = [0.2300253764, 0.7699746236]
        ! mu^3 weight
        !mu(1:2) = [0.2349482374, 0.7812446456]
        !weight(1:2) = [0.2101393363, 0.7898606637]
        ! mu^4 weight
        !mu(1:2) = [0.2249148266, 0.7750368438]
        !weight(1:2) = [0.1978663403, 0.8021336597]
        ! mu^5 weight
        !mu(1:2) = [0.2180469863, 0.7706935658]
        !weight(1:2) = [0.1895464137, 0.8104535863]
      else if (nang_local == 3) then
        mu(1:3) = [0.1024922169, 0.4417960320, 0.8633751621]
        weight(1:3) = [0.0437820218, 0.3875796738, 0.5686383044]
      else if (nang_local == 4) then
        mu(1:4) = [0.0454586727, 0.2322334416, 0.5740198775, 0.9030775973]
        weight(1:4) = [0.0092068785, 0.1285704278, 0.4323381850, 0.4298845087]
      else if (nang_local == 5) then
        mu(1:5) = [0.0218799020, 0.1240348649, 0.3523952142, 0.6662598336, 0.9274122772]
        weight(1:5) = [0.0022182930, 0.0403781426, 0.2026064698, 0.4204280382, 0.3343690564]
      else if (nang_local == 6) then
        mu(1:6) = [0.0113101689, 0.0686859364, 0.2136209604, &
             &     0.4537170880, 0.7322872118, 0.9434907860]
        weight(1:6) = [0.0006079549, 0.0131275539, 0.0838179185, &
             &         0.2489180819, 0.3867531550, 0.2667753358]
      else if (nang_local == 7) then
        mu(1:7) = [0.0062121285, 0.0395847398, 0.1310504716, 0.3012526699, &
             &     0.5365521763, 0.7808646396, 0.9547040091]
        weight(1:7) = [0.0001866276, 0.0045323030, 0.0339624822, 0.1257557891, &
             &         0.2706298895, 0.3474456812, 0.2174872274]
      else if (nang_local == 8) then
        mu(1:8) = [0.0035921679, 0.0237052993, 0.0821015442, 0.1997077560, &
             &     0.3812774921, 0.6037571134, 0.8175144862, 0.9628512704]
        weight(1:8) = [0.0000631899, 0.0016701517, 0.0140120490, 0.0602676863, &
             &         0.1589118067, 0.2753139192, 0.3092053430, 0.1805558542]
      else
        mu(1:16) = [0.0001457809, 0.0010714784, 0.0042566315, 0.0122621677, &
             &      0.0286976381, 0.0578386189, 0.1039041907, 0.1701025130, &
             &      0.2576328357, 0.3648726840, 0.4869664791, 0.6159660528, &
             &      0.7415668747, 0.8523592307, 0.9373996714, 0.9878422725]
        weight(1:16) = [0.0000001081, 0.0000036923, 0.0000432005, 0.0002831227, &
             &          0.0012656394, 0.0042743008, 0.0115817439, 0.0261542136, &
             &          0.0504546422, 0.0844669595, 0.1237887203, 0.1591145057, &
             &          0.1781461152, 0.1698212979, 0.1292922091, 0.0613095288]
      end if
#elif defined(GAUSS_JACOBI_3)
      ! Gauss-Jacobi weighted by mu^3
      if (nang_local == 1) then
        mu(1:1) = [0.6242950770]
        weight(1:1) = [1.0000000000]
      else if (nang_local == 2) then
        mu(1:2) = [0.2349482374, 0.7812446456]
        weight(1:2) = [0.2101393363, 0.7898606637]
      else if (nang_local == 3) then
        mu(1:3) = [0.0874707892, 0.4151791875, 0.8534188521]
        weight(1:3) = [0.0337794282, 0.3671025488, 0.5991180230]
      else if (nang_local == 4) then
        mu(1:4) = [0.0347818334, 0.2022426456, 0.5433165128, 0.8938137691]
        weight(1:4) = [0.0057725036, 0.1060509098, 0.4246726272, 0.4635039594]
      else if (nang_local == 5) then
        mu(1:5) = [0.0149092234, 0.0986620567, 0.3125617509, 0.6350459475, 0.9190800496]
        weight(1:5) = [0.0011118620, 0.0280708244, 0.1768871544, 0.4263324843, 0.3675976749]
      else if (nang_local == 6) then
        mu(1:6) = [0.0068546680, 0.0495503941, 0.1753082863, 0.4087850488, 0.7022620909, 0.9360772414]
        weight(1:6) = [0.0002424788, 0.0075640894, 0.0632176307, 0.2279083544, 0.4029301896, 0.2981372571]
      else if (nang_local == 7) then
        mu(1:7) = [0.0033536451, 0.0258276321, 0.0988138915, 0.2536346871, 0.4898128245, 0.7527432969, 0.9481207154]
        weight(1:7) = [0.0000593269, 0.0021488443, 0.0217770553, 0.1013633762, 0.2578344957, 0.3703483283, 0.2464685734]
      else if (nang_local == 8) then
        mu(1:8) = [0.0017324901, 0.0139859972, 0.0566875572, 0.1560660871, &
             &     0.3278765825, 0.5572544978, 0.7915221980, 0.9569945526]
        weight(1:8) = [0.0000160890, 0.0006504888, 0.0075758567, 0.0421198924, &
             &         0.1353164566, 0.2711932428, 0.3360480756, 0.2070798980]
      else
        mu(1:16) = [0.0000326885, 0.0003073222, 0.0014975432, 0.0051369730, &
             &      0.0139811654, 0.0321263133, 0.0646812661, 0.1168909532, &
             &      0.1927872045, 0.2936423434, 0.4166485035, 0.5542763449, &
             &      0.6946500175, 0.8230297196, 0.9241861162, 0.9851895899]
        weight(1:16) = [0.0000000060, 0.0000003457, 0.0000061909, 0.0000582387, &
             &          0.0003553366, 0.0015709036, 0.0053760948, 0.0148591346, &
             &          0.0341109566, 0.0662374608, 0.1099664363, 0.1566225998, &
             &          0.1902844545, 0.1929330102, 0.1532362599, 0.0743825708]
      end if
#else
      ! Gauss-Jacobi weighted by mu^4
      if (nang_local == 1) then
        mu(1:1) = [0.6209213231]
        weight(1:1) = [1.0000000000]
      else if (nang_local == 2) then
        mu(1:2) = [0.2249148266, 0.7750368438]
        weight(1:2) = [0.1978663403, 0.8021336597]
      else if (nang_local == 3) then
        mu(1:3) = [0.0784317576, 0.3979502136, 0.8467120319]
        weight(1:3) = [0.0282265034, 0.3528809560, 0.6188925406]
      else if (nang_local == 4) then
        mu(1:4) = [0.0287603127, 0.1834478492, 0.5226891912, 0.8873514770]
        weight(1:4) = [0.0041364167, 0.0922355842, 0.4173953904, 0.4862326087]
      else if (nang_local == 5) then
        mu(1:5) = [0.0112768811, 0.0835879686, 0.2865680766, 0.6133656689, 0.9130934048]
        weight(1:5) = [0.0006708555, 0.0214648056, 0.1593367688, 0.4276816587, 0.3908459113]
      else if (nang_local == 6) then
        mu(1:6) = [0.0047254145, 0.0389089797, 0.1514424378, 0.3783879142, 0.6807956269, 0.9306149581]
        weight(1:6) = [0.0001221237, 0.0049989268, 0.0509126775, 0.2116666960, 0.4115766859, 0.3207228902]
      else if (nang_local == 7) then
        mu(1:7) = [0.0021051025, 0.0187194658, 0.0798570837, 0.2227216609, 0.4571845590, 0.7321283708, 0.9431640320]
        weight(1:7) = [0.0000248683, 0.0012156973, 0.0154356905, 0.0852618670, 0.2458799433, 0.3843213601, 0.2678605736]
      else if (nang_local == 8) then
        mu(1:8) = [0.0009908313, 0.0093390210, 0.0426771277, 0.1291466768, &
             &     0.2919494780, 0.5238946383, 0.7720490819, 0.9525012389]
        weight(1:8) = [0.0000056157, 0.0003135776, 0.0046812155, 0.0316256515, &
             &         0.1181113268, 0.2646092424, 0.3535749056, 0.2270784649]
      else
        mu(1:16) = [0.0000096083, 0.0001102151, 0.0006331044, 0.0025012328, &
             &      0.0076975959, 0.0196928664, 0.0435503543, 0.0854111740, &
             &      0.1512042671, 0.2447132168, 0.3654702650, 0.5071893381, &
             &      0.6574558043, 0.7991033388, 0.9131931340, 0.9829586317]
        weight(1:16) = [0.0000000006, 0.0000000488, 0.0000012282, 0.0000154523, &
             &          0.0001213096, 0.0006680949, 0.0027695777, 0.0090443184, &
             &          0.0239833922, 0.0526916350, 0.0970782165, 0.1506688823, &
             &          0.1960367989, 0.2093378910, 0.1722896229, 0.0852935307]
      end if
#endif
      ! Adjust the weights so that sum(mu*weight)==0.5, including
      ! dividing by mu since the weights will be multiplied by mu
      ! later
      weight(1:nang_local) = weight(1:nang_local) / mu(1:nang_local)
      weight = weight * 0.5_jprb / sum(mu*weight)
    end if

    ! ! If nang is negative then perform an adjustment to favour smaller
    ! ! zenith angles
    ! if (nang == -1) then
    !   mu(1:1) = 1.0_jprb / 1.6238_jprb
    !   weight(1:1) = 0.5_jprb / mu(1:1)
    ! else if (nang == -2) then
    !   ! Optimized to minimize error in transmission across a wide
    !   ! range of path lengths
    !   mu(1:2) = [0.8018_jprb, 0.2704_jprb]
    !   ! ...and forced to be a factor of 3 apart for optimization
    !   !mu(1:2) = [0.8034_jprb, 0.2678_jprb]*0.99_jprb
    !   weight = weight * 0.5_jprb / sum(mu*weight)
    ! else if (nang < 0) then
    !   !mu = (2.0_jprb/acos(-1.0_jprb))*(sqrt(1.0_jprb-mu*mu)*mu + asin(mu))
    !   mu(1:nang_local) = 1.0_jprb - 0.96_jprb * (1.0_jprb*mu(1:nang_local))

    !   !mu(1) = 1.0_jprb/1.66_jprb
    !   ! Adjust the weights so that sum(mu*weight)==0.5
    !   weight = weight * 0.5_jprb / sum(mu*weight)
    ! end if

    ! Work down from top of atmosphere
    flux_dn(istartspec:iendspec,:) = 0.0_jprb
    
    ! Loop over angles, adding contribution from each angle to the total
    do jang = 1,nang_local

      secant = 1.0_jprb / mu(jang)

      ! Planck function coming in, planck_hl, is emission through a
      ! horizontal surface in W m-2.  Radiance calculations would use
      ! the Planck function per steradian, which is a factor of 1/PI
      ! times this.  To integrate radiances over a hemisphere to get
      ! the flux through a horizontal surface, we integrate over
      ! azimuth (factor of 2*PI) and over mu (from 0 to 1 dealt with
      ! by Legendre-Gauss quadrature) and multiply by mu because
      ! glancing angles contribute less to a flux than zenith/nadir
      ! angles. Thus we have 2*mu.
      planck_scaling = 2.0_jprb * mu(jang) * weight(jang)

      ! Zero flux coming in at top of atmosphere
      flux_after = 0.0_jprb

      ! Down through atmosphere
      do jlev = 1,nlev

        emissivity = 1.0_jprb - exp(-secant*od(istartspec:iendspec,jlev))
        ! Below 1.0e-5 the full calculation is subject to errors,
        ! which are large below 1.0e-6, corrected by the following
        where (emissivity > 1.0e-5)
          factor = 1.0_jprb - emissivity * mu(jang) / od(istartspec:iendspec,jlev)
        elsewhere
          factor = 0.5 * emissivity
        end where

        flux_before = flux_after
        flux_after = flux_before * (1.0_jprb - emissivity) &
             &  + (planck_hl(istartspec:iendspec,jlev)*(emissivity-factor) &
             &  +  planck_hl(istartspec:iendspec,jlev+1) * factor) * planck_scaling
        flux_dn(istartspec:iendspec,jlev+1) = flux_dn(istartspec:iendspec,jlev+1) &
             &  + flux_after
      end do

    end do

    ! Work up from surface - first zero the upwelling fluxes
    flux_up(istartspec:iendspec,:) = 0.0_jprb

    ! Surface reflection and emission
    flux_up(istartspec:iendspec,nlev+1) = surf_emission(istartspec:iendspec) &
         &  + (1.0_jprb - surf_emissivity(istartspec:iendspec)) &
         &  * flux_dn(istartspec:iendspec,nlev+1)

    ! Loop over angles, adding contribution from each angle to the total
    do jang = 1,nang_local

      secant = 1.0_jprb / mu(jang)

      planck_scaling = 2.0_jprb * mu(jang) * weight(jang)

      ! Assuming Lambertian reflection (if there is any reflection),
      ! the surface emission at this angle is simply the flux times
      ! the scaling, remembering that we include the weight by which
      ! this will contribute to the flux at all heights above
      flux_after = flux_up(istartspec:iendspec,nlev+1) * planck_scaling

      ! Up through atmosphere
      do jlev = nlev,1,-1

        emissivity = 1.0_jprb - exp(-secant*od(istartspec:iendspec,jlev))
        where (emissivity > 1.0e-5)
          factor = 1.0_jprb - emissivity * mu(jang) / od(istartspec:iendspec,jlev)
        elsewhere
          factor = 0.5 * emissivity
        end where

        flux_before = flux_after
        flux_after = flux_before * (1.0_jprb - emissivity) &
             &  + (planck_hl(istartspec:iendspec,jlev+1)*(emissivity-factor) &
             &  +  planck_hl(istartspec:iendspec,jlev) * factor) * planck_scaling
        flux_up(istartspec:iendspec,jlev) = flux_up(istartspec:iendspec,jlev) &
             &  + flux_after
      end do

    end do

  end subroutine calc_longwave_fluxes_n

end module longwave_fluxes
