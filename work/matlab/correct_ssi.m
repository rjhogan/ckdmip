function ckd_out = correct_ssi(ref, ckd_in, specdef)
  nband = length(ref.band_wavenumber1_sw);

  iband_g = zeros(nband,size(specdef.gpoint_fraction,2));

  for iband = 1:nband
    index = find(specdef.wavenumber1 >= ref.band_wavenumber1_sw(iband) ...
		 & specdef.wavenumber2 <= ref.band_wavenumber2_sw(iband));
    iband_g(iband,:) = sum(specdef.gpoint_fraction(index,:),1);
  end

  ckd_out.iband_g = iband_g;
  ssi_g_orig = squeeze(ckd_in.spectral_flux_dn_sw(:,1,end,1)) ./ ckd_in.mu0(end);
  ssi_band_orig = iband_g * ssi_g_orig;
  ssi_band_ref  = ref.band_flux_dn_sw(:,1,end,1) ./ ref.mu0(end);

  scaling_band = ssi_band_ref ./ ssi_band_orig

  for iband = 1:nband
    disp(sprintf('Band %05d-%05d: scaling = %f',ref.band_wavenumber1_sw(iband), ...
		 ref.band_wavenumber2_sw(iband), scaling_band(iband)));
  end

  scaling_g    = iband_g' * scaling_band; % Does this work if not 0s and 1s?
  ckd_out = ckd_in;
  for ig = 1:size(ckd_out.spectral_flux_dn_sw,1)
    ckd_out.spectral_flux_dn_sw(ig,:,:,:) = ckd_out.spectral_flux_dn_sw(ig,:,:,:).*scaling_g(ig);
    ckd_out.spectral_flux_up_sw(ig,:,:,:) = ckd_out.spectral_flux_up_sw(ig,:,:,:).*scaling_g(ig);
    ckd_out.spectral_flux_dn_direct_sw(ig,:,:,:) = ckd_out.spectral_flux_dn_direct_sw(ig,:,:,:).*scaling_g(ig);
  end
  ckd_out.flux_dn_sw = squeeze(sum(ckd_out.spectral_flux_dn_sw,1));
  ckd_out.flux_up_sw = squeeze(sum(ckd_out.spectral_flux_up_sw,1));
  ckd_out.flux_dn_direct_sw = squeeze(sum(ckd_out.spectral_flux_dn_direct_sw,1));


