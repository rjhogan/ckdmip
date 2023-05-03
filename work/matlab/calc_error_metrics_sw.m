function [metric, metric_names] = calc_error_metrics(ref_files, ckd_files)

  ref.pressure_hl = [];
  ref.flux_up_sw = [];
  ref.flux_dn_sw = [];
  ckd = ref;

  for ifile = 1:length(ref_files)
    ref0 = flatten_sza(loadnc(ref_files{ifile}));
    ckd0 = flatten_sza(loadnc(ckd_files{ifile}));
    ref.pressure_hl = [ref.pressure_hl ref0.pressure_hl];
    ref.flux_up_sw  = [ref.flux_up_sw  ref0.flux_up_sw];
    ref.flux_dn_sw  = [ref.flux_dn_sw  ref0.flux_dn_sw];
    ckd.pressure_hl = [ckd.pressure_hl ckd0.pressure_hl];
    ckd.flux_up_sw  = [ckd.flux_up_sw  ckd0.flux_up_sw];
    ckd.flux_dn_sw  = [ckd.flux_dn_sw  ckd0.flux_dn_sw];
  end

  p_fl = 0.5.*(ref.pressure_hl(2:end,:)+ref.pressure_hl(1:end-1,:));
  %dp   = ref.pressure_hl(2:end,:)-ref.pressure_hl(1:end-1,:);

  hr_ref = calc_hr(ref,'sw')';
  hr_ckd = calc_hr(ckd,'sw')';

  toa_up_err = ckd.flux_up_sw(1,:)-ref.flux_up_sw(1,:);
  metric(1) = mean(toa_up_err);
  %metric(2) = std(toa_up_err);
  metric(2) = sqrt(mean(toa_up_err.^2));
  surf_dn_err = ckd.flux_dn_sw(end,:)-ref.flux_dn_sw(end,:);
  metric(3) = mean(surf_dn_err);
  %metric(4) = std(surf_dn_err);
  metric(4) = sqrt(mean(surf_dn_err.^2));

  weight = ref.pressure_hl(2:end,:).^(1./3)-ref.pressure_hl(1:end-1,:).^(1./3);
  index = find(p_fl >= 400);
  metric(5) = sqrt(sum((hr_ckd(index)-hr_ref(index)).^2.*weight(index))./sum(weight(index)));

  index = find(p_fl < 400 & p_fl >= 2);
  metric(6) = sqrt(sum((hr_ckd(index)-hr_ref(index)).^2.*weight(index))./sum(weight(index)));

  metric_names = {'Bias in TOA upwelling (W m^{-2})',...
		  'RMSE in TOA upwelling (W m^{-2})',...
		  'Bias in surface downwelling (W m^{-2})',...
		  'RMSE in surface downwelling (W m^{-2})',...
		  'RMSE in 4-1100 hPa heating rate (K d^{-1})',...
		  'RMSE in 0.02-4 hPa heating rate (K d^{-1})'};
