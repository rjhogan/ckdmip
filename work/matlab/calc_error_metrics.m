function [metric, metric_names] = calc_error_metrics(ref_files, ckd_files)

  ref.pressure_hl = [];
  ref.flux_up_lw = [];
  ref.flux_dn_lw = [];
  ckd = ref;

  for ifile = 1:length(ref_files)
    ref0 = loadnc(ref_files{ifile});
    ckd0 = loadnc(ckd_files{ifile});
    ref.pressure_hl = [ref.pressure_hl ref0.pressure_hl];
    ref.flux_up_lw  = [ref.flux_up_lw  ref0.flux_up_lw];
    ref.flux_dn_lw  = [ref.flux_dn_lw  ref0.flux_dn_lw];
    ckd.pressure_hl = [ckd.pressure_hl ckd0.pressure_hl];
    ckd.flux_up_lw  = [ckd.flux_up_lw  ckd0.flux_up_lw];
    ckd.flux_dn_lw  = [ckd.flux_dn_lw  ckd0.flux_dn_lw];
  end

  p_fl = 0.5.*(ref.pressure_hl(2:end,:)+ref.pressure_hl(1:end-1,:));
  %dp   = ref.pressure_hl(2:end,:)-ref.pressure_hl(1:end-1,:);

  hr_ref = calc_hr(ref,'lw')';
  hr_ckd = calc_hr(ckd,'lw')';

  toa_up_err = ckd.flux_up_lw(1,:)-ref.flux_up_lw(1,:);
  metric(1) = mean(toa_up_err);
  metric(2) = std(toa_up_err);
  surf_dn_err = ckd.flux_dn_lw(end,:)-ref.flux_dn_lw(end,:);
  metric(3) = mean(surf_dn_err);
  metric(4) = std(surf_dn_err);
  
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
