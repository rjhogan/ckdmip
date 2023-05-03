function [hr, spectral_hr] = calc_hr(data, band, icol)

eval(['flux_up = data.flux_up_' band ';']);
eval(['flux_dn = data.flux_dn_' band ';']);
flux_net = flux_dn-flux_up;
g=9.81;
scaling = 3600.*24;
hr_all = -scaling.*(diff(flux_net).*g./diff(data.pressure_hl)./1004)';
if nargin > 2
  hr = hr_all(icol,:);
else
  hr = hr_all;
end

if isfield(data, ['spectral_flux_up_' band])
  spectral_flux_net = data.(['spectral_flux_dn_' band]) ...
		      - data.(['spectral_flux_up_' band]);
  for ig = 1:size(spectral_flux_net,1)
    spectral_hr_all(ig,:,:) = -scaling.*(diff(squeeze(spectral_flux_net(ig,:,:))).*g./diff(data.pressure_hl)./1004)';
  end
  if nargin > 2
    spectral_hr = squeeze(spectral_hr_all(:,icol,:));
  else
    spectral_hr = spectral_hr_all;
  end
end
