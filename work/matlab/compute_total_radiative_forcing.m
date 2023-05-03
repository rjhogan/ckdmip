% This script uses CKDMIP's line-by-line radiative transfer
% calculations in which the concentration of individual gases have
% been reduced to zero, to compute the present-day radiative forcing of
% those gases compared to zero concentration.

directory='/hugetmp/parr/ckdmip/evaluation1';

domains = {'sw','lw'};
flstr = {'fluxes','fluxes-4angle'};
%gases = {'co2','ch4','n2o','o2','n2','cfc11'}%,'cfc12'};
gases = {'co2','ch4','n2o','cfc11','cfc12','n2','o2'};
lowval = {'0','0','0','0','0','0','0'};
lowval = {'280','700','270','0','0','0','0'};

%gases = {'co2'}

in = loadnc(['/hugetmp/parr/ckdmip/evaluation1/conc/ckdmip_evaluation1_concentrations_present.nc']);

for id = 1:length(domains)
  now = loadnc([directory '/' domains{id} '_fluxes/ckdmip_evaluation1_' domains{id} '_' flstr{id} '_present.h5']);

  prefix = [directory '/' domains{id} '_fluxes/ckdmip_evaluation1_' domains{id} '_' flstr{id} '_'];
  for igas = 1:length(gases)
    z = loadnc([prefix gases{igas} '-' lowval{igas} '.h5']);
    if id == 1
      toa_r(:,igas,id) = squeeze(mean(z.flux_up_sw(1,:,:)-now.flux_up_sw(1,:,:),2));
      surf_r(:,igas,id) = squeeze(mean(now.flux_dn_sw(end,:,:)-z.flux_dn_sw(end,:,:) ...
		     -now.flux_up_sw(end,:,:)+z.flux_up_sw(end,:,:),2));
      toa_f(igas,id) = mean(toa_r(:,igas,id));
      surf_f(igas,id) = mean(surf_r(:,igas,id));
      toa_std(igas,id) = std(toa_r(:,igas,id));
      surf_std(igas,id) = std(surf_r(:,igas,id));
    else
      toa_r(:,igas,id) = z.flux_up_lw(1,:)-now.flux_up_lw(1,:);
      surf_r(:,igas,id) = now.flux_dn_lw(end,:)-z.flux_dn_lw(end,:) ...
			  -now.flux_up_lw(end,:)+z.flux_up_lw(end,:);
      toa_f(igas,id) = mean(toa_r(:,igas,id));
      surf_f(igas,id) = mean(surf_r(:,igas,id));
      toa_std(igas,id) = std(toa_r(:,igas,id));
      surf_std(igas,id) = std(surf_r(:,igas,id));
    end
  end
end

for id = 1:length(domains)
  disp(['Domain: ' domains{id}])    
  for igas = 1:length(gases)
    disp(['  Gas: ' gases{igas}]);
    disp(['    TOA:     ' num2str(toa_f(igas,id)) ' +/- ' num2str(toa_std(igas,id))]);
    disp(['    surface: ' num2str(surf_f(igas,id)) ' +/- ' num2str(surf_std(igas,id))]);
  end
end
