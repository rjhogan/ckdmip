function evaluate_ckd_lw_bands(ref_file, ckd_file_in, model_name, scenario_title, specdef)
%ref_file = '/hugetmp/parr/ckdmip/evaluation1/lw_fluxes/ckdmip_evaluation1_lw_fluxes_present.h5';
%ckd_file_in = '/hugetmp/parr/ckdmip/results/ecrad-rrtmg/lw_fluxes/ecrad-rrtmg_evaluation1_lw_climate_narrow-112_fluxes_present.nc';
%specdef = loadnc(['/hugetmp/parr/ckdmip/results/ecrad-rrtmg/ecrad-rrtmg_lw_climate_narrow-112_spectral-definition.nc']);
disp(['evaluate_ckd_lw_fluxes ' ref_file ' ' ckd_file_in ' ' model_name ' ' scenario_title])
nsza = 5;

cols = {'r','b','g','m','c'};
alpha= [0.2 0.3 0.25 0.25 0.1 0.1];

titles = model_name;

if iscell(ckd_file_in)
  ckd_file = ckd_file_in;
else
  ckd_file = {ckd_file_in};
end

ckd_file


[ref0, attr_ref] = loadnc(ref_file);
ref = ref0;
if strcmp(specdef.setting(1:4),'wide')
  band_wavenumber1 = [0, 500, 820, 1180, 1800];
  band_wavenumber2 = [500, 820, 1180, 1800, 3260];   

  % Interpolate reference spectral fluxes to the wide bands
  nband = length(band_wavenumber1);
  ndim = size(ref.band_flux_dn_lw);
  ndim(1) = nband;
  ref.band_flux_dn_lw = zeros(ndim);
  ref.band_flux_up_lw = zeros(ndim);
  for iband = 1:nband
    index = find(ref0.band_wavenumber1_lw >= band_wavenumber1(iband) & ref0.band_wavenumber2_lw <= band_wavenumber2(iband));
    ref.band_flux_dn_lw(iband,:,:) = sum(ref0.band_flux_dn_lw(index,:,:),1)
    ref.band_flux_up_lw(iband,:,:) = sum(ref0.band_flux_up_lw(index,:,:),1)
  end
  ref.band_wavenumber1 = band_wavenumber1;
  ref.band_wavenumber2 = band_wavenumber2;
  paperpos = [0.5 0.5 35 20];
elseif strcmp(specdef.setting(1:4),'psla')
  band_wavenumber1 = [2500, 2200, 1900, 1700, 1400, 1250, 1100, 980,  800, 670, 540, 400, 280, 10];
  band_wavenumber2 = [2850, 2500, 2200, 1900, 1700, 1400, 1250, 1100, 980, 800, 670, 540, 400, 280];
  paperpos = [0.5 0.5 70 20];
else
  band_wavenumber1 = [0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080];
  band_wavenumber2 = [350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260];
  paperpos = [0.5 0.5 70 20];
end

flux_ref_mean{1} = mean(ref.flux_up_lw,2);
flux_ref_mean{2} = mean(ref.flux_dn_lw,2);
flux_ref_mean_band{1} = mean(ref.band_flux_up_lw,3);
flux_ref_mean_band{2} = mean(ref.band_flux_dn_lw,3);

dirstr = {'up','down'};

hr_ref = calc_hr(ref,'lw')';
hr_ref_mean = mean(hr_ref,2);
hr_ref_band = calc_spectral_hr(ref,'lw',1:50,'band');
hr_ref_mean_band = squeeze(mean(hr_ref_band,2));

nband = length(band_wavenumber1);
bandlist = 1:nband;

for ie = 1:length(ckd_file)
  ckd{ie} = loadnc(ckd_file{ie});
  hr_ckd{ie} = calc_hr(ckd{ie},'lw')';
  hr_ckd_gpt{ie} = calc_spectral_hr(ckd{ie},'lw',1:50,'spectral');
  [nhl,ncol] = size(ckd{ie}.pressure_hl);
  % Calculate fluxes and heating rates in each band
  ckd{ie}.band_flux_up_lw = calc_band_average(band_wavenumber1, band_wavenumber2,...
					      specdef, ckd{ie}.spectral_flux_up_lw);
  ckd{ie}.band_flux_dn_lw = calc_band_average(band_wavenumber1, band_wavenumber2,...
					      specdef, ckd{ie}.spectral_flux_dn_lw);
  hr_ckd_mean{ie} = mean(hr_ckd{ie},2);
  std_hr_err{ie} = 1.96.*std(hr_ckd{ie}-hr_ref,0,2);
  hr_ckd_band{ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
				      specdef, hr_ckd_gpt{ie});
  hr_ckd_mean_band{ie} = squeeze(mean(hr_ckd_band{ie},2));
  std_hr_err_band{ie}  = 1.96.*squeeze(std(hr_ckd_band{ie}-hr_ref_band,0,2));

  flux_ckd_mean{1,ie} = mean(ckd{ie}.flux_up_lw,2);
  std_flux_err{1,ie}  = 1.96.*std(ckd{ie}.flux_up_lw-ref.flux_up_lw,0,2);
  flux_ckd_band{1,ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
					  specdef, ckd{ie}.spectral_flux_up_lw);
  flux_ckd_mean_band{1,ie} = squeeze(mean(ckd{ie}.band_flux_up_lw,3));
  std_flux_err_band{1,ie} = 1.96.*squeeze(std(ckd{ie}.band_flux_up_lw-ref.band_flux_up_lw,0,3));

  flux_ckd_mean{2,ie} = mean(ckd{ie}.flux_dn_lw,2);
  std_flux_err{2,ie}  = 1.96.*std(ckd{ie}.flux_dn_lw-ref.flux_dn_lw,0,2);
  flux_ckd_band{2,ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
					  specdef, ckd{ie}.spectral_flux_dn_lw);
  flux_ckd_mean_band{2,ie} = squeeze(mean(ckd{ie}.band_flux_dn_lw,3));
  std_flux_err_band{2,ie} = 1.96.*squeeze(std(ckd{ie}.band_flux_dn_lw-ref.band_flux_dn_lw,0,3));
end

pmid = 0.5.*(ref.pressure_hl(1:end-1,1) + ref.pressure_hl(2:end,1))./100;
phl = ref.pressure_hl(:,1)./100;
nx = nband+1;

%ssi_lbl = lbl.band_flux_dn_lw(:,1,1);
ref_lw = 1.5;
%ref_col = [0.5 0.5 0.5];


figure(4)
clf
set(gcf,'paperposition',paperpos);


for idir = 1:2


  %%%%%%%%%%%%%%%%%%% BAND FLUX %%%%%%%%%%%%%%%%%%%
  for iband = bandlist
    subplot(3,nx,nx*(2-idir)+1+iband)
    for ie = 1:length(ckd_file)
      h=fill([flux_ckd_mean_band{idir,ie}(iband,:)+std_flux_err_band{idir,ie}(iband,:) ...
		    flip(flux_ckd_mean_band{idir,ie}(iband,:)-std_flux_err_band{idir,ie}(iband,:))],...
	     [phl' flip(phl')],'k');
      set(h,'facecolor',cols{ie},'edgecolor','none','facealpha',alpha(ie));
      hold on
      semilogy(flux_ckd_mean_band{idir,ie}(iband,:), phl, cols{ie},'linewidth',ref_lw);
      ylim([0.01 1013]);
    end
    semilogy(flux_ref_mean_band{idir}(iband,:), phl, 'k--','linewidth',ref_lw);
    set(gca,'yscale','log','ydir','reverse','ytick',10.^[-2:3],'layer','top');
    xx = xlim;
    xlabel(['Flux ' dirstr{idir} ' (W m^{-2})']);
    if idir == 2
      title([num2str(band_wavenumber1(iband)) '-' num2str(band_wavenumber2(iband)) ' cm^{-1}']);
    end

  end
  %%%%%%%%%%%%%%%%%%% BROADBAND FLUX %%%%%%%%%%%%%%%%%%%
  subplot(3,nx,1+nx.*(2-idir))
  if idir == 2
    plot(0,-1,'k--','linewidth',1.5);
    hold on
    plot(0,-1,cols{1},'linewidth',1.5);
  end

  for ie = 1:length(ckd_file)
    h=fill([flux_ckd_mean{idir,ie}'+std_flux_err{idir,ie}' ...
	    flip(flux_ckd_mean{idir,ie}'-std_flux_err{idir,ie}')],...
	   [phl' flip(phl')],'k');
    set(h,'facecolor',cols{ie},'edgecolor','none','facealpha',alpha(ie));
    hold on
    semilogy(flux_ckd_mean{idir,ie}, phl, cols{ie},'linewidth',ref_lw);
  end
  semilogy(flux_ref_mean{idir}, phl,'k--','linewidth',ref_lw);
  ylim([0.01 1013]);
  set(gca,'yscale','log','layer','top');
  set(gca,'ydir','reverse','ytick',10.^[-2:3]);
  if idir == 2
    title('Longwave')
    legend({'Reference',titles},'location','northwest');
  end
  ylabel('Pressure (hPa)');
  xlabel(['Flux ' dirstr{idir} ' (W m^{-2})']);
end

%%%%%%%%%%%%%%%%%%% BROADBAND HEATING RATE %%%%%%%%%%%%%%%%%%%
subplot(3,nx,nx*2+1)
for ie = 1:length(ckd_file)
  h=fill([hr_ckd_mean{ie}'+std_hr_err{ie}' ...
	  flip(hr_ckd_mean{ie}'-std_hr_err{ie}')],...
	 [pmid' flip(pmid')],'k');
  set(h,'facecolor',cols{ie},'edgecolor','none','facealpha',alpha(ie));
  hold on
  semilogy(hr_ckd_mean{ie}, pmid, cols{ie},'linewidth',ref_lw);
end
set(gca,'yscale','log','layer','top');
semilogy(hr_ref_mean, pmid,'k--','linewidth',ref_lw);
ylim([0.01 1013]);
set(gca,'ydir','reverse','ytick',10.^[-2:3]);
ylabel('Pressure (hPa)');
xlabel('Heating rate (K d^{-1})');
plot([0 0],[0.01 1013],'k:');

%%%%%%%%%%%%%%%%%%% BAND HEATING RATE %%%%%%%%%%%%%%%%%%%
for iband = bandlist
  subplot(3,nx,nx*2+1+iband)
  for ie = 1:length(ckd_file)
    h=fill([hr_ckd_mean_band{ie}(iband,:)+std_hr_err_band{ie}(iband,:) ...
	    flip(hr_ckd_mean_band{ie}(iband,:)-std_hr_err_band{ie}(iband,:))],...
	   [pmid' flip(pmid')],'k');
    set(h,'facecolor',cols{ie},'edgecolor','none','facealpha',alpha(ie));
    hold on
    semilogy(hr_ckd_mean_band{ie}(iband,:), pmid, cols{ie},'linewidth',ref_lw);
  end
  semilogy(hr_ref_mean_band(iband,:), pmid, 'k--','linewidth',ref_lw);

  ylim([0.01 1013]);
  set(gca,'yscale','log','layer','top');
  set(gca,'ydir','reverse','ytick',10.^[-2:3]);
  %  ylabel('Pressure (hPa)');
  xlabel('Heating rate (K d^{-1})');
  %set(gca,'yticklabel','');
  %title([num2str(band_wavenumber1(iband)) '-' num2str(band_wavenumber2(iband)) ' cm^{-1}']);
  plot([0 0],[0.01 1013],'k:');
end
drawnow
