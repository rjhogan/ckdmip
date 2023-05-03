function evaluate_ckd_sw_bands(ref_file, ckd_file_in, model_name, scenario_title, specdef, fix_ssi)
%ref_file = '/hugetmp/parr/ckdmip/evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_present.h5';
%ckd_file_in = '/hugetmp/parr/ckdmip/results/ecrad-rrtmg/sw_fluxes/ecrad-rrtmg_evaluation1_sw_climate_narrow-112_fluxes_present.nc';
%specdef = loadnc(['/hugetmp/parr/ckdmip/results/ecrad-rrtmg/ecrad-rrtmg_sw_climate_narrow-112_spectral-definition.nc']);
disp(['evaluate_ckd_sw_fluxes ' ref_file ' ' ckd_file_in ' ' model_name ' ' scenario_title])
nsza = 5;isza_list = 1:nsza;
%nsza = 1;isza_list = 5;
%nsza = 5;isza_list = [3 3 3 3 3];

cols = {'r','b','g','m','c'};
alpha= [0.2 0.3 0.25 0.25 0.1 0.1];

titles = model_name;
if fix_ssi == 1
  if ~iscell(ckd_file_in)
    %ckd_file_in = {ckd_file_in, ckd_file_in};
    do_fix_ssi = [1];
    titles = [model_name ' (fix SSI)'];
    %titles = {model_name, [model_name ' (fix SSI)']};
  else
    do_fix_ssi = ones(size(ckd_file_in));
  end
else
  do_fix_ssi = zeros(size(ckd_file_in));
  if fix_ssi == -1
    cols{1} = 'b';
    cols{2} = 'r';
  end
end

if iscell(ckd_file_in)
  ckd_file = ckd_file_in;
else
  ckd_file = {ckd_file_in};
end

ckd_file

dirstr = {'up','down'};

vars = {'pressure_hl','temperature_hl','flux_up_sw','flux_dn_sw','flux_dn_direct_sw',...
	'band_wavenumber1_sw','band_wavenumber2_sw','band_flux_up_sw','band_flux_dn_sw',...
	'band_flux_dn_direct_sw','mu0'};

[ref0, attr_ref] = loadnc(ref_file, vars);

% ecckd 0.6
%ref = flatten_sza(ref0,isza_list);

% ecckd 0.7, needed for regridding
ref = flatten_sza(ref0,isza_list);
ref0 = ref;

paperpos = [0.5 0.5 35 20];
if strcmp(specdef.setting(1:4),'wide')
  band_wavenumber1 = [250 4000, 8050, 16000 29000];
  band_wavenumber2 = [4000 8050 16000 29000 50000];  
  do_reband = 1;
elseif strcmp(specdef.setting(1:4),'doub')
  %band_wavenumber1 = [250 12850];
  %band_wavenumber2 = [12850 50000];  
  band_wavenumber1 = [250 16000];
  band_wavenumber2 = [16000 50000];  
  do_reband = 1;
elseif strcmp(specdef.setting(1:4),'psla')
  % PSLACKD band definitions
  band_wavenumber1 = [44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000, 800];
  band_wavenumber2 = [57000, 44500, 41000, 35000, 33500, 31008, 27972, 22857, 20101, 16807, 14500, 12600, 11250, 9600, 7090, 5250, 4000, 2850, 2500, 2000];
  paperpos = [0.5 0.5 100 20];
  do_reband = 0;
  ref_lw = 1;
elseif strcmp(specdef.setting(1:3),'rgb')
  %band_wavenumber1 = [250 14286 16667 20000 25000];
  %band_wavenumber2 = [14286 16667 20000 25000 50000];
  band_wavenumber1 = [250 14300 16650 20000 25000];
  band_wavenumber2 = [14300 16650 20000 25000 50000];
  do_reband = 1;
elseif strcmp(specdef.setting(1:2),'gb')
  band_wavenumber1 = [250 8000 16650 20000 25000];
  band_wavenumber2 = [8000 16650 20000 25000 50000];
  do_reband = 1;
elseif strcmp(model_name,'MSTRN')
  band_wavenumber1 = [250, 2600, 3250, 6150, 8050, 12850];
  band_wavenumber2 = [2600, 3250, 6150, 8050, 12850, 50000];
  do_reband = 1;
  paperpos = [0.5 0.5 45 20];
elseif strcmp(specdef.setting(1:4),'narr')
  band_wavenumber1 = [250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000];
  band_wavenumber2 = [2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000];
  paperpos = [0.5 0.5 70 20];
  do_reband = 0;
else
  band_wavenumber1 = specdef.wavenumber1_band;
  band_wavenumber2 = specdef.wavenumber2_band;
  if band_wavenumber2(end) >= 49999
     band_wavenumber2(end) = 50000
  end
  paperpos = [0.5 0.5 (length(specdef.wavenumber1_band)+2).*4.55 20];
  do_reband = 1;
end

if (do_reband)
  % Interpolate reference spectral fluxes to the wide bands
  nband = length(band_wavenumber1);
  ndim = size(ref.band_flux_dn_sw);
  ndim(1) = nband;
  ref.band_flux_dn_sw = zeros(ndim);
  ref.band_flux_up_sw = zeros(ndim);
  ref0.band_wavenumber1_sw
  ref0.band_wavenumber2_sw
  band_wavenumber1
  band_wavenumber2
  for iband = 1:nband
    index = find(ref0.band_wavenumber1_sw >= band_wavenumber1(iband) & ref0.band_wavenumber2_sw <= band_wavenumber2(iband));
    ref.band_flux_dn_sw(iband,:,:) = sum(ref0.band_flux_dn_sw(index,:,:),1);
    ref.band_flux_up_sw(iband,:,:) = sum(ref0.band_flux_up_sw(index,:,:),1);
  end
  ref.band_wavenumber1 = band_wavenumber1;
  ref.band_wavenumber2 = band_wavenumber2;
end

ref_orig = ref0;

flux_ref_mean{1} = mean(ref.flux_up_sw,2);
flux_ref_mean{2} = mean(ref.flux_dn_sw,2);
flux_ref_mean_band{1} = mean(ref.band_flux_up_sw,3);
flux_ref_mean_band{2} = mean(ref.band_flux_dn_sw,3);

nband = length(band_wavenumber1);
bandlist = 1:nband;

hr_ref = calc_hr(ref,'sw')';
hr_ref_mean = mean(hr_ref,2);
%hr_ref_band = calc_spectral_hr(ref,'sw',1:250,'band');
hr_ref_band = calc_spectral_hr(ref,'sw',1:size(ref.band_flux_up_sw,3),'band');
hr_ref_mean_band = squeeze(mean(hr_ref_band,2));

for ie = 1:length(ckd_file)
  ckd_tmp = loadnc(ckd_file{ie});
  if do_fix_ssi(ie)
    ckd_tmp = correct_ssi(ref_orig, ckd_tmp, specdef);
  end
  ckd{ie} = flatten_sza(ckd_tmp,isza_list);
  hr_ckd{ie} = calc_hr(ckd{ie},'sw')';
  hr_ckd_gpt{ie} = calc_spectral_hr(ckd{ie},'sw',1:size(ckd{ie}.spectral_flux_up_sw,3),'spectral');
  [nhl,ncol] = size(ckd{ie}.pressure_hl);
  % Calculate fluxes and heating rates in each band
  ckd{ie}.band_flux_up_sw = calc_band_average(band_wavenumber1, band_wavenumber2,...
					      specdef, ckd{ie}.spectral_flux_up_sw);
  ckd{ie}.band_flux_dn_sw = calc_band_average(band_wavenumber1, band_wavenumber2,...
					      specdef, ckd{ie}.spectral_flux_dn_sw);
  ckd{ie}.band_flux_dn_direct_sw = calc_band_average(band_wavenumber1, band_wavenumber2,...
						     specdef, ckd{ie}.spectral_flux_dn_direct_sw);
  hr_ckd_mean{ie} = mean(hr_ckd{ie},2);
  std_hr_err{ie} = 1.96.*std(hr_ckd{ie}-hr_ref,0,2);
  hr_ckd_band{ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
				      specdef, hr_ckd_gpt{ie});
  hr_ckd_mean_band{ie} = squeeze(mean(hr_ckd_band{ie},2));
  std_hr_err_band{ie}  = 1.96.*squeeze(std(hr_ckd_band{ie}-hr_ref_band,0,2));

  flux_ckd_mean{1,ie} = mean(ckd{ie}.flux_up_sw,2);
  std_flux_err{1,ie}  = 1.96.*std(ckd{ie}.flux_up_sw-ref.flux_up_sw,0,2);
  flux_ckd_band{1,ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
					  specdef, ckd{ie}.spectral_flux_up_sw);
  flux_ckd_mean_band{1,ie} = squeeze(mean(ckd{ie}.band_flux_up_sw,3));
  std_flux_err_band{1,ie} = 1.96.*squeeze(std(ckd{ie}.band_flux_up_sw-ref.band_flux_up_sw,0,3));

  flux_ckd_mean{2,ie} = mean(ckd{ie}.flux_dn_sw,2);
  std_flux_err{2,ie}  = 1.96.*std(ckd{ie}.flux_dn_sw-ref.flux_dn_sw,0,2);
  flux_ckd_band{2,ie} = calc_band_average(band_wavenumber1, band_wavenumber2,...
					  specdef, ckd{ie}.spectral_flux_dn_sw);
  flux_ckd_mean_band{2,ie} = squeeze(mean(ckd{ie}.band_flux_dn_sw,3));
  std_flux_err_band{2,ie} = 1.96.*squeeze(std(ckd{ie}.band_flux_dn_sw-ref.band_flux_dn_sw,0,3));
end

pmid = 0.5.*(ref.pressure_hl(1:end-1,1) + ref.pressure_hl(2:end,1))./100;
phl = ref.pressure_hl(:,1)./100;
nx = nband+1;

%ssi_lbl = lbl.band_flux_dn_sw(:,1,1);
ref_lw = 1.5;
%ref_col = [0.5 0.5 0.5];


figure(4)
clf
set(gcf,'paperposition', paperpos);


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
    if xx(1) < 0
      xlim(xx.*[0 1]);
    end
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
    title('Shortwave')
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
%xlim(xlim.*[0 1]);
xlim([0 20]);
ylim([0.01 1013]);
set(gca,'ydir','reverse','ytick',10.^[-2:3]);
ylabel('Pressure (hPa)');
xlabel('Heating rate (K d^{-1})');

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

  xlim(xlim.*[0 1]);
  ylim([0.01 1013]);
  set(gca,'yscale','log','layer','top');
  set(gca,'ydir','reverse','ytick',10.^[-2:3]);
  %  ylabel('Pressure (hPa)');
  xlabel('Heating rate (K d^{-1})');
  set(gca,'yticklabel','');
  %title([num2str(band_wavenumber1(iband)) '-' num2str(band_wavenumber2(iband)) ' cm^{-1}']);
  if iband == 1
    xlim([0 5]);
  end
end
drawnow
