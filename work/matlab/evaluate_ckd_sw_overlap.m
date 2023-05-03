
gases={'co2','ch4','n2o'};
gas_names = {'CO_2','CH_4','N_2O'};
units={'ppmv','ppbv','ppbv'};
concs={[180 415 2240],...
       [350 1921 3500],...
       [190 332 540]};
conc_axes = {[100 3000],[0 4000],[100 600]};
rf_axis = {[-0.8 0.8],[-0.4 0.4],[-0.1 0.1]};
concs_present = [415 1921 332];
islog = [1 0 0];
my_xlim = [0.4 3.6];
barcols = {'k','b','r'};


ngas = length(gases)
ref_dir = '/hugetmp/parr/ckdmip/evaluation1/sw_fluxes';
ckd_prefixes = {'/hugetmp/parr/ckdmip/results/ecrad-rrtmg/sw_fluxes/ecrad-rrtmg_evaluation1_sw_climate_narrow-112_fluxes_',...
		'/hugetmp/parr/ckdmip/results/RTE-RRTMGP-181204/sw_fluxes/RTE-RRTMGP-181204_evaluation1_sw_climate_narrow-256_fluxes_'};
titles = {'ecRad-RRTMG','RTE-RRTMGP'};
%ckd_prefixes = ckd_prefixes(2);
%titles = titles(1);

ickd = 1:length(ckd_prefixes);
ickd = 1;

col_ref = 'k';
cols = {'b-.','r--','g','m','c'};
ngas = length(gases);
combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end

gas1 = [1 1 2];
gas2 = [2 3 3];

ncase = length(gas1)

figure(1)
clf
set(gcf,'defaultlinelinewidth',1);
set(gcf,'paperposition',[0.5 0.5 27 15]);

for igas = 1:ncase

  igas1 = gas1(igas);
  igas2 = gas2(igas);

  clear flux_* ipresent

  for iconc1 = 1:3
    str1 = '';
    if iconc1 ~= 2
      str1 = [gases{igas1} '-' num2str(concs{igas1}(iconc1))];
    end
    for iconc2 = 1:3
      str2 = '';
      if iconc2 ~= 2
	str2 = [gases{igas2} '-' num2str(concs{igas2}(iconc2))];
      end
      if isempty(str1)
	if isempty(str2)
	  scenario = 'present';
	else
	  scenario = str2;
	end
      else
	if isempty(str2)
	  scenario = str1
	else
	  scenario = [str1 '-' str2];
	end
      end
      ref_orig = loadnc([ref_dir '/ckdmip_evaluation1_sw_fluxes_' scenario '.h5']);
      ref = flatten_sza(ref_orig);
      flux_toa_ref(iconc1,iconc2)  = -mean(ref.flux_up_sw(1,:));
      flux_surf_ref(iconc1,iconc2) = mean(ref.flux_dn_sw(end,:));
      for ie = ickd
	ckd_orig = loadnc([ckd_prefixes{ie} scenario '.nc']);
	ckd = flatten_sza(ckd_orig);
	flux_toa_ckd{ie}(iconc1,iconc2)  = -mean(ckd.flux_up_sw(1,:));
	flux_surf_ckd{ie}(iconc1,iconc2) = mean(ckd.flux_dn_sw(end,:));
      end
    end
  end

  if 0
  subplot(2,ngas,igas)
  plot(concs{igas2}, flux_toa_ref'-flux_toa_ref(2,2), [col_ref 'o-']);
  hold on
  for ie = ickd
    plot(concs{igas2}, flux_toa_ckd{ie}'-flux_toa_ckd{ie}(2,2), cols{ie});
  end
  if islog(igas2)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas2})
  xlabel([gas_names{igas2} ' volume mixing ratio (' units{igas2} ')']);
  ylabel('TOA radiative forcing w.r.t. present (W m^{-2})');
  if igas2 == 1
    legend(['Reference LBL' titles],'location','northwest');
  end

  end

  if 0

  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  title(['\bf' gas_names{igas}]);
  axis(ax);

  subplot(2,ngas,igas+ngas)
  plot(concs{igas}, flux_surf_ref-flux_surf_ref(ipresent), [col_ref 'o-']);
  hold on
  for ie = ickd
    plot(concs{igas}, flux_surf_ckd{ie}-flux_surf_ckd{ie}(ipresent), cols{ie});
  end
  if islog(igas)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas})
  xlabel([gas_names{igas} ' volume mixing ratio (' units{igas} ')']);
  ylabel('Surface radiative forcing w.r.t. present (W m^{-2})');
  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  axis(ax);

  else

  ibaseconc = 2; % present
  %ibaseconc = 1; % glacial

  subplot(2,ngas,igas)
  bars(1:3,1) = flux_toa_ref(3,:)-flux_toa_ref(ibaseconc,:);
  bars_dn(1:3,1) = flux_toa_ref(1,:)-flux_toa_ref(ibaseconc,:);

  overlap_effect_ref(igas1,igas2) = bars(3,1)-bars(1,1);
  overlap_effect_pc_ref(igas1,igas2) = 100.*(bars(3,1)-bars(1,1))./bars(2,1);
  overlap_effect_dn_ref(igas1,igas2) = bars_dn(3,1)-bars_dn(1,1);
  overlap_effect_dn_pc_ref(igas1,igas2) = 100.*(bars_dn(3,1)-bars_dn(1,1))./bars_dn(2,1);

  for ie = ickd
    bars(1:3,ie+1) = flux_toa_ckd{ie}(3,:)-flux_toa_ckd{ie}(ibaseconc,:);
    bars_dn(1:3,ie+1) = flux_toa_ckd{ie}(1,:)-flux_toa_ckd{ie}(ibaseconc,:);
  end

  %bar(concs{igas2},bars)
  h=bar(1:3,bars);
  for ii = 1:length(h)
    h(ii).FaceColor = barcols{ii};
  end
  set(gca,'xticklabel',concs{igas2})
  hold on
  h=bar(1:3,bars_dn);
  for ii = 1:length(h)
    h(ii).FaceColor = barcols{ii};
  end
  xlabel([gas_names{igas2} ' concentration (' units{igas2} ')']);
  ylabel([gas_names{igas1} ' radiative forcing (W m^{-1})']);
  ylim(rf_axis{igas1})
  xlim(my_xlim);

  if igas == 1
    legend('LBLRTM','RRTMG','RRTMGP','location','east');
  end 

  subplot(2,ngas,igas+ngas)
  bars(1:3,1) = (flux_toa_ref(:,3)-flux_toa_ref(:,ibaseconc))';
  bars_dn(1:3,1) = (flux_toa_ref(:,1)-flux_toa_ref(:,ibaseconc))';

  overlap_effect_ref(igas2,igas1) = bars(3,1)-bars(1,1);
  overlap_effect_pc_ref(igas2,igas1) = 100.*(bars(3,1)-bars(1,1))./bars(2,1);
  overlap_effect_dn_ref(igas2,igas1) = bars_dn(3,1)-bars_dn(1,1);
  overlap_effect_dn_pc_ref(igas2,igas1) = 100.*(bars_dn(3,1)-bars_dn(1,1))./bars_dn(2,1);

  for ie = ickd
    bars(1:3,ie+1) = (flux_toa_ckd{ie}(:,3)-flux_toa_ckd{ie}(:,ibaseconc))';
    bars_dn(1:3,ie+1) = (flux_toa_ckd{ie}(:,1)-flux_toa_ckd{ie}(:,ibaseconc))';
  end

  h=bar(1:3,bars);
  for ii = 1:length(h)
    h(ii).FaceColor = barcols{ii};
  end
  set(gca,'xticklabel',concs{igas1})
  hold on
  h=bar(1:3,bars_dn);
  for ii = 1:length(h)
    h(ii).FaceColor = barcols{ii};
  end
  xlabel([gas_names{igas1} ' concentration (' units{igas1} ')']);
  ylabel([gas_names{igas2} ' radiative forcing (W m^{-1})']);
  ylim(rf_axis{igas2})
  xlim(my_xlim);

  end

end

%print_png([combined_titles 'evaluation1_forcing.png'],100)
