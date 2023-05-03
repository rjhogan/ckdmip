function evaluate_ckd_sw_forcing(ref_base, ckd_base_in, model_name)

if iscell(ckd_base_in)
  ckd_base = ckd_base_in;
else
  ckd_base = {ckd_base_in};
end

if iscell(model_name)
  titles = model_name;
else
  titles = {model_name}
end

gases={'co2','ch4','n2o'};%,'cfc11','cfc12'};
%gases={'co2','ch4'};%,'n2o'};%,'cfc11','cfc12'};
hrxlim={[-2 5],[-0.02 0.02],[-0.02 0.02]};
hrxlim={[-2 5],[-0.03 0.03],[-0.02 0.02]};
gas_names = {'CO_2','CH_4','N_2O','CFC11','CFC12'};
units={'ppmv','ppbv','ppbv','pptv','pptv'};
concs={[180 280 415 560 1120 2240],...
       [350 700 1200 1921 2600 3500],...
       [190 270 332 405 540],...
       [0 861 2000],...
       [0 495 550]};
%concs{1} = [180 280 415 1120 2240];
conc_axes = {[100 3000],[0 4000],[100 600],[0 2000],[0 600]};
concs_present = [415 1921 332 861 495];
islog = [1 0 0 0 0];

isza = [1:5]; % All solar zenith angles
%isza = [3];   % Just 60 degrees

%etminan_conc{1} = [180 389 700 2000];
%etminan_forcing{1} = -[0.16 0 -0.14 -0.43];
%etminan_conc{2} = [340 750 1800 3500];
%etminan_forcing{2} = [-0.05 -0.03 0 0.03];
%etminan_conc{3} = [200 323 535];
%etminan_forcing{3} = [0 0 -0.1];

% Only have daytime solar zenith angles, so scale to include night
% time
dayfactor = 1;

plot_hr = 1;

ngas = length(gases)
do_fix_ssi = [0 1];
do_fix_ssi = [0 0];
do_flatten = [1 1];
ickd = 1:length(ckd_base);
%ickd = 1;
col_ref = 'k';
col_future = 'm';
col_glacial = 'g';
cols = {'r--','b-.','g','m','c'};
if length(ickd) > 1
  cols = {'b-.','r--','g','m','c'};
end
ngas = length(gases);
combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end

%specdef = loadnc('/hugetmp/parr/ckdmip/results/ecrad-rrtmg/sw_spectral-definition/ecrad-rrtmg_sw_climate_narrow-112_spectral-definition.nc');

clf
set(gcf,'defaultlinelinewidth',1);
%set(gcf,'paperposition',[0.5 0.5 ngas*10+1 (2+plot_hr).*8+1]);
%set(gcf,'paperposition',[0.5 0.5 ngas*10+1 20+plot_hr.*6]);
%set(gcf,'paperposition',[0.5 0.5 ngas*9+1 18+plot_hr.*5]);
set(gcf,'paperposition',[0.5 0.5 ngas*7+1 14+plot_hr.*4]);
for igas = 1:ngas

  clear flux_* ipresent
  clear ckd ref

  for iconc = 1:length(concs{igas})
    conc = concs{igas}(iconc);
    scenario = [gases{igas} '-' num2str(conc)];
    if conc == concs_present(igas)
      scenario = 'present';
      ipresent = iconc;
    end
    ref_orig = loadnc([ref_base scenario '.h5']);
    ref{iconc} = flatten_sza(ref_orig, isza);
    flux_toa_ref(iconc)  = dayfactor.*-mean(ref{iconc}.flux_up_sw(1,:));
    flux_surf_ref(iconc) = dayfactor.*mean(ref{iconc}.flux_dn_sw(end,:)...
					   -ref{iconc}.flux_up_sw(end,:));
    
    for ie = ickd
      ckd_tmp = loadnc([ckd_base{ie} scenario '.nc']);
      if do_fix_ssi(ie)
	ckd_tmp = correct_ssi(ref_orig, ckd_tmp, specdef);
      end
      if do_flatten(ie)
	ckd{ie,iconc} = flatten_sza(ckd_tmp, isza);
      else
	ckd{ie,iconc} = ckd_tmp;
      end
      flux_toa_ckd{ie}(iconc)  = dayfactor.*-mean(ckd{ie,iconc}.flux_up_sw(1,:));
      flux_surf_ckd{ie}(iconc) = dayfactor.*mean(ckd{ie,iconc}.flux_dn_sw(end,:)...
						-ckd{ie,iconc}.flux_up_sw(end,:));
    end
  end

  subplot(2+plot_hr,ngas,igas)
  if plot_hr
    plot(concs{igas}(2:end-1), flux_toa_ref(2:end-1)-flux_toa_ref(ipresent), [col_ref '-o']);
  else
    plot(concs{igas}, flux_toa_ref-flux_toa_ref(ipresent), [col_ref '-o']);
  end
  hold on
  for ie = ickd
    plot(concs{igas}, flux_toa_ckd{ie}-flux_toa_ckd{ie}(ipresent), cols{ie});
  end
  plot(concs{igas}, flux_toa_ref-flux_toa_ref(ipresent), [col_ref '-']);
  if plot_hr
    plot(concs{igas}(1), flux_toa_ref(1)-flux_toa_ref(ipresent), [col_glacial 's']);
    plot(concs{igas}(end), flux_toa_ref(end)-flux_toa_ref(ipresent), [col_future 'd']);
  end
  if islog(igas)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas})
  %xlabel([gas_names{igas} ' volume mixing ratio (' units{igas} ')']);
  xlabel([gas_names{igas} ' mole fraction (' units{igas} ')']);
  if igas == 1
    legend(['Reference' titles],'location','north','box','off');
    %ylabel('TOA forcing w.r.t. present (W m^{-2})');
    ylabel('TOA forcing (W m^{-2})');
  end
  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  title(['\bf' gas_names{igas}]);
  axis(ax);
  text(0.025, 0.95, ['\bf(' 'a'-1+igas ')'],'units','normalized');

%  plot(etminan_conc{igas}, etminan_forcing{igas}, 'r-o');

  subplot(2+plot_hr,ngas,igas+ngas)
  if plot_hr
    plot(concs{igas}(2:end-1), flux_surf_ref(2:end-1)-flux_surf_ref(ipresent), [col_ref 'o-']);
  else
    plot(concs{igas}, flux_surf_ref-flux_surf_ref(ipresent), [col_ref 'o-']);
  end
  hold on
  for ie = ickd
    plot(concs{igas}, flux_surf_ckd{ie}-flux_surf_ckd{ie}(ipresent), cols{ie});
  end
  plot(concs{igas}, flux_surf_ref-flux_surf_ref(ipresent), [col_ref '-']);
  if plot_hr
    plot(concs{igas}(1), flux_surf_ref(1)-flux_surf_ref(ipresent), [col_glacial 's']);
    plot(concs{igas}(end), flux_surf_ref(end)-flux_surf_ref(ipresent), [col_future 'd']);
  end
  if islog(igas)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas})
  %xlabel([gas_names{igas} ' volume mixing ratio (' units{igas} ')']);
  xlabel([gas_names{igas} ' mole fraction (' units{igas} ')']);
  if igas == 1
    %ylabel('Surface forcing w.r.t. present (W m^{-2})');
    ylabel('Surface forcing (W m^{-2})');
  end
  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  axis(ax);
  text(0.025, 0.95, ['\bf(' 'a'-1+igas+ngas ')'],'units','normalized');

  if plot_hr
    subplot(3,ngas,igas+ngas*2)
    if length(ickd) > 1
      semilogy(-1,-1,'g-s');
      hold on
      plot(-1,-1,'m-d')
      for ie = ickd
	plot(-1,-1,cols{ie})
      end
    else
      semilogy(-1,-1,'g-s');
      hold on
      plot(-1,-1,'g--');
      plot(-1,-1,'m-d');
      plot(-1,-1,'m--');
    end

    phl = median(ref{1}.pressure_hl,2); % Half-level median pressure
    pfl = 0.5.*(phl(1:end-1) + phl(2:end)); % Full-level median pressure
    for ie = ickd
      for iconc = [1 length(concs{igas})]
	hr_ref = mean(calc_hr(ref{iconc},'sw'),1);
	hr_ckd = mean(calc_hr(ckd{ie,iconc},'sw'),1);
	hr_ref_minus_present = mean(calc_hr(ref{iconc},'sw')-calc_hr(ref{ipresent},'sw'),1);
	hr_ckd_minus_present = mean(calc_hr(ckd{ie,iconc},'sw')-calc_hr(ckd{ie,ipresent},'sw'),1);
	hr_err = hr_ckd - hr_ref;
	hr_err_minus_present = hr_ckd_minus_present - hr_ref_minus_present;
%	semilogy(hr_err_minus_present,pfl./100)
%	semilogy(hr_ref,pfl./100)
	if iconc == 1
	  semilogy(hr_ref_minus_present,pfl./100,col_glacial);
	  hold on
	  if length(ickd) > 1
	    semilogy(hr_ckd_minus_present,pfl./100,[cols{ie}]);
	  else
	    semilogy(hr_ckd_minus_present,pfl./100,[col_glacial '--']);
	  end
	else
	  semilogy(hr_ref_minus_present,pfl./100,col_future);
	  hold on
	  if length(ickd) > 1
	    semilogy(hr_ckd_minus_present,pfl./100,[cols{ie}]);
	  else
	    semilogy(hr_ckd_minus_present,pfl./100,[col_future '--']);
	  end
	end
	%semilogy(hr_ref_minus_present,pfl./100,col_ref);
%	hold on
%	semilogy(hr_ckd_minus_present,pfl./100,cols{ie});

	imarker = [6 20];
	if iconc == 1
	  semilogy(hr_ref_minus_present(imarker),pfl(imarker)./100,[col_glacial 's']);
	else
	  semilogy(hr_ref_minus_present(imarker),pfl(imarker)./100,[col_future 'd']);
	end
	plot([0 0],[0.01 1000],'k:','linewidth',0.5)
	set(gca,'ydir','reverse');
	ylim([0.01 1000]);
	xlabel('Heating rate change (K d^{-1})');
	ylabel('Pressure (hPa)');
	set(gca,'ytick',10.^[-2:3]);
	text(0.025, 0.95, ['\bf(' 'a'-1+igas+2.*ngas ')'],'units','normalized');
	xlim(hrxlim{igas});

      end
    end

    if igas == 1
      if length(ickd) > 1
	h=legend(['Reference min'],...
	       ['Reference max'],...
	       [titles{1}],[titles{2}],...
	       'location','southeast');
      else
	h=legend(['Reference min'],[titles{1} ' min'],...
	       ['Reference max'],[titles{1} ' max'],...
	       'location','southeast');
      end
      pos = get(h,'position')
      set(h,'position',pos+[0.04 0 0 0])
      %set(h,'fontsize',8,'box','on');
    end

  end

end
drawnow
%print_png([combined_titles 'evaluation1_forcing.png'],100)
