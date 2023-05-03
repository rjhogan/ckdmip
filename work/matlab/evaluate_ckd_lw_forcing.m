function evaluate_ckd_lw_forcing(ref_base, ckd_base_in, model_name)

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

invert_subplots = 0;

plot_hr = 0;

gases={'co2','ch4','n2o','cfc11','cfc12'};
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

ngas = length(gases)

ickd = 1:length(ckd_base);
col_ref = 'k';
%col_future = 'm';
col_future = [1 0 1];
%col_future = [240, 128, 128]./255;
%col_future = [186, 85, 211]./255;
col_future = [0.9 0 0.9];
%col_glacial = 'g';
col_glacial = [0 1 0];
%col_glacial = [102, 221, 170]./255;
col_glacial = [50, 205, 50]./255;

cols = {'r--','b-.','g','m','c'};
if length(ickd) > 1
  cols = {'b-.','r--','g','m','c'};
end
ngas = length(gases);
combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end

clf
set(gcf,'defaultlinelinewidth',1);
%set(gcf,'paperposition',[0.5 0.5 45 20+plot_hr.*6]);
if invert_subplots
  set(gcf,'paperposition',[0.5 0.5 18+plot_hr.*5 ngas*7]);
  set(gcf,'paperposition',[0.5 0.5 18+plot_hr.*6.5 ngas*7]);
else
  set(gcf,'paperposition',[0.5 0.5 ngas*9 18+plot_hr.*5]);
  % Compressed, 13 Aug 2021
  set(gcf,'paperposition',[0.5 0.5 ngas*7 14+plot_hr.*4]);
end
for igas = 1:ngas

  clear flux_*ref flux_*ckd ipresent

  for iconc = 1:length(concs{igas})
    conc = concs{igas}(iconc);
    scenario = [gases{igas} '-' num2str(conc)];
    if conc == concs_present(igas)
      scenario = 'present';
      ipresent = iconc;
    end
    ref{iconc} = loadnc([ref_base scenario '.h5']);
    flux_toa_ref(iconc)  = -mean(ref{iconc}.flux_up_lw(1,:));
    flux_surf_ref(iconc) = mean(ref{iconc}.flux_dn_lw(end,:));
    
    for ie = ickd
      ckd{ie,iconc} = loadnc([ckd_base{ie} scenario '.nc']);
      flux_toa_ckd{ie}(iconc)  = -mean(ckd{ie,iconc}.flux_up_lw(1,:));
      flux_surf_ckd{ie}(iconc) = mean(ckd{ie,iconc}.flux_dn_lw(end,:));
    end
  end

  if invert_subplots
    ipanel = (igas-1)*(2+plot_hr)+1;
    subplot(ngas,2+plot_hr,ipanel)
  else
    ipanel = igas;
    subplot(2+plot_hr,ngas,ipanel)
  end
  if plot_hr
    plot(concs{igas}(2:end-1), flux_toa_ref(2:end-1)-flux_toa_ref(ipresent), [col_ref 'o-']);
  else
    plot(concs{igas}, flux_toa_ref-flux_toa_ref(ipresent), [col_ref 'o-']);
  end
  hold on
  for ie = ickd
    plot(concs{igas}, flux_toa_ckd{ie}-flux_toa_ckd{ie}(ipresent), cols{ie});
  end
  plot(concs{igas}, flux_toa_ref-flux_toa_ref(ipresent), [col_ref '-']);
  if plot_hr
    set(plot(concs{igas}(1), flux_toa_ref(1)-flux_toa_ref(ipresent), ['ks']),'color',col_glacial);
    set(plot(concs{igas}(end), flux_toa_ref(end)-flux_toa_ref(ipresent), ['kd']),'color',col_future);
  end

  if islog(igas)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas})
  xlabel([gas_names{igas} ' mole fraction (' units{igas} ')']);
  if invert_subplots
    ylabel(['{\bf' gas_names{igas} '}' 10 'TOA forcing w.r.t. present (W m^{-2})']);
  else
    if igas == 1
      %ylabel('TOA forcing w.r.t. present (W m^{-2})');
      ylabel('TOA forcing (W m^{-2})');
    end
    title(['\bf' gas_names{igas}]);
  end

  if igas == 1
    legend(['Reference' titles],'location','north','box','off');
  end
  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  axis(ax);
  text(0.025, 0.95, ['\bf(' 'a'-1+ipanel ')'],'units','normalized');

  if invert_subplots
     ipanel = (igas-1)*(2+plot_hr)+2;
    subplot(ngas,2+plot_hr,ipanel)
  else
    ipanel = igas+ngas;
    subplot(2+plot_hr,ngas,ipanel)
  end
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
    set(plot(concs{igas}(1), flux_surf_ref(1)-flux_surf_ref(ipresent), ['ks']),'color',col_glacial);
    set(plot(concs{igas}(end), flux_surf_ref(end)-flux_surf_ref(ipresent), ['kd']),'color',col_future);
  end

  if islog(igas)
    set(gca,'xscale','log');
  end
  xlim(conc_axes{igas})
  xlabel([gas_names{igas} ' mole fraction (' units{igas} ')']);
  if invert_subplots
    ylabel('Surf. forcing w.r.t. present (W m^{-2})');
  else
    %ylabel('Surf. forcing w.r.t. present (W m^{-2})');
    if igas == 1
      ylabel('Surface forcing (W m^{-2})');
    end
  end
  ax = axis;
  plot(concs{igas}(ipresent)+[0 0],ylim,'k:','linewidth',0.5);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  axis(ax);
  text(0.025, 0.95, ['\bf(' 'a'-1+ipanel ')'],'units','normalized');

  if plot_hr
    if invert_subplots
      ipanel = igas*3;
      subplot(ngas,3,ipanel)
    else
      ipanel = igas+ngas*2;
      subplot(3,ngas,ipanel)
    end
    set(semilogy(-1,-1,'-s'),'color',col_glacial);
    hold on
    set(plot(-1,-1,'--'),'color',col_glacial);
    set(plot(-1,-1,'-d'),'color',col_future);
    set(plot(-1,-1,'--'),'color',col_future);

    phl = median(ref{1}.pressure_hl,2); % Half-level median pressure
    pfl = 0.5.*(phl(1:end-1) + phl(2:end)); % Full-level median pressure
    hrscale = 1;
    if igas > 3
       hrscale = 1000;
    end
    for ie = ickd
      for iconc = [1 length(concs{igas})]
	hr_ref = mean(calc_hr(ref{iconc},'lw'),1);
	hr_ckd = mean(calc_hr(ckd{ie,iconc},'lw'),1);
	hr_ref_minus_present = mean(calc_hr(ref{iconc},'lw')-calc_hr(ref{ipresent},'lw'),1);
	hr_ckd_minus_present = mean(calc_hr(ckd{ie,iconc},'lw')-calc_hr(ckd{ie,ipresent},'lw'),1);
	hr_err = hr_ckd - hr_ref;
	hr_err_minus_present = hr_ckd_minus_present - hr_ref_minus_present;
%	semilogy(hr_err_minus_present,pfl./100)
%	semilogy(hr_ref,pfl./100)
	if iconc == 1
	  set(semilogy(hr_ref_minus_present.*hrscale,pfl./100,'-'),'color',col_glacial);
	  hold on
	  if length(ickd) > 1
	    semilogy(hr_ckd_minus_present.*hrscale,pfl./100,cols{ie})
	  else
	    set(semilogy(hr_ckd_minus_present.*hrscale,pfl./100,'--'),'color',col_glacial)
	  end
	else
	  set(semilogy(hr_ref_minus_present.*hrscale,pfl./100,'-'),'color',col_future);
	  hold on
	  if length(ickd) > 1
	    semilogy(hr_ckd_minus_present.*hrscale,pfl./100,cols{ie});
	  else
	    set(semilogy(hr_ckd_minus_present.*hrscale,pfl./100,'--'),'color',col_future);
	  end
	end

	imarker = [6 20];
	if iconc == 1
	  set(semilogy(hr_ref_minus_present(imarker).*hrscale,pfl(imarker)./100, ['ks']),'color',col_glacial);
	else
	  set(semilogy(hr_ref_minus_present(imarker).*hrscale,pfl(imarker)./100, ['kd']),'color',col_future);
	end
	plot([0 0],[0.01 1000],'k:','linewidth',0.5)
	set(gca,'ydir','reverse','ytick',10.^[-2:3]);
	ylim([0.01 1000]);
	if hrscale > 1
	  xlabel('Heating rate change (mK d^{-1})');
	else
	  xlabel('Heating rate change (K d^{-1})');
	end
	ylabel('Pressure (hPa)');
      end
    end
    text(0.025, 0.95, ['\bf(' 'a'-1+ipanel ')'],'units','normalized');
    if igas == 1
      if invert_subplots
	legend(['Ref min'],[titles{1} ' min'],...
	       ['Ref max'],[titles{1} ' max'],...
	       'location','southeast');
      elseif length(ickd) > 1
	legend(['Reference min'],...
	       ['Reference max'],[titles{1}],[titles{2}],...
	       'location','southwest');
      else
	legend(['Reference min'],[titles{1} ' min'],...
	       ['Reference max'],[titles{1} ' max'],...
	       'location','northeast');
      end
    end
  end


end

%print_png([combined_titles 'evaluation1_forcing.png'],100)
