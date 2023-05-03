function evaluate_ckd_lw_overlap(ref_base, ckd_base_in, model_name)

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

gases={'co2','ch4','n2o'};
gas_names = {'CO_2','CH_4','N_2O'};
units={'ppmv','ppbv','ppbv'};
concs={[180 415 2240],...
       [350 1921 3500],...
       [190 332 540]};
conc_axes = {[100 3000],[0 4000],[100 600]};
rf_axis = {[-3 6],[-1 1],[-1 1]};
rf_axis = {[-5 7],[-1.2 0.8],[-1 1]};
concs_present = [415 1921 332];
islog = [1 0 0];
my_xlim = [0.4 3.6];
my_xlim = [0 4];
barcols = {'k','r','b'};
pspace = [0.1 0.02 0.02];
%ibar = 1:3;
ibar = [1 3];
is_horiz = 1;


ngas = length(gases)

%nx = ngas; ny = 2;
ny = ngas; nx = 2;

nckd = length(ckd_base);
titles = titles(1:nckd);
if nckd == 1
  textloc = 3+([1 2]-1.5).*0.5;
else
  textloc = 3+([1:nckd]-2)*0.45;
end


ickd = 1:nckd
col_ref = 'k';
cols = {'r--','b-.','g','m','c'};
ngas = length(gases);
combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end

gas1 = [1 1 2 2 3 3];
gas2 = [2 3 1 3 1 2];

ncase = length(gas1)

clf
set(gcf,'defaultlinelinewidth',1);
set(gcf,'paperposition',[0.5 0.5 27 15]);
set(gcf,'paperposition',[0.5 0.5 27 19]);

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
	  if igas1 > igas2
	    scenario = [str2 '-' str1];
	  else
	    scenario = [str1 '-' str2];
	  end
	end
      end
      ref = loadnc([ref_base scenario '.h5']);
      flux_toa_ref(iconc1,iconc2)  = -mean(ref.flux_up_lw(1,:));
      flux_surf_ref(iconc1,iconc2) = mean(ref.flux_dn_lw(end,:));
      for ie = ickd
	ckd = loadnc([ckd_base{ie} scenario '.nc']);
	flux_toa_ckd{ie}(iconc1,iconc2)  = -mean(ckd.flux_up_lw(1,:));
	flux_surf_ckd{ie}(iconc1,iconc2) = mean(ckd.flux_dn_lw(end,:));
      end
    end
  end

  ibaseconc = 2; % present
  %ibaseconc = 1; % glacial

  subplot(ny,nx,igas)
  bars(1:3,1) = flux_toa_ref(3,:)-flux_toa_ref(ibaseconc,:)
  bars_dn(1:3,1) = flux_toa_ref(1,:)-flux_toa_ref(ibaseconc,:)
  for ie = ickd
    bars(1:3,ie+1) = flux_toa_ckd{ie}(3,:)-flux_toa_ckd{ie}(ibaseconc,:);
    bars_dn(1:3,ie+1) = flux_toa_ckd{ie}(1,:)-flux_toa_ckd{ie}(ibaseconc,:);
  end

  if is_horiz
    h=barh(-ibar,-10.*ones(length(ibar),1+nckd));
    hold on
    for ii = 1:length(h)
      h(ii).FaceColor = barcols{2+nckd-ii};
    end

    %bar(concs{igas2},bars)
    h=barh(ibar,bars(ibar,:));
    perc = 100.*(bars(ibar(end),:)-bars(ibar(1),:))./bars(ibar(1),:);
    for ii = 1:length(h)
      h(ii).FaceColor = barcols{ii};
      hold on
      text(bars(ibar(end),ii)+pspace(gas1(igas)), textloc(ii), num2str(perc(ii),'%+0.1f%%'),'color',barcols{ii})
    end
    set(gca,'ytick',ibar,'yticklabel',concs{igas2}(ibar))
    hold on
    h=barh(ibar,bars_dn(ibar,:));
    perc = 100.*(bars_dn(ibar(end),:)-bars_dn(ibar(1),:))./bars_dn(ibar(1),:);
    for ii = 1:length(h)
      h(ii).FaceColor = barcols{ii};
      hold on
      text(bars_dn(ibar(end),ii)-pspace(gas1(igas)), textloc(ii), num2str(perc(ii),'%+0.1f%%'),...
	  'horizontalalignment','right','color',barcols{ii})

    end
    ylabel([gas_names{igas2} ' conc. (' units{igas2} ')']);
    xlabel([gas_names{igas1} ' radiative forcing (W m^{-2})']);
    xlim(rf_axis{igas1})
    ylim(my_xlim);
    %set(gca,'ydir','reverse');
    text(0.015,0.92,['\bf(' 'a'-1+igas ')'],'units','normalized');
  else
    %bar(concs{igas2},bars)
    h=bar(ibar,bars(ibar,:));
    for ii = 1:length(h)
      h(ii).FaceColor = barcols{ii};
    end
    set(gca,'xticklabel',concs{igas2}(ibar))
    hold on
    h=bar(ibar,bars_dn(ibar,:));
    for ii = 1:length(h)
      h(ii).FaceColor = barcols{ii};
    end
    xlabel([gas_names{igas2} ' conc. (' units{igas2} ')']);
    ylabel([gas_names{igas1} ' radiative forcing (W m^{-2})']);
    ylim(rf_axis{igas1})
    xlim(my_xlim);
  end

  if igas == 1
    %legend('LBLRTM','RRTMG','RRTMGP','location','south');
    legend(titles{1},'Reference','location','south');
  end 

end
drawnow
%print_png([combined_titles 'evaluation1_forcing.png'],100)
