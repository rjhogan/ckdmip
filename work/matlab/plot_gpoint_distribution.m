function plot_gpoint_distribution(file, default_bands)
    
  if nargin < 2
    default_bands = 0;
  end

  [d,a]=loadnc(file);
  ng = length(d.band_number);
  nband = length(d.wavenumber1_band);

  for iband = 1:length(d.wavenumber1_band)
    disp(['Band ' num2str(iband-1) ': ' num2str(d.wavenumber1_band(iband)) ...
	  '-' num2str(d.wavenumber2_band(iband)) ' cm-1, ' num2str(length(find(d.band_number==iband-1))) ' g points']);
  end

  if max(d.wavenumber2_band) > 10000
    is_sw = 1;
  else 
    is_sw = 0;
  end

  % Check if gpoint file available
  loc = strfind(a.global.history, 'output=');
  if ~isempty(loc)
     subs = a.global.history(loc(1)+7:end);
     loc = strfind(subs,' ');
     gpfile = subs(1:loc-1);
     disp(['Loading g-point file ' gpfile]);
     gp = loadnc(gpfile);
  end

  interval_size = num2str(median(d.wavenumber2-d.wavenumber1));

  if ~default_bands
     bb_wide = d.wavenumber1_band(2:end)';
     bb_narrow = [];
    if nband < 4
      cmax = 0.1
    elseif nband > 8
      cmax = 0.2;
    else
      cmax = 0.15;
    end
  elseif is_sw
    bb_wide = [4000, 8050, 16000, 29000];
    bb_narrow = [2600, 3250, 4650, 5150, 6150, 12850, 22650, 38000];
    %interval_size = '50';
    if nband < 4
      cmax = 0.1
    elseif nband > 8
      cmax = 0.2;
    else
      cmax = 0.15;
    end
  else
    bb_wide = [500, 820, 1180, 1800];
    bb_narrow = [0, 350, 630, 700, 980, 1080, 1390, 1480, 2080];
    %interval_size = '10';
    if nband < 4
      cmax = 0.2
    elseif nband > 8
      cmax = 0.4;
    else
      cmax = 0.3;
    end
  end

  clf
  set(gcf,'paperposition',[0.5 0.5 25 10]);
  axes('position',[0.15 0.2 0.65 0.7]) 
  data = d.gpoint_fraction([1:end end],[1:end end])';
  data(find(data <= 0)) = NaN;
  pcolor([d.wavenumber1;d.wavenumber2(end)], 0.5+[0:ng], data);
  shading flat
  chiljet
  caxis([0 cmax]);
  hold on
  plot([bb_wide;bb_wide],[0.5 ng+0.5]'*ones(1,length(bb_wide)),'k-','linewidth',1)
  plot([bb_narrow;bb_narrow],[0.5 ng+0.5]'*ones(1,length(bb_narrow)),'k:')
  xlabel('Wavenumber (cm^{-1})');
  ylabel('{\itk} term');
  set(gca,'layer','top')
  if is_sw
     xt=[250 5000:5000:50000];
    set(gca,'xtick',xt,'xticklabel',xt)
  end

  if exist('gp','var')
     gases = {'h2o','co2','ch4','n2o','o3','o2n2'};
     gasstr= {'H_2O','CO_2','CH_4','N_2O','O_3','O_2+N_2'};
     for igas = 1:length(gases)
       if isfield(gp,[gases{igas} '_g_max']);
	 loc = find(diff(gp.([gases{igas} '_g_max'])) ...
		    & ~diff(gp.band_number));
	 for ii = 1:length(loc)
	   text(d.wavenumber2(end), loc(ii)+1, [' ' gasstr{igas}])
	 end
       end
     end
  end

  %h=axes('position',[0.9 0.2 0.025 0.3]);
  h=colorbar('location','eastoutside');
  %set(h,'position',[0.77 0.20 0.02 0.5]);
  pp=get(h,'position');
  set(h,'position',[0.83+0.03 0.3 0.025 0.5]);%pp+[0.15 0.1 0 -0.2])
  h.Label.String = ['Fraction of {\itk} term in ' interval_size ' cm^{-1} interval'];
  h.FontSize=10;
  h.Label.FontSize=10;
