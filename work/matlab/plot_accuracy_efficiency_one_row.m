function metrics = plot_accuracy_efficiency(model, model_legend, domain, application, ...
					    bandstructures, ngpoints,dataset)

  ndomain = 2;

  idomain = 1;
  if ndomain == 2
    if strcmp(domain,'sw')
      idomain = 2;
    end
  end

  ckdmip_dir = '/hugetmp/parr/ckdmip';
  ckd_dir = [ckdmip_dir '/results/' model '/' domain '_fluxes'];

  if strcmp(application, 'climate')
    scenarios = {'present','preindustrial','future','glacialmax',...
		 'co2-180','co2-280','co2-560','co2-1120','co2-2240',...
		 'ch4-350','ch4-700','ch4-1200','ch4-2600','ch4-3500',...
		 'n2o-190','n2o-270','n2o-405','n2o-540'};
    %scenarios = {'present','preindustrial','future','glacialmax'};
    %scenarios = {'present','future'};
    %scenarios = {'present'};
  else
    scenarios = {'present'};
  end

  if strcmp(domain,'lw')
    rt_mode = 'fluxes-4angle';
    Domain = 'Longwave';
  else
    rt_mode = 'fluxes';
    Domain = 'Shortwave';
  end

  for iband = 1:length(bandstructures)

    ref_dir = [ckdmip_dir '/' dataset{iband} '/' domain '_fluxes'];
    for iscen = 1:length(scenarios)
      ref_files{iscen} = [ref_dir '/ckdmip_' dataset{iband} '_' domain '_' rt_mode '_' scenarios{iscen} '.h5'];
    end

    bs = bandstructures{iband};
    if strcmp(bs,'fsck')
      leg{iband} = 'FSCK';
    else
      leg{iband} = [upper(bs(1)) bs(2:end)];
    end

    for ig = 1:length(ngpoints)
      if strcmp(bs,'rgb') & ngpoints(ig,iband) == 8
	bs_tmp = 'double';
      else
	bs_tmp = bs;
      end

      for iscen = 1:length(scenarios)
	model_id = [bs_tmp '-' num2str(ngpoints(ig,iband))];
	ckd_files{iscen} = [ckd_dir '/' model '_' dataset{iband} '_' domain '_' application '_' model_id '_' rt_mode '_' scenarios{iscen} '.nc'];
      end

      if strcmp(domain,'lw')
	[metrics(:,ig,iband), metric_names] = calc_error_metrics_lw(ref_files, ckd_files);
      else
	[metrics(:,ig,iband), metric_names] = calc_error_metrics_sw(ref_files, ckd_files);
      end
    end
  end

%  clf
  set(gcf,'paperposition',[0.5 0.5 32 10+8*(ndomain-1)]);
  nmet = size(metrics,1);

  metmap = [1 2 6 3 4 5];

  if length(bandstructures) < 3
    %cols = 'gb';
    stys = {'--','-'};
  else
    %cols = 'rgb';
    stys = {'--','-','s:'};
  end
  
  maxx = ceil(max(ngpoints(:))./10).*10;

  imet = 1;
  ylabs = {[Domain ' irradiance bias (W m^{-2})'],...
	   [Domain ' irradiance RMSE (W m^{-2})'],...
	   [Domain ' heating rate RMSE (K d^{-1})']};
  dataset_str = {'(training)','(independent)'};
  met_names = {'TOA up','TOA up','Surface down','Surface down','4-1100 hPa','0.02-4 hPa'};

  for ipan = 1:nmet./2
    subplot(ndomain,3,ipan+(idomain-1)*3)
    for ibs = 1:length(bandstructures)
      plot(ngpoints(:,ibs), metrics(metmap(imet),:,ibs),['b' stys{ibs}],'linewidth',1)
      hold on
      leg{ibs}=[met_names{metmap(imet)} ' ' dataset_str{ibs}];
    end
    for ibs = 1:length(bandstructures)
      plot(ngpoints(:,ibs), metrics(metmap(imet+3),:,ibs),['r' stys{ibs}],'linewidth',1)
      leg{ibs+2}    =[met_names{metmap(imet+3)} ' ' dataset_str{ibs}];
    end
    %ylabel(metric_names{metmap(imet)});
    ylabel(ylabs{ipan})
    xlabel('Number of {\itk} terms');
    set(gca,'xtick',[0:16:64]);
    hold on
    plot([0 maxx],[0 0],'k-','linewidth',0.5);
    xlim([0 maxx]);
    %if imet == 2
    legend(leg,'location','best');
    %end
    text(1-0.975,0.95,['\bf(' 'a'-1+imet+(idomain-1)*3 ')'],'units','normalized','horizontalalignment','left');
    imet = imet+1;
    grid on
  end

metric_names
