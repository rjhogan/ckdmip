function metrics = plot_accuracy_efficiency(model, model_legend, domain, applications, ...
					    bandstructure, ngpoints, maxg)

  dataset = 'evaluation1';
  ckdmip_dir = '/hugetmp/parr/ckdmip';
  ref_dir = [ckdmip_dir '/' dataset '/' domain '_fluxes'];
  ckd_dir = [ckdmip_dir '/results/' model '/' domain '_fluxes'];

%    scenarios = {'present','preindustrial','future','glacialmax',...
%		 'co2-180','co2-280','co2-560','co2-1120','co2-2240',...
%		 'ch4-350','ch4-700','ch4-1200','ch4-2600','ch4-3500',...
%		 'n2o-190','n2o-270','n2o-405','n2o-540'};
  scenarios = {'present'};

  if strcmp(domain,'lw')
    rt_mode = 'fluxes-4angle';
  else
    rt_mode = 'fluxes';
  end

  for iscen = 1:length(scenarios)
    ref_files{iscen} = [ref_dir '/ckdmip_' dataset '_' domain '_' rt_mode '_' scenarios{iscen} '.h5'];
  end

  for iapp = 1:length(applications)
    application = applications{iapp}
    
    leg{iapp} = [upper(application(1)) application(2:end)];

    for ig = 1:size(ngpoints,1)
      for iscen = 1:length(scenarios)
	model_id = [bandstructure '-' num2str(ngpoints(ig,iapp))];
	ckd_files{iscen} = [ckd_dir '/' model '_' dataset '_' domain '_' application '_' model_id '_' rt_mode '_' scenarios{iscen} '.nc'];
      end

      ckd_files

      if strcmp(domain,'lw')
	[metrics(:,ig,iapp), metric_names] = calc_error_metrics_lw(ref_files, ckd_files);
      else
	[metrics(:,ig,iapp), metric_names] = calc_error_metrics_sw(ref_files, ckd_files);
      end
    end
  end

  clf
  set(gcf,'paperposition',[0.5 0.5 25 19]);
  nmet = size(metrics,1);

  metmap = [1 2 6 3 4 5];

  cols = 'rgb';
  %cols = 'gb';
  
  maxx = ceil(max(ngpoints(:))./10).*10;
  if exist('maxg','var')
    maxx = maxg;
  end

  for imet = 1:nmet
    subplot(2,3,imet)
    for iapp = 1:length(applications)
      plot(ngpoints(:,iapp), metrics(metmap(imet),:,iapp),[cols(iapp) '-o'],'linewidth',1.5)
      hold on
    end
    ylabel(metric_names{metmap(imet)});
    xlabel('Number of {\itk} terms');
    hold on
    plot([0 maxx],[0 0],'k:');
    xlim([0 maxx]);
    if imet == 2
      legend(leg,'location','north');
    end
    text(0.975,0.95,['\bf(' 'a'-1+imet ')'],'units','normalized','horizontalalignment','right');
    grid on
  end

