function metrics = plot_accuracy_efficiency(model, model_legend, domain, application, ...
					    bandstructures, ngpoints)

  dataset = 'evaluation1';
  ckdmip_dir = '/hugetmp/parr/ckdmip';
  ref_dir = [ckdmip_dir '/' dataset '/' domain '_fluxes'];
  ckd_dir = [ckdmip_dir '/results/' model '/' domain '_fluxes'];

  if strcmp(application, 'climate')
    scenarios = {'present','preindustrial','future','glacialmax',...
		 'co2-180','co2-280','co2-560','co2-1120','co2-2240',...
		 'ch4-350','ch4-700','ch4-1200','ch4-2600','ch4-3500',...
		 'n2o-190','n2o-270','n2o-405','n2o-540'};
    scenarios = {'present','preindustrial','future','glacialmax'}
    scenarios = {'present'};
  else
    scenarios = {'present'};
  end

  if strcmp(domain,'lw')
    rt_mode = 'fluxes-4angle';
  else
    rt_mode = 'fluxes';
  end

  for iscen = 1:length(scenarios)
    ref_files{iscen} = [ref_dir '/ckdmip_' dataset '_' domain '_' rt_mode '_' scenarios{iscen} '.h5'];
  end

  for iband = 1:length(bandstructures)
    bs = bandstructures{iband};
    if strcmp(bs,'fsck')
      leg{iband} = 'FSCK';
    else
      leg{iband} = [upper(bs(1)) bs(2:end)];
    end

    for ig = 1:length(ngpoints)
      for iscen = 1:length(scenarios)
	model_id = [bandstructures{iband} '-' num2str(ngpoints(ig,iband))];
	ckd_files{iscen} = [ckd_dir '/' model '_' dataset '_' domain '_' application '_' model_id '_' rt_mode '_' scenarios{iscen} '.nc'];
      end

      if strcmp(domain,'lw')
	[metrics(:,ig,iband), metric_names] = calc_error_metrics_lw(ref_files, ckd_files);
      else
	[metrics(:,ig,iband), metric_names] = calc_error_metrics_sw(ref_files, ckd_files);
      end
    end
  end

  clf
  set(gcf,'paperposition',[0.5 0.5 25 19]);
  nmet = size(metrics,1);

  metmap = [1 2 6 3 4 5];

  if length(bandstructures) < 3
    cols = 'gb';
  else
    cols = 'rgb';
  end
  
  maxx = ceil(max(ngpoints(:))./10).*10;

  for imet = 1:nmet
    subplot(2,3,imet)
    for ibs = 1:length(bandstructures)
      plot(ngpoints(:,ibs), metrics(metmap(imet),:,ibs),[cols(ibs) '-o'],'linewidth',1.5)
      hold on
    end
    ylabel(metric_names{metmap(imet)});
    xlabel('Number of {\itk} terms');
    hold on
    plot([0 maxx],[0 0],'k:');
    xlim([0 maxx]);
    if imet == 2
       legend(leg,'location','best');
    end
    text(0.975,0.95,['\bf(' 'a'-1+imet ')'],'units','normalized','horizontalalignment','right');
  end

