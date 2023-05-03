function evaluate_ckd2(model, model_legend, domain, application, setting, plot_type, specdef_setting, setting2)
% Matlab function to produce evaluation plots for CKDMIP
% Usage:
%   evaluate_ckd(model, model_legend, domain, application, setting, plot_type)
% e.g.:
%   evaluate_ckd ecrad-rrtmg ecRad-RRTMG lw climate narrow-140 png

dataset = 'evaluation2';

ref_dir = ['/hugetmp/parr/ckdmip/' dataset '/' domain '_fluxes'];
ckd_top_dir = ['/hugetmp/parr/ckdmip/results/' model];
ckd_dir = [ckd_top_dir '/' domain '_fluxes'];
dest_dir= ['/home/rd/parr/src/ckdmip-plots'];

read_specdef = 0;
plot_band_error = 0;
if nargin < 5
  help evaluate_ckd
  return
elseif nargin < 6
  plot_type = '';
end
if nargin < 7
  specdef_setting = setting
end

if strcmp(application, 'climate')
  scenarios = {'present','preindustrial','future','glacialmax'}; 
  scenario_titles = {'Present-day (2020)','Preindustrial','Future','Glacial Maximum'}; 
else
  scenarios = {'present'};
  scenario_titles = {'Present'};
end
scenarios = {'present'};

if strcmp(domain,'lw')
%  fluxscheme = 'fluxes'; % Classic 2-stream
  fluxscheme = 'fluxes-4angle'; % 8-stream
else
  fluxscheme = 'fluxes';
end

domain_suffix = '';
if strcmp(setting(1:4),'psla')
  domain_suffix = '-pslackd20'
elseif strcmp(setting(1:3),'rgb') | strcmp(setting(1:2),'gb')
  domain_suffix = '-rgb';
end

if read_specdef
specdef_file = [ckd_top_dir '/' domain '_spectral-definition/' model '_' domain '_' application '_' specdef_setting '_spectral-definition.nc'];

%specdef = loadnc([ckd_top_dir '/' model '_' domain '_' application '_' setting '_spectral-definition.nc']);
specdef = loadnc(specdef_file);

specdef.setting=setting;
end

%%% EVALUATE BROADBAND FLUX PROFILES AND HEATING RATES
for iscen = 1:length(scenarios)
  figure(1)
  scenario = scenarios{iscen};
  if strcmp(domain, 'lw')
    evaluate_ckd_lw_fluxes([ref_dir '/ckdmip_' dataset '_' domain '_' fluxscheme '_' scenario '.h5'], ...
			   {[ckd_dir '/' model '_' dataset '_' domain '_' ...
				    application '_' setting '_' fluxscheme '_' scenario '.nc'],...
			   [ckd_dir '/' model '_' dataset '_' domain '_' ...
				    application '_' setting2 '_' fluxscheme '_' scenario '.nc']},...
			   model_legend, scenario_titles{iscen});
  else
    evaluate_ckd_sw_fluxes([ref_dir '/ckdmip_' dataset '_' domain '_' fluxscheme '_' scenario '.h5'], ...
			   {[ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_' scenario '.nc'],...
			    [ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting2 '_' fluxscheme '_' scenario '.nc']},...
			   model_legend, scenario_titles{iscen});
  end
  plot_base = [dest_dir '/' model '_' dataset '_' domain '_' application '_' setting '_' fluxscheme '_' scenario];
  if strcmp(plot_type, 'png')
    print_png([plot_base '.png'],'120');
  elseif strcmp(plot_type, 'pdf');
    print_pdf([plot_base '.pdf']);
  end
end

%%% EVALUATE FLUXES AND HEATING RATES IN BANDS
if 0

for iscen = 1:length(scenarios)
  if ~strcmp(setting(1:4),'fsck')
    figure(4)
    scenario = scenarios{iscen};
    plot_base = [dest_dir '/' model '_' dataset '_' domain '_' application '_' setting '_bands_' scenario];
    if strcmp(domain, 'lw')
      evaluate_ckd_lw_bands([ref_dir domain_suffix '/ckdmip_' dataset '_' domain '_' fluxscheme domain_suffix '_' scenario '.h5'], ...
			    [ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_' scenario '.nc'],...
			    model_legend, scenario_titles{iscen}, specdef);
    else

      evaluate_ckd_sw_bands([ref_dir domain_suffix '/ckdmip_' dataset '_' domain '_' fluxscheme domain_suffix '_' scenario '.h5'], ...
			    [ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_' scenario '.nc'],...
			    model_legend, scenario_titles{iscen}, specdef);
      
    end
    if strcmp(plot_type, 'png')
      print_png([plot_base '.png'],'120');
    elseif strcmp(plot_type, 'pdf');
      print_pdf([plot_base '.pdf']);
    end
  end
end

end

if plot_band_error

for iscen = 1:length(scenarios)
  if ~strcmp(setting(1:4),'fsck')
    figure(6)
    scenario = scenarios{iscen};
    plot_base = [dest_dir '/' model '_' dataset '_' domain '_' application '_' setting '_bands-err_' scenario];
    if strcmp(domain, 'lw')
      evaluate_ckd_lw_bands_err([ref_dir domain_suffix '/ckdmip_' dataset '_' domain '_' fluxscheme domain_suffix '_' scenario '.h5'], ...
			    [ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_' scenario '.nc'],...
			    model_legend, scenario_titles{iscen}, specdef);
    else

      evaluate_ckd_sw_bands_err([ref_dir domain_suffix '/ckdmip_' dataset '_' domain '_' fluxscheme domain_suffix '_' scenario '.h5'], ...
				[ckd_dir '/' model '_' dataset '_' domain '_' ...
					 application '_' setting '_' fluxscheme '_' scenario '.nc'],...
				model_legend, scenario_titles{iscen}, specdef);
    end
    if strcmp(plot_type, 'png')
      print_png([plot_base '.png'],'120');
    elseif strcmp(plot_type, 'pdf');
      print_pdf([plot_base '.pdf']);
    end
  end
end

end

%%% EVALUATE RADIATIVE FORCING OF GREENHOUSE GASES
if strcmp(application,'climate')
  figure(2)
  if strcmp(domain, 'lw')
    evaluate_ckd_lw_forcing([ref_dir '/ckdmip_' dataset '_' domain '_' fluxscheme '_'], ...
			    {[ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_'],...
			    [ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting2 '_' fluxscheme '_']},...
			    model_legend);
  else
    evaluate_ckd_sw_forcing([ref_dir '/ckdmip_' dataset '_' domain '_' fluxscheme '_'], ...
			    {[ckd_dir '/' model '_' dataset '_' domain '_' ...
				     application '_' setting '_' fluxscheme '_'],...
			     [ckd_dir '/' model '_' dataset '_' domain '_' ...
				      application '_' setting2 '_' fluxscheme '_']},...
			    model_legend);

  end
  plot_base = [dest_dir '/' model '_' dataset '_' domain '_' application '_' setting '_forcing'];
  if strcmp(plot_type, 'png')
    print_png([plot_base '.png'],'120');
  elseif strcmp(plot_type, 'pdf');
    print_pdf([plot_base '.pdf']);
  end
end

%%% EVALUATE FORCING OVERLAP OF GREENHOUSE GASES
if strcmp(application,'climate') & strcmp(domain, 'lw')
  figure(3);
  evaluate_ckd_lw_overlap([ref_dir '/ckdmip_' dataset '_' domain '_' fluxscheme '_'], ...
			  [ckd_dir '/' model '_' dataset '_' domain '_' ...
				   application '_' setting '_' fluxscheme '_'],...
			  model_legend);
  plot_base = [dest_dir '/' model '_' dataset '_' domain '_' application '_' setting '_overlap'];
  if strcmp(plot_type, 'png')
    print_png([plot_base '.png'],'120');
  elseif strcmp(plot_type, 'pdf');
    print_pdf([plot_base '.pdf']);
  end
end

if read_specdef
  figure(5)
  plot_gpoint_distribution(specdef_file);
  plot_base = [dest_dir '/' model '_' domain '_' application '_' setting '_spectral-definition'];
  if strcmp(plot_type, 'png')
    print_png([plot_base '.png'],'150');
  elseif strcmp(plot_type, 'pdf');
    print_pdf([plot_base '.pdf']);
  end
end
