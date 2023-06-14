% Plot the coefficients from various water vapour continuum models
DIR='../../data';
MODELS={'mt-ckd-2.5','mt-ckd-3.2','mt-ckd-4.1.1'};
for imod=1:length(MODELS)
  data{imod}=loadnc([DIR '/wv-continuum_' MODELS{imod} '.nc']);
end
vars = {'self_absco_ref','for_absco_ref','self_texp'};
nvar = length(vars);

stys = {'k','g--','r-.'};

for pr=1:2

  figure(pr)
  clf
  for ivar = 1:nvar
    subplot(nvar,1,ivar)
    for imod = pr:length(MODELS)
      if pr == 1
	semilogy(data{imod}.wavenumbers, data{imod}.(vars{ivar}),stys{imod});
      else
	semilogy(data{imod}.wavenumbers, data{imod}.(vars{ivar})./data{1}.(vars{ivar}),stys{imod});
      end
      hold on
    end
    xlim([0 20000]);
    legend(MODELS{pr:end});
    h=ylabel(vars{ivar});
    set(h,'interpreter','none')
    if ivar == 1
      if pr == 1
	title('Water vapour continuum coefficients');
      else
	title(['Water vapour continuum coefficients divided by ' MODELS{1}]);
      end
    end
  end
  xlabel('Wavenumber (cm-1)')
  set(gcf,'paperposition',[0.5 0.5 20 20])
end
