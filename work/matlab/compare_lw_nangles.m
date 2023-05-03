% Compare the longwave fluxes with different numbers of quadrature angles
BASEDIR='/hugetmp/parr/ckdmip/mmm/lw_fluxes';
CODE = {'fluxes','fluxes1ang','fluxes2ang','fluxes3ang','fluxes4ang','fluxes-8ang'};
PREFIX='ckdmip_mmm_lw_';
SUFFIX='_present_1-1.h5';
clear leg;
for ic = 1:length(CODE)
  data{ic} = loadnc([BASEDIR '/' PREFIX CODE{ic} SUFFIX]);
  hr{ic}   = calc_hr(data{ic},'lw');
  if length(CODE{ic}) < 7
    leg{ic} = 'Classic';
  else
    leg{ic} = [CODE{ic}(7) ' angles'];
  end
end
iref = length(CODE)

phl = data{1}.pressure_hl./100;
pfl = 0.5.*(phl(1:end-1) + phl(2:end));

myylim = [0.01 max(phl)];

cols = {'r','b','c','g','m','k'};

clf
subplot(2,3,1)
for ic = 1:length(data)
  semilogy(data{ic}.flux_dn_lw, phl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Downwelling longwave flux (W m^{-2})');
ylabel('Pressure (hPa)');
ylim(myylim);
legend(leg,'location','northeast');

subplot(2,3,2)
for ic = 1:length(data)
  semilogy(data{ic}.flux_up_lw, phl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Upwelling longwave flux (W m^{-2})');
ylabel('Pressure (hPa)');
ylim(myylim);

subplot(2,3,3)
for ic = 1:length(data)
  semilogy(hr{ic}, pfl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Heating rate (K d^{-1})');
ylabel('Pressure (hPa)');
ylim(myylim);

subplot(2,3,4)
for ic = 1:length(data)-1
  semilogy(data{ic}.flux_dn_lw-data{iref}.flux_dn_lw, phl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Error in downwelling longwave flux (W m^{-2})');
ylabel('Pressure (hPa)');
ylim(myylim);

subplot(2,3,5)
for ic = 1:length(data)-1
  semilogy(data{ic}.flux_up_lw-data{iref}.flux_up_lw, phl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Error in upwelling longwave flux (W m^{-2})');
ylabel('Pressure (hPa)');
ylim(myylim);

subplot(2,3,6)
for ic = 1:length(data)-1
  semilogy(hr{ic}-hr{iref}, pfl, cols{ic});
  hold on
  set(gca,'ydir','reverse');
end
xlabel('Error in heating rate (K d^{-1})');
ylabel('Pressure (hPa)');
ylim(myylim);
