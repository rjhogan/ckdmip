% Convert CAVIAR text file to NetCDF
[data, attr, dims] = loadnc('../../data/wv-continuum_mt-ckd-4.1.1.nc');

use_caviar_foreign = 1;

cav_self296 = load('SOC_new_296K.dat');
wn=cav_self296(:,1);
cav_self260 = load('SOC_new_260K.dat');
cav_self296=cav_self296(:,2);
cav_self260=cav_self260(:,2);

clear *_new;

data_new = data;
data_new.wavenumbers = wn;


% Reproduce MT_CKD calculation of radiation term
radiation_term = wn;
xviokt = wn .* 1.4387752 / data.ref_temp;
index = find(xviokt <= 10.0);
expvkt = exp(-xviokt(index));
radiation_term(index) = wn(index).*(1.0-expvkt)./(1.0+expvkt);
index = find(xviokt <= 0.01);
radiation_term(index) = 0.5 * xviokt(index).*wn(index);
foreign_absco = data_new.for_absco_ref(3:end).*radiation_term;
foreign_absco_mtckd = foreign_absco;
isnew = zeros(size(foreign_absco));
if use_caviar_foreign
  % Overwrite with CAVIAR foreign continuum were available
  dataf1 = load('Fore_296+310K_1300-1800cm-1.dat');
  dataf2 = load('Fore_400K_1800-8970cm-1.dat');

  wnx = dataf1(:,1);
  index = find(wnx < 1800);
  for ii = index(1):index(end)
    index = find(abs(wn - wnx(ii)) < 5);
    if ~isempty(index)
      foreign_absco(index(1)) = dataf1(ii,2);
      isnew(index(1)) = 1;
    end
  end
  wnx = dataf2(:,1);
  index = find(wnx >= 1800);
  for ii = index(1):index(end)
    index = find(abs(wn - wnx(ii)) < 5);
    if ~isempty(index)
      foreign_absco(index(1)) = dataf2(ii,2);
      isnew(index(1)) = 2;
    end
  end
end

  
data_new.self_absco_ref = cav_self296;
data_new.self_texp = log(cav_self296./cav_self260) / log(260/296);
data_new.self_texp(find(~isfinite(data_new.self_texp))) = 0;
data_new.for_absco_ref = foreign_absco;
data_new.use_radiation_term = int8(0);
data_new.for_source = int8(isnew);

attr_new.wavenumbers.dimensions = {'wavenumbers'};
attr_new.wavenumbers.long_name = 'Wavenumber';
attr_new.wavenumbers.units = 'cm-1';

attr_new.self_absco_ref.dimensions = {'wavenumbers'};
attr_new.self_absco_ref.long_name = 'Water vapour self continuum absorption coefficient at 296 K';
attr_new.self_absco_ref.units = 'cm-2 molec-1';

attr_new.for_absco_ref.dimensions = {'wavenumbers'};
attr_new.for_absco_ref.long_name = 'Water vapour foreign continuum absorption coefficient at 296 K';
attr_new.for_absco_ref.units = 'cm-2 molec-1';

attr_new.self_texp.dimensions = {'wavenumbers'};
attr_new.self_texp.long_name = 'Temperature exponent for water vapour self continuum';
				%attr_new.self_texp.units = '1';

attr_new.ref_press = attr.ref_press;
attr_new.ref_press.units = 'hPa';

attr_new.ref_temp = attr.ref_temp;

attr_new.for_source.dimensions = {'wavenumbers'};
attr_new.for_source.long_name = 'Source of foreign continuum';
attr_new.for_source.definition = ['0: MT_CKD 4.1.1' 10 '1: CAVIAR measurements at 296 & 310 K' 10 '2: CAVIAR measurements at 400 K'];

attr_new.use_radiation_term.dimensions = {};
attr_new.use_radiation_term.long_name = 'Should the radiation term be used?';

%attr_new = attr;
attr_new.global.title='CAVIAR water vapour continuum';
%attr_new.self_absco_ref.units = 'cm**2/molecule';
%attr_new.for_absco_ref.units = 'cm**2/molecule';
attr_new.global.contact = 'Robin Hogan (r.j.hogan@ecmwf.int)';
attr_new.global.mtckd_copyright = attr.global.copyright;
attr_new.global.source = 'absco-ref_wv-mt-ckd.nc SOC_new_260K.dat SOC_new_296K.dat Fore_296+310K_1300-1800cm-1.dat Fore_400K_1800-8970cm-1.dat';
attr_new.global.comment = wrap_text('This file contains "CAVIAR" water vapour continuum absorption coefficients obtained from Keith Shine, with some filling of missing data using MT_CKD 4.1.1. The absorption coefficients in this file do not need to be multiplied by the radiation term.', 100);
dims_new = dims;
dims_new.wavenumbers = dims.wavenumbers-2;
write_nc_struct('wv-continuum_caviar.nc', dims_new, data_new, attr_new);

clf
subplot(2,1,1)
semilogy(data.wavenumbers(3:end),data.self_absco_ref(3:end).*radiation_term,'k')
hold on
semilogy(data_new.wavenumbers,data_new.self_absco_ref,'r--')
xlim([0 20000])
xlabel('Wavenumber (cm-1)')
ylabel('Reference absorption at 296 K (cm2 molec-1)')
legend('MT-CKD-4.1.1','CAVIAR','location','southwest')

subplot(2,1,2)
semilogy(data.wavenumbers,data.self_texp,'k')
hold on
semilogy(data_new.wavenumbers,data_new.self_texp,'r--')
xlim([0 20000])
xlabel('Wavenumber (cm-1)')
ylabel('Temperature exponent');
