% Convert CAVIAR text file to NetCDF
[data, attr, dims] = loadnc('../../data/wv-continuum_mt-ckd-4.1.1.nc');

cav_self296 = load('SOC_new_296K.dat');
wn=cav_self296(:,1);
cav_self260 = load('SOC_new_260K.dat');
cav_self296=cav_self296(:,2);
cav_self260=cav_self260(:,2);

data_new = data;
data_new.wavenumbers = wn;

% Reproduce MT_CKD calculation
radiation_term = wn;
xviokt = wn .* 1.4387752 / data.ref_temp;
index = find(xviokt <= 10.0);
expvkt = exp(-xviokt(index));
radiation_term(index) = wn(index).*(1.0-expvkt)./(1.0+expvkt);
index = find(xviokt <= 0.01);
radiation_term(index) = 0.5 * xviokt(index).*wn(index);

data_new.self_absco_ref = cav_self296;
data_new.self_texp = log(cav_self296./cav_self260) / log(260/296);
data_new.for_absco_ref = data_new.for_absco_ref(3:end).*radiation_term;
data_new.use_radiation_term = int8(0);

attr_new = attr;
attr_new.global.title='CAVIAR self and MT_CKD-4.1.1 foreign water vapour continua';
attr_new.self_absco_ref.units = 'cm**2/molecule';
attr_new.for_absco_ref.units = 'cm**2/molecule';
attr_new.use_radiation_term.dimensions = {};
attr_new.use_radiation_term.long_name = 'Should the radiation term be used?';
attr_new.global.contact = 'Robin Hogan (r.j.hogan@ecmwf.int)';
attr_new.global.copyright = ['Self continuum: University of Reading' 10 ...
			     'Foreign continuum: ' attr.global.copyright];
attr_new.global.notes = 'The absorption coefficients do not need to be multiplied by the radiation term.';
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
