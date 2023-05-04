% Convert CAVIAR continuum observations to a Fortran statement to be
% added to caviar_continuum.F90
data = load('SOC_new_296K.dat');
wn=data(:,1);
self=data(:,2);
% Add zero at the end 
wn(end+1) = wn(end)+10;
self(end+1)=0;
nwav = length(wn);

dataf1 = load('Fore_296+310K_1300-1800cm-1.dat');
dataf2 = load('Fore_400K_1800-8970cm-1.dat');

foreign = zeros(size(self));

% Merge two sources of foreign continuum
index = find(wn < 1800);
foreign(index) = interp1(dataf1(:,1),dataf1(:,2),wn(index));
index = find(wn >= 1800);
foreign(index) = interp1(dataf2(:,1),dataf2(:,2),wn(index)).*400./296; % Convert to 296 K
foreign(find(~isfinite(foreign))) = 0;

fid = fopen('caviar.F90','w');
fprintf(fid,'    integer,    parameter :: nwav_lut = %d\n', nwav);
fprintf(fid,'    real(jprb), parameter :: temperature_lut = 296.0_jprb\n');
fprintf(fid,'    real(jprb), parameter :: wavenumber_lut_cm1(nwav_lut) = [ &\n');
str=sprintf('      &  %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, &\n',wn);
str = strrep(str,' 0_jprb',' 0.0_jprb');
fprintf(fid,str(1:end-2));
fprintf(fid,' ]\n');
fprintf(fid,'    real(jprb), parameter :: self_lut(nwav_lut) = [ &\n');
str=sprintf('      &  %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, &\n',self);
str = strrep(str,' 0_jprb',' 0.0_jprb');
fprintf(fid,str(1:end-2));
fprintf(fid,' ]\n');
fprintf(fid,'    real(jprb), parameter :: foreign_lut(nwav_lut) = [ &\n');
str=sprintf('      &  %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, %g_jprb, &\n',foreign);
str = strrep(str,'0_jprb','0.0_jprb');
fprintf(fid,str(1:end-2));
fprintf(fid,' ]\n');
fclose(fid)
