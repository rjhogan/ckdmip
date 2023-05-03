PLOTDIR='/home/rd/parr/src/ckdmip-plots';
bs = {'double','wide','narrow'};
application='climate';
%application='global-nwp'
%application='limited-area-nwp'
if strcmp(application, 'climate')
  ng = [9 14 21 27 40 59;
	11 19 31 52 74 102;
	15 20 29 42 57 81]';
  ng = [9 14 21 00 27 40 59;
	11 19 31 38 52 74 102;
	15 20 29 32 42 57 81]';
  ng = [9 14 21 00 27 40 59;
	11 14 22 27 38 60 84
        14 16 23 25 35 47 66]';
  ng = [0 0 0 0 0 0 0;
	11 14 22 27 38 60 84
        14 16 23 25 35 47 66]';
elseif strcmp(application, 'global-nwp');
  ng = [24 0  0  0  0 0;
       11 19 29 50 70 97;
       15 20 29 42 56 80]';
  ng = [27 45  72 111 163 244;
	8 15 21 37 58 81;
	14 16 23 35 47 65]';
  ng = [0 0 0 0 0 0 0;
	11 15 21 27 37 58 81;
	14 16 23 25 35 47 65]';
else
  ng = [8 12 16 23 35 55;
	10 16 28 45 66 91;
	14 18 26 37 52 74]';
  ng = [0 0 0 0 0 0 0;
	10 11 19 25 34 52 75;
        13 14 20 24 30 41 61]';
end

%bs = {'rgb'};ng = [19,27,38]';
bs = {'rgb'};ng = [16 17 20 21 24]';
bs = {'wide'};ng = [37]';
%bs = {'gb'};ng = [32]';
bs = {'rgb'};ng = [16 24 32]';
bs = {'rgb'};ng = [64]';
bs = {'window'};ng = [48]';
%bs = {'double'};ng = [18]';
%bs = {'fine'};ng=65;
bs = {'narrow'};ng=64;
%bs = {'rgb'};ng=[8 12 16 20 24 28 32 36 40 48 64]';
%bs = {'vfine'};ng=96;
%bs ={'rgb'};ng=32;
%bs = {'double'};ng = 8';

%bs = {'double'};
%ng = [45 76 115]';

%    suf = {'-scaled','-raw2'};
    suf = {'-scaled','-raw','-raw2'};
    suf = {'-raw',''};
%    suf = {'-scaled'};
    suf = {''};

if 1
for ibs = 1:length(bs)
  gpt = ng(:,ibs)
  for ig = 1:length(gpt)
    %v = [bs{ibs} '-' num2str(gpt(ig)) 'b'];
    v = [bs{ibs} '-' num2str(gpt(ig)) ''];
    for is = 1:length(suf)
      s = suf{is};
      evaluate_ckd('ecckd-1.2',['ecCKD-' upper(v)],'sw',application,[v s],'png',[v]);
    end
  end
end
end
%evaluate_ckd('ecckd-0.5','ecCKD-double-32','sw','global-nwp','double-32','png','double-32');

if 0
%  figure(13)
%  metric= plot_accuracy_efficiency('ecckd-1.2','ecCKD','sw',application,bs,ng)
  figure(14)
  metric= plot_accuracy_efficiency_one_row('ecckd-1.2','ecCKD','sw',application,bs([1 1]),ng(:,[1 1]),{'evaluation1','evaluation2'})
  %print_pdf([PLOTDIR '/ecckd-0.6_evaluation1_sw_' application '_accuracy-efficiency.pdf']);
  %print_png([PLOTDIR '/ecckd-0.6_evaluation1_sw_' application '_accuracy-efficiency.png'],'150');
end

if 0
od=loadnc('/hugetmp/parr/ckdmip/results/ecckd-0.5/sw_optical-depth/ecckd-0.5_evaluation1_sw_global-nwp_narrow-27-raw_optical-depth_present.nc');
raw=loadnc('/hugetmp/parr/ckdmip/results/ecckd-0.5/sw_fluxes/ecckd-0.5_evaluation1_sw_global-nwp_narrow-27-raw_fluxes_present.nc');
ckd=loadnc('/hugetmp/parr/ckdmip/results/ecckd-0.5/sw_fluxes/ecckd-0.5_evaluation1_sw_global-nwp_narrow-27_fluxes_present.nc');
lbl=loadnc('/hugetmp/parr/ckdmip/evaluation1/sw_fluxes/ckdmip_evaluation1_sw_fluxes_present.h5');
figure(100)
clf
plot(squeeze(lbl.flux_dn_direct_sw(end,3,:)), ...
     squeeze(ckd.flux_dn_direct_sw(end,3,:)-lbl.flux_dn_direct_sw(end,3,:)),'ko');
hold on
plot(squeeze(lbl.flux_dn_direct_sw(end,3,:)), ...
     squeeze(raw.flux_dn_direct_sw(end,3,:)-lbl.flux_dn_direct_sw(end,3,:)),'ro');
plot(squeeze(lbl.flux_dn_direct_sw(end,3,:)), ...
     squeeze(od.flux_dn_direct_sw(end,:))'-squeeze(lbl.flux_dn_direct_sw(end,3,:)),'bx');
end
