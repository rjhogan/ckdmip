

%ver = {'fsck-tol0.04','fsck-tol0.01','wide-tol0.04','wide-tol0.01','narrow-tol0.04','narrow-tol0.01'};
PLOTDIR='/home/rd/parr/src/ckdmip-plots';
bs = {'fsck','wide','narrow'};
application='climate';
%application='global-nwp'
%application='limited-area-nwp'
if strcmp(application, 'climate')
  ng = [9 14 21 27 40 59;
	14 19 30 43 57 87;
	22 28 37 48 69 101]';
elseif strcmp(application, 'global-nwp');
  ng = [11 16 22 30 44 75;
       11 16 28 38 49 74;
       20 27 36 48 69 93]';
else
  ng = [8 12 16 23 35 55;
	10 15 23 31 42 61;
	17 25 28 39 58 77]';
end

bs = 'rgb';ng = 29;
bs = {'fsck'};ng=[24 28 32]';
bs = {'fsck'};ng=[16]';
%bs = {'narrow'};ng = 64;
%bs = {'narrow'};ng=37;
%bs = {'narrow'};ng=64;
bs = {'fsck'};ng=[8 12 16 20 24 28 32 36 40 48 64]';
bs = {'narrow'};ng=64;
bs = {'fsck'};ng=[48]';

if 1
for ibs = 1:length(bs)
  gpt = ng(:,ibs)
  for ig = 1:length(gpt)
    v = [bs{ibs} '-' num2str(gpt(ig)) ''];
    suf = {'b'};
    %suf = {''};
    for is = 1:length(suf)
      s = suf{is};
      evaluate_ckd('ecckd-1.4',['ecCKD-' v],'lw',application,[v s],'png',[v s]);
    end
  end
end
end

if 0
  %figure(10)
  %metric= plot_accuracy_efficiency('ecckd-0.5','ecCKD','lw',application,bs,ng)
  figure(14)
  clf
  metric= plot_accuracy_efficiency_one_row('ecckd-1.2','ecCKD','lw',application,bs([1 1]),ng(:,[1 1]),{'evaluation1','evaluation2'})
%print_pdf([PLOTDIR '/ecckd-0.5_evaluation1_lw_' application '_accuracy-efficiency.pdf']);
%print_png([PLOTDIR '/ecckd-0.5_evaluation1_lw_' application '_accuracy-efficiency.png'],'150');
end
