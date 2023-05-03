PLOTDIR='/home/rd/parr/src/ckdmip-plots';
bs = 'wide';
applications = {'climate','global-nwp','limited-area-nwp'};
ng = [9 14 21 27 40 59;
      11 16 22 30 44 75;
      8 12 16 23 35 55]';
figure(10)
clf
metric= plot_accuracy_efficiency_applications('ecckd-0.5','ecCKD','lw',applications,bs,ng,50)
print_png([PLOTDIR '/ecckd-0.5_evaluation1_lw_' bs '_accuracy-efficiency_present.png'],'150');

bs = 'wide';
ng = [11 14 22 27 38 60 84;
      11 15 21 27 37 58 81;
      10 11 19 25 34 52 75]';
figure(11)
clf
metric= plot_accuracy_efficiency_applications('ecckd-0.6','ecCKD','sw',applications,bs,ng(2:end,:),70)
print_png([PLOTDIR '/ecckd-0.6_evaluation1_sw_' bs '_accuracy-efficiency_present.png'],'150');
