DIR='/hugetmp/parr/ckdmip/results';
PLOTDIR='/home/rd/parr/src/ckdmip-plots';
TOOL='ecckd-0.5';
%TOOL='ecrad-rrtmg';
APP='climate';
APP='global-nwp';
%APP='limited-area-nwp';

if strcmp(TOOL(1:5), 'ecckd')
  if strcmp(APP, 'climate')
    ngtable = [9 14 21 27 40 59;
	       14 19 30 43 57 87;
	       22 28 37 48 69 101]';
  elseif strcmp(APP, 'global-nwp');
    ngtable = [11 16 22 30 44 75;
	       11 16 28 38 49 74;
	       20 27 36 48 69 93]';
  else
    ngtable = [8 12 16 23 35 55;
	       10 15 23 31 42 61;
	       17 25 28 39 58 77]';
  end
  bs = {'fsck','wide','narrow'};
  cmax = [0.2 0.3 0.4];
else
  ngtable=140;
  bs = {'narrow'};
  cmax=0.4;
end

bb_wide = [500, 820, 1180, 1800];
bb_narrow = [0, 350, 630, 700, 980, 1080, 1390, 1480, 2080];

for ibs = 1;%1:length(bs)
  for ig = 4;%1:size(ngtable,1)
    ng = ngtable(ig,ibs);
    file=[DIR '/' TOOL '/lw_spectral-definition/' TOOL '_lw_' APP '_' bs{ibs} '-' num2str(ng) '_spectral-definition.nc'];
    d=loadnc(file);

    clf
    set(gcf,'paperposition',[0.5 0.5 25 10]);
    axes('position',[0.15 0.2 0.65 0.7]) 
    data = d.gpoint_fraction([1:end end],[1:end end])';
    data(find(data <= 0)) = NaN;
    pcolor([d.wavenumber1;d.wavenumber2(end)], 0.5+[0:ng], data);
    shading flat
    chiljet
    caxis([0 cmax(ibs)]);
%    cm = colormap;
%    cm(1,:) = [1 1 1];
%    colormap(cm);
    set(gca,'layer','top');
    hold on
    plot([bb_wide;bb_wide],[0.5 ng+0.5]'*ones(1,length(bb_wide)),'k-','linewidth',1)
    plot([bb_narrow;bb_narrow],[0.5 ng+0.5]'*ones(1,length(bb_narrow)),'k:')
    xlabel('Wavenumber (cm^{-1})');
    ylabel('{\itk} term');

    %h=axes('position',[0.9 0.2 0.025 0.3]);
    h=colorbar('location','eastoutside');
    %set(h,'position',[0.77 0.20 0.02 0.5]);
    pp=get(h,'position');
    set(h,'position',[0.83 0.3 0.025 0.5]);%pp+[0.15 0.1 0 -0.2])
    h.Label.String = 'Fraction of {\itk} term in 10 cm^{-1} interval';
    h.FontSize=10;
    h.Label.FontSize=10;

    drawnow
    disp([PLOTDIR '/' TOOL '_lw_' APP '_' bs{ibs} '-' num2str(ng) '_spectral-definition.pdf'])
%    print_pdf([PLOTDIR '/' TOOL '_lw_' APP '_' bs{ibs} '-' num2str(ng) '_spectral-definition.pdf']);
    print_png([PLOTDIR '/' TOOL '_lw_' APP '_' bs{ibs} '-' num2str(ng) '_spectral-definition.png'],'150');
  end
end
