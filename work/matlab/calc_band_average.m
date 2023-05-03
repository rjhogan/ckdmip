function out = calc_band_average(wn1, wn2, specdef, in)
nband = length(wn1);
[nspec,nz,ncol] = size(in);
out = zeros(nband,nz,ncol);
for iband = 1:nband
  index = find(specdef.wavenumber1 >= wn1(iband) & specdef.wavenumber2 <= wn2(iband));
  weight = sum(specdef.gpoint_fraction(index,:),1);
  for icol = 1:ncol
    out(iband,:,icol) = in(:,:,icol)'*weight(:);
  end
end
