function r = untunedcomponent(stimtmp, simpara, filtsize_surround, gparams)
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
dx = 0:(1/(filtsize_surround-1)):1;
dy = 0:(1/(filtsize_surround-1)):1;
cx = gparams(1);
cy = gparams(2);
[iy, ix] = ndgrid(dx, dy);
senv1 = gparams(6)*filtsize/filtsize_surround;
senv2 = gparams(6)*filtsize/filtsize_surround*1.2;
G1 = 1/(2*pi*senv1)*exp( - 1/2*((ix-cx).^2+(iy-cy).^2)/(2*senv1^2));
G2 = 1/(2*pi*senv2)*exp( - 1/2*((ix-cx).^2+(iy-cy).^2)/(2*senv2^2));
dt = 0:(1/(t_intgrat-1)):1;
tenv = gparams(7);
G3 = exp( - (dt-0.5).^2/(2*tenv^2)  );

if length(gparams)<9
    ii = 1;
elseif abs(gparams(9))~=1
    ii = 1;
else
    ii = gparams(9);
end
h = ii*(G1-G2);
h = bsxfun(@times, repmat(h,1,1,length(G3)), reshape(G3,1,1,length(G3)));
norm = sqrt(sum(sum(sum(h.^2),2),3));
h = h/norm;
if length(size(stimtmp))==4
    r = bsxfun(@times, single(stimtmp),  single(h));
    r = squeeze(nansum(nansum(nansum(r,1),2),3));
elseif length(size(stimtmp))==2
    r = sameconv(stimtmp, reshape(h, [], size(h,3))');
end

