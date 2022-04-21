function [gauss_combine, gauss_components] = make_AngleKernel(xytsize, params)
%%%% generate kernels with two gaussian components with certain angle

cx = params(1);
cy = params(2);
ct = 0.5;
dirlist(1) = mod(params(3), 360); % orientation of the first components
ang = params(4); % angle between two gaussian kernels
dirlist(2) = mod(dirlist(1)+ang, 360); % orientation of the second components
senvlist(1) = params(5); % spatial width of the first gaussian kernel
SpatialARlist(1) = params(6); % aseptic ratio of the first gaussian kernel
senvlist(2) = params(7); % spatial width of the second gaussian kernel
SpatialARlist(2) = params(8); % aseptic ratio of the second gaussian kernel
tenv = params(9); % temporal width of the first gaussian kernel

xx = 0:(1/(xytsize(1)-1)):1;
yy = 0:(1/(xytsize(2)-1)):1;
if length(xytsize) >= 3 & xytsize(3)>1
	dt = 0:(1/(xytsize(3)-1)):1;
else  % only 1 time slice
	xytsize(3) = 1;
	dt = 0.5;
end
[iy, ix, it] = ndgrid(xx, yy, dt);

ii = [1 1; -1 -1; 1 -1; -1 1];
clear gauss_components
for k = 1:size(ii,1)
    gauss = zeros(size(ix));
    for i = 1:2
        dir = dirlist(i);
        senv = senvlist(i);        
        SpatialAR = SpatialARlist(i);
        if SpatialAR>1
            dx = -senv*SpatialAR*sin(deg2rad(dir));
            dy = senv*SpatialAR*cos(deg2rad(dir));
        else
            dx = senv*cos(deg2rad(dir));
            dy = -senv*sin(deg2rad(dir));            
        end
        cx1 = cx+dx*sqrt(2*log(2)); 
        cy1 = cy+dy*sqrt(2*log(2));
        a = cos(deg2rad(dir))^2/(2*senv^2) + sin(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        b = -sin(2*deg2rad(dir))/(4*senv^2) + sin(2*deg2rad(dir))/(4*(senv*SpatialAR)^2);
        c = sin(deg2rad(dir))^2/(2*senv^2) + cos(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        gausstmp = ii(k,i)*exp(-(a*(ix-cx1).^2 + 2*b*(ix-cx1).*(iy-cy1) + c*(iy-cy1).^2)-(it-ct).^2/(2*tenv^2));
        gauss_components{k,i} = gausstmp;
        gauss = gauss+gausstmp;
    end
    gauss_combine{k} = gauss;
end


