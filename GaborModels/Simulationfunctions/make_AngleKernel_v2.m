function [gabor_combine, gauss_components] = make_AngleKernel_v2(xytsize, monitorresolution, params)
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

sf1 = 1/(params(5)*(2*sqrt(2*log(2))*xytsize(1)))/monitorresolution;
sf = sf1*monitorresolution*xytsize(1); % cycle per bin

ft = 0; 

xx = 0:(1/(xytsize(1)-1)):1;
yy = 0:(1/(xytsize(2)-1)):1;
if length(xytsize) >= 3 & xytsize(3)>1
	dt = 0:(1/(xytsize(3)-1)):1;
else  % only 1 time slice
	xytsize(3) = 1;
	dt = 0.5;
end
[iy, ix, it] = ndgrid(xx, yy, dt);

k1 = 1;
ii = [1 1; -1 -1];
clear gauss_components gabor_combine
for k = 1:size(ii,1)
    gauss = zeros(size(ix));
    gauss_surround = gauss;
    gabor_comp = zeros(size(ix));
    gabor90_comp = zeros(size(ix));    
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
        
        fx = -sf*cos(dir/180*pi)*2*pi; % sf: frequency per bin
        fy = sf*sin(dir/180*pi)*2*pi;
        cx1 = cx+dx*sqrt(2*log(2)); 
        cy1 = cy+dy*sqrt(2*log(2));

        a = cos(deg2rad(dir))^2/(2*senv^2) + sin(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        b = -sin(2*deg2rad(dir))/(4*senv^2) + sin(2*deg2rad(dir))/(4*(senv*SpatialAR)^2);
        c = sin(deg2rad(dir))^2/(2*senv^2) + cos(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        gausstmp = ii(k,i)*exp(-(a*(ix-cx1).^2 + 2*b*(ix-cx1).*(iy-cy1) + c*(iy-cy1).^2)-(it-ct).^2/(2*tenv^2));

        grat = sin( (ix-cx)*fx + (iy-cy)*fy + (it-0.5)*ft);
        gabor_comp = gabor_comp + gausstmp.*grat;
        
        grat = cos( (ix-cx)*fx + (iy-cy)*fy + (it-0.5)*ft);
        gabor90_comp = gabor90_comp + gausstmp.*grat;
        gauss = gauss+gausstmp;
        
        
        a = cos(deg2rad(dir))^2/(2*(senv*1.5)^2) + sin(deg2rad(dir))^2/(2*(senv*1.5*SpatialAR)^2);
        b = -sin(2*deg2rad(dir))/(4*(senv*1.5)^2) + sin(2*deg2rad(dir))/(4*(senv*1.5*SpatialAR)^2);
        c = sin(deg2rad(dir))^2/(2*(senv*1.5)^2) + cos(deg2rad(dir))^2/(2*(senv*1.5*SpatialAR)^2);
        gausstmp = ii(k,i)*exp(-(a*(ix-cx1).^2 + 2*b*(ix-cx1).*(iy-cy1) + c*(iy-cy1).^2)-(it-ct).^2/(2*tenv^2));
        gauss_surround = gauss_surround+gausstmp;
    end
    
    gauss_components{k1} = gauss_surround;
    gabor_combine{k1} = gabor_comp;
    k1 = k1+1;
    gauss_components{k1} = gauss_surround;    
    gabor_combine{k1} = gabor90_comp;
    k1 = k1+1;
end
