function [Gauss_pos, Gauss_neg, gauss_components, gauss_surround] = make_GaussKernel(xytsize, params)

% function [gabor, gabor90] = make3dgabor(xytsize, params)
% 
%   returns a gabor functions of size X-by-Y-by-T, specified by a vector PARAMS.
%
% INPUT:
%     [xytsize] = vector of x, y, and t size, i.e. [64 64 5]
% [params(1:2)] = center_x, center_y
%                 The spatial center of the kernel. The axes are normalized
%                 to 0 (lower left corner) to 1(upper right corner).
%                 e.g., [0.5 0.5] put the kernel at the center of the matrix.
%   [params(3)] = The direction of the Gabor function in degree (0-360).
% [params(4:5)] = Spatial frequency and temporal frequency
%                 They determine how many cycles in XYTSIZE pixels for each dimension.
% [params(6:8)] = Spatial and Temporal envelope size in standard deviation
%   [params(9)] = number of gaussian components
% OUTPUT:
%       [gabor] = a gabor function of size X-by-Y-by-T, specified by a vector PARAMS.
%     [gabor90] = the quadrature pair Gabor function
%

cx = params(1);
cy = params(2);
dir = params(3); % spatial angle of gaussian kernels
sf = params(4); % defines the spatial distance between individual gaussian components
tf = params(5); % defines the time distance bewteen individual gaussian components
senv = params(6); 
SpatialAR = params(7);
tenv = params(8);
Nkernel = params(9); % number of Gaussian components

xx = 0:(1/(xytsize(1)-1)):1;
yy = 0:(1/(xytsize(2)-1)):1;
if length(xytsize) >= 3 & xytsize(3)>1
	dt = 0:(1/(xytsize(3)-1)):1;
else  % only 1 time slice
	xytsize(3) = 1;
	dt = 0.5;
end
[iy, ix, it] = ndgrid(xx, yy, dt);

% spatial frequency
spatial_threshold = 1/xytsize(1)*Nkernel;
dd = (1/sf)/2;
if dd > max([senv, senv*SpatialAR])*2
    dd = max([senv, senv*SpatialAR])*2;    
end
dx = dd*cos(deg2rad(-dir));
dy = dd*sin(deg2rad(-dir));
    
    
% % wx = sf*cos((-dir)/180*pi); % cycle per filter window size, counter-clockwise 
% if abs(wx)>spatial_threshold && Nkernel>1
%     % distance bewteen two peak (half cycle)
%     dx = (1/wx)/2; % distance bewteen individual gaussian components
% else
%     dx = 0;
% end
% wy = sf;
% % wy = sf*sin((-dir)/180*pi);
% if abs(wy)>spatial_threshold && Nkernel >1
%     dy = (1/wy)/2; % distance bewteen individual gaussian components
% else
%     dy = 0;
% end
temporal_threshold = 1/xytsize(3)*Nkernel;
wt = tf; % cycle per integration window size
if wt>temporal_threshold && Nkernel>1
    dz = (1/wt)/2; % distance bewteen individual gaussian components
else
    dz = 0;
end

if Nkernel==1
    ii = [1; -1]; % on, off kernel
elseif Nkernel==2
    ii = [1 -1; -1 1];
elseif Nkernel==3
    ii = [1 -1 1; -1 1 -1];
end

cx0 = cx-(dx*(Nkernel-1))/2;
cy0 = cy-(dy*(Nkernel-1))/2;
ct0 = 0.5-(dz*(Nkernel-1))/2;
clear gauss_components

for k = 1:size(ii,1)
    gauss = zeros(size(ix));
    for ni = 1:Nkernel
        cx1 = cx0+dx*(ni-1);
        cy1 = cy0+dy*(ni-1);
        ct1 = ct0+dz*(ni-1);        
        a = cos(deg2rad(dir))^2/(2*senv^2) + sin(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        b = -sin(2*deg2rad(dir))/(4*senv^2) + sin(2*deg2rad(dir))/(4*(senv*SpatialAR)^2);
        c = sin(deg2rad(dir))^2/(2*senv^2) + cos(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
        gausstmp = ii(k,ni)*exp(-(a*(ix-cx1).^2 + 2*b*(ix-cx1).*(iy-cy1) + c*(iy-cy1).^2)-(it-ct1).^2/(2*tenv^2));
        % in the general gaussian function, the direction is counter clockwise
        gauss = gauss + gausstmp;
        gauss_components{k,ni} = gausstmp;
    end
    if k==1
        Gauss_pos = gauss;
%         norm = sqrt(sum(sum(sum(gauss.^2),2),3));
%         Gauss_pos = gauss/norm;
    else
        Gauss_neg = gauss;
%         norm = sqrt(sum(sum(sum(gauss.^2),2),3));
%         Gauss_neg = gauss/norm;
    end
end

senv_x_surround = abs(dy)*(Nkernel-1) + senv;
senv_y_surround = abs(dx)*(Nkernel-1) + senv*SpatialAR;
tenv_surround = tenv + abs(dz)*(Nkernel-1);
% if senv_x_surround/senv_y_surround < 1
%     dir = dir-90;
% end
gauss_surround = Gauss_SurroundSuppressKernel(xytsize, cx, cy, senv_x_surround, senv_y_surround, tenv_surround, dir);
