function r = untunedcomponent_v2(stimtmp, simpara, gparams)
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
dx = 0:(1/(filtsize-1)):1;
dy = 0:(1/(filtsize-1)):1;
cx = gparams(1);
cy = gparams(2);
dir = gparams(3);
SpatialAR = gparams(8);
if t_intgrat > 1
	dt = 0:(1/(t_intgrat-1)):1;
else  % only 1 time slice
	dt = 0.5;
end
[iy, ix, it] = ndgrid(dx, dy, dt);

senv = gparams(6)*1.5;
dt = 0:(1/(t_intgrat-1)):1;
tenv = gparams(7);

if length(gparams)<9
    ii = sign(rand(1)-0.5);
    if ii==0
        ii=1;
    end
    gparams(9) = ii;
elseif abs(gparams(9))~=1
    ii = sign(rand(1)-0.5);
    if ii==0
        ii=1;
    end
    gparams(9) = ii;    
else
    ii = gparams(9);
end

a = cos(deg2rad(dir))^2/(2*senv^2) + sin(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
b = -sin(2*deg2rad(dir))/(4*senv^2) + sin(2*deg2rad(dir))/(4*(senv*SpatialAR)^2);
c = sin(deg2rad(dir))^2/(2*senv^2) + cos(deg2rad(dir))^2/(2*(senv*SpatialAR)^2);
gauss = ii*exp(-(a*(ix-cx).^2 + 2*b*(ix-cx).*(iy-cy) + c*(iy-cy).^2)-(it-0.5).^2/(2*tenv^2));

% gauss = ii*exp( - ((ix-cx).^2+(iy-cy).^2)/(2*senv^2) - (it-0.5).^2/(2*tenv^2)  );
r = sameconv(stimtmp, reshape(gauss, [], size(gauss,3))');



