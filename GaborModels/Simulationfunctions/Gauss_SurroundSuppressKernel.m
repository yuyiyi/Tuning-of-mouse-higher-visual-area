function gauss_surround = Gauss_SurroundSuppressKernel(xytsize, cx, cy, senv_x, senv_y, tenv, dir)
xx = 0:(1/(xytsize(1)-1)):1;
yy = 0:(1/(xytsize(2)-1)):1;
if length(xytsize) >= 3 & xytsize(3)>1
	dt = 0:(1/(xytsize(3)-1)):1;
else  % only 1 time slice
	xytsize(3) = 1;
	dt = 0.5;
end
[iy, ix, it] = ndgrid(xx, yy, dt);

senv_x = senv_x*1.5;
senv_y = senv_y*1.5;
tenv = tenv*1;

a = cos(deg2rad(dir))^2/(2*senv_x^2) + sin(deg2rad(dir))^2/(2*senv_y^2);
b = -sin(2*deg2rad(dir))/(4*senv_x^2) + sin(2*deg2rad(dir))/(4*senv_y^2);
c = sin(deg2rad(dir))^2/(2*senv_x^2) + cos(deg2rad(dir))^2/(2*senv_y^2);
gausstmp = exp(-(a*(ix-cx).^2 + 2*b*(ix-cx).*(iy-cy) + c*(iy-cy).^2)-(it-0.5).^2/(2*tenv^2));
% norm = sqrt(sum(sum(sum(gausstmp.^2),2),3));

gauss_surround{1} = gausstmp;
gauss_surround{2} = -gausstmp;
% gauss_surround{1} = exp( - ((ix-cx).^2)/(2*senv_x^2) - ((iy-cy).^2)/(2*senv_y^2) - (it-0.5).^2/(2*tenv^2));
% gauss_surround{2} = -exp( - ((ix-cx).^2)/(2*senv_x^2) - ((iy-cy).^2)/(2*senv_y^2) - (it-0.5).^2/(2*tenv^2));