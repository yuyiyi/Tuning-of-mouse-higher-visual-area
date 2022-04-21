function resp = generate_Resp_Angfilter(stim, gparams, simpara)
resp = [];
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
[gauss_combine, gauss_components] = make_AngleKernel([filtsize,filtsize,t_intgrat], gparams);
for i = 1:size(gauss_components,1)
    h1 = gauss_components{i,1};
    h2 = gauss_components{i,2};
    resptmp1 = max(0, sameconv(stim, reshape(h1, filtsize.^2, [])'));
    resptmp2 = max(0, sameconv(stim, reshape(h2, filtsize.^2, [])'));
    
    h = gauss_combine{i};
    resptmp = sameconv(stim, reshape(h, filtsize.^2, [])');    
    resp = cat(2, resp, resptmp); 
    resptmp = resptmp - (rand(1)*0.1+0.9)*(resptmp1+resptmp2);
    resp = cat(2, resp, resptmp); 
end