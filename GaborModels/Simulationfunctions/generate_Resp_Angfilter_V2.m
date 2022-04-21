function resp = generate_Resp_Angfilter_V2(stim, gparams, simpara)
resp = [];
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
[gabor_combine, gauss_components] = make_AngleKernel_v2([filtsize,filtsize,t_intgrat], simpara.monitorresolution, gparams);
for i = 1:length(gabor_combine)
    h1 = gauss_components{i};
    resptmp1 = max(0, sameconv(stim, reshape(h1, filtsize.^2, [])'));
    
    h = gabor_combine{i};
    resptmp = sameconv(stim, reshape(h, filtsize.^2, [])');    
    resp = cat(2, resp, resptmp); 
    resptmp = resptmp - (rand(1)*0.1+0.9)*(resptmp1);
    resp = cat(2, resp, resptmp); 
end