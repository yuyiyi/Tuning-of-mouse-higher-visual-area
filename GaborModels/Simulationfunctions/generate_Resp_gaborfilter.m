function [resp, untuneGain] = generate_Resp_gaborfilter(stim, gparams, simpara, simplecell, complexcell, dirSelectiveCell, withLocalNormalization)
resp = [];
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
[SimpleUnit0, SimpleUnit90] = make3dgabor_V2([filtsize, filtsize, t_intgrat], gparams);
if simplecell
    resp_simple = sameconv(stim, reshape(SimpleUnit0, filtsize.^2, [])');
    resp = cat(2, resp, resp_simple); 
end
if complexcell
    resp_simple90 = sameconv(stim, reshape(SimpleUnit90, filtsize.^2, [])');
    resp_complex = sqrt(max(0,resp_simple).^2+max(0,resp_simple90).^2)/2; %motion energy model
    resp = cat(2, resp, resp_complex); 
end
if dirSelectiveCell
    gparams(3) = mod(gparams(3)+180,360);
    [SimpleUnit180, SimpleUnit270] = make3dgabor_V2([filtsize, filtsize, t_intgrat], gparams);
    resp_simple180 = sameconv(stim, reshape(SimpleUnit180, filtsize.^2, [])');    
    resp_simple270 = sameconv(stim, reshape(SimpleUnit270, filtsize.^2, [])');
    resp_complex180 = sqrt(max(0,resp_simple180).^2+max(0,resp_simple270).^2)/2; %motion energy model
    resp_DirSel = resp_complex-resp_complex180;
    resp = cat(2, resp, resp_DirSel); 
end
if nargin>6
    resp_untuned = untunedcomponent(stim, simpara, filtsize, gparams);
    untuneGain = rand(1,3)*0.1+0.4;
    if withLocalNormalization && simplecell
        simple_untune = simpara.nlfun(double(resp_simple))-simpara.nlfun(double(resp_untuned))*untuneGain(1); 
        resp = cat(2, resp, simple_untune); 
    end
    if withLocalNormalization && complexcell
        complex_untune = simpara.nlfun(double(resp_complex))-simpara.nlfun(double(resp_untuned))*untuneGain(2);
        resp = cat(2, resp, complex_untune); 
    end
    if withLocalNormalization && dirSelectiveCell    
        DirSel_untune = simpara.nlfun(double(resp_DirSel))-simpara.nlfun(double(resp_untuned))*untuneGain(3);
        resp = cat(2, resp, DirSel_untune); 
    end
end
