function [resp, resp_untuned] = generate_Resp_gaborfilter_v3(stim, gparams, simpara, stim_mask, display)
resp = [];
if nargin < 5
    display = 0;
end
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
[SimpleUnit0, SimpleUnit90, gparams, gauss] = make3dgabor_V3([filtsize, filtsize, t_intgrat], gparams);
SimpleUnit0 = SimpleUnit0/sqrt(sum(sum(sum(SimpleUnit0.^2))));
if display
    figure(1)
    subplot(1,3,1), imshow(SimpleUnit0(:,:,round(t_intgrat/2)),[])
    subplot(1,3,2), imshow(squeeze(SimpleUnit0(filtsize/2,:,:)),[])
    subplot(1,3,3), imshow(squeeze(SimpleUnit0(:,filtsize/2,:)),[])
end
SimpleUnit90 = SimpleUnit90/sqrt(sum(sum(sum(SimpleUnit90.^2))));
    resp_simple = sameconv(stim, reshape(SimpleUnit0, filtsize.^2, [])');
    resp = cat(2, resp, resp_simple); 
    resp_simple90 = sameconv(stim, reshape(SimpleUnit90, filtsize.^2, [])');
    resp_complex = sqrt(max(0,resp_simple).^2+max(0,resp_simple90).^2)/2; %motion energy model
    resp = cat(2, resp, resp_complex); 
    gparams(3) = mod(gparams(3)+180,360);
    [SimpleUnit180, SimpleUnit270] = make3dgabor_V2([filtsize, filtsize, t_intgrat], gparams);
    resp_simple180 = sameconv(stim, reshape(SimpleUnit180, filtsize.^2, [])');    
    resp_simple270 = sameconv(stim, reshape(SimpleUnit270, filtsize.^2, [])');
    resp_complex180 = sqrt(max(0,resp_simple180).^2+max(0,resp_simple270).^2)/2; %motion energy model
    resp_DirSel = resp_complex-resp_complex180;
    resp = cat(2, resp, resp_DirSel); 
    resp_untuned = sameconv(stim, reshape(gauss, [], size(gauss,3))');

