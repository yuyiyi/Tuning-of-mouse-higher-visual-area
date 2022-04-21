function [resp, untuneGain] = generate_Resp_Gaussfilter(stim, gparams, simpara, simplecell, complexcell, dirSelectiveCell, withLocalNormalization)
resp = [];
filtsize = simpara.filtsize;
t_intgrat = simpara.t_intgrat;
[Gauss_pos, Gauss_neg, ~, gauss_surround] = make_GaussKernel([filtsize, filtsize, t_intgrat], gparams);
if gparams(9)>1
    if simplecell
        resp_pos = sameconv(stim, reshape(Gauss_pos, filtsize.^2, [])');
        resp = cat(2, resp, resp_pos); 
        resp_neg = sameconv(stim, reshape(Gauss_neg, filtsize.^2, [])');
        resp = cat(2, resp, resp_neg); 
    end
    if complexcell
        resp_complex = sqrt(max(0,resp_pos).^2+max(0,resp_neg).^2)/2; %motion energy model
        resp = cat(2, resp, resp_complex); 
    end    
    if dirSelectiveCell
        gparams(3) = mod(gparams(3)+90,360);
        [Gauss_pos90, Gauss_neg90, ~, gauss_surround90] = make_GaussKernel([filtsize, filtsize, t_intgrat], gparams);
        resp_pos90 = sameconv(stim, reshape(Gauss_pos90, filtsize.^2, [])');
        %     resp = cat(2, resp, resp_pos90); 
        resp_neg90 = sameconv(stim, reshape(Gauss_neg90, filtsize.^2, [])');
        %     resp = cat(2, resp, resp_neg90); 
        resp_complex90 = sqrt(max(0,resp_pos90).^2+max(0,resp_neg90).^2)/2; %motion energy model
        %     resp = cat(2, resp, resp_complex90); 
        resp_DirSel = resp_complex-resp_complex90;
        resp = cat(2, resp, resp_DirSel); 
    end    
    if nargin>6
        untuneGain = rand(1)*0.1+0.4;
        resp_surround = [];
        for k = 1:length(gauss_surround)
            resp_surround(:,k) = untuneGain * sameconv(stim, reshape(gauss_surround{k}, filtsize.^2, [])');
        end
        resp_untuned = max(0,double(resp_surround(:,1)))+max(0,double(resp_surround(:,2)));
        % simple/complex cell - on+off untuned supression
        resp_suppress = max(0,double(resp(:,1:3)))-repmat(resp_untuned,1,3);
        resp = cat(2, resp, resp_suppress); 
        
%         resp_surround = [];
%         for k = 1:length(gauss_surround)
%             resp_surround(:,k) = untuneGain * sameconv(stim, reshape(gauss_surround90{k}, filtsize.^2, [])');
%         end
%         resp_untuned = max(0,double(resp_surround(:,1)))+max(0,double(resp_surround(:,2)));
%         resp_suppress = max(0,double(resp(:,4)))- resp_untuned;
%         resp = cat(2, resp, resp_suppress); 
    else
        untuneGain = [];
    end
    
elseif gparams(9) == 1
    if simplecell
        resp_pos = sameconv(stim, reshape(Gauss_pos, filtsize.^2, [])');
        resp = cat(2, resp, resp_pos); 
        resp_neg = sameconv(stim, reshape(Gauss_neg, filtsize.^2, [])');
        resp = cat(2, resp, resp_neg); 
    end
    if nargin>6
        untuneGain = rand(1)*0.1+0.9;
        resp_surround = [];
        for k = 1:length(gauss_surround)
            resp_surround(:,k) = untuneGain*sameconv(stim, reshape(gauss_surround{k}, filtsize.^2, [])');
            resp_untuned = max(0,double(resp_surround(:,k)));
            resp_suppress = max(0,double(resp(:,k)))-resp_untuned;
            resp = cat(2, resp, resp_suppress); 
        end
    else
        untuneGain = [];
    end
end