function [simSpiketrain_bin, simSpiketrain, rb, S_lf, simRate, spiketime] = simuByinputV3_batch(resp, simpara, T, sps_bin, display)
if nargin < 5
    display = 0;
end
clear spiketime
simSpiketrain = []; simSpiketrain_bin = []; rb = []; S_lf = []; simRate = [];
for i = 1:size(resp, 2)
    resptmp = interp1([1:size(resp,1)]*simpara.dtStim, resp(:,i)', [1:T/simpara.Simdt]*simpara.Simdt);
    resptmp(isnan(resptmp)) = 0;
    [sptime, simSpiketrain(:,:,i), simRate(:,i)] = simuByinputV3(resptmp', simpara.Simdt, simpara.SimRep, simpara.nlfun, simpara.ih);
    sps = reshape(sum(reshape(simSpiketrain(:,:,i), simpara.SimRep, sps_bin,[]),2), simpara.SimRep, []);
    
    % reliability
    rb1 = corr(sps');
    A = triu(ones(size(rb1)),1);
    rb(i) = nanmean(rb1(A==1));
    % lifetime sparseness
    r = nanmean(sps,1)';
    N = size(r,1);
    S_lf(i) = (1-1/N*((sum(r)).^2./sum(r.^2)))./(1-1/N);

    simSpiketrain_bin(:,:,i) = sps;
    
    % rester plots
    if display
        subplot(size(resp,2),1,i)
        createRaster3(sptime', 0, T);
    end
    spiketime{i} = sptime;
end
