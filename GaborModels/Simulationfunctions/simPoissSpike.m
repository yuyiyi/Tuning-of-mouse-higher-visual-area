function [spsum, meanFR] = simPoissSpike(resp,simpara)
dtStim = simpara.dtStim;
Simdt = simpara.Simdt;
SimRep = simpara.SimRep;
binWidth = simpara.binWidth;
% b = simpara.bias;
filtsize = simpara.filtsize;
t = simpara.t_intgrat;
if ~isfield(simpara, 'nlfun')
    simpara.nlfun = @exp;
end
nlfun = simpara.nlfun;

resp = repmat(resp', round(dtStim/Simdt), 1);
b = max(resp(:))*0.8;
resp = resp(:)-b;
[spiketime, simSpiketrain] = simuByinputV2(resp, Simdt, SimRep, nlfun);
% bin simulated spike train at 100 ms
simSpiketrain = squeeze(sum(reshape(simSpiketrain', round(dtStim/Simdt)*(binWidth/dtStim), [],SimRep)));
spsum = sum(simSpiketrain(:))/SimRep;
meanFR = mean(simSpiketrain,2);
% a = sum(reshape(simSpiketrain,1/dtStim,[],SimRep));
% a = squeeze(a);
% meanFR = mean(a,2);

