function [resp, spsum, meanFR] = simOpponentMotionCell(stimtmp, gparams, simpara)
dtStim = simpara.dtStim;
Simdt = simpara.Simdt;
binWidth = simpara.binWidth;
% b = simpara.bias;
filtsize = simpara.filtsize;
t = simpara.t_intgrat;
if ~isfield(simpara, 'nlfun')
    simpara.nlfun = @exp;
end
nlfun = simpara.nlfun;

[SimpleUnit0, SimpleUnit90] = make3dgabor_V2([filtsize, filtsize, t], gparams);
[resp1, ~, ~] = simMotionEnergyCell(stimtmp, SimpleUnit0, SimpleUnit90,simpara);

gparams(3) = mod(gparams(3)+180,360);
[SimpleUnit0, SimpleUnit90] = make3dgabor_V2([filtsize, filtsize, t], gparams);
[resp2, ~, ~] = simMotionEnergyCell(stimtmp, SimpleUnit0, SimpleUnit90,simpara);
resp = resp1-resp2;
r = resp;
r = repmat(r', round(dtStim/Simdt), 1);
if isfield(simpara, 'bias')
    b = simpara.bias;
else
    b = max(r(:))*0.75;
end
% b = max(r(:))*0.75;
r = r(:)-b;

if nargout>1
    SimRep = simpara.SimRep;
    [spiketime, simSpiketrain] = simuByinputV2(r, Simdt, SimRep, nlfun);
    % bin simulated spike train at 100 ms
    simSpiketrain = squeeze(sum(reshape(simSpiketrain', round(dtStim/Simdt)*(binWidth/dtStim), [],SimRep)));
    spsum = sum(simSpiketrain(:))/SimRep;
    meanFR = mean(simSpiketrain,2);
    % a = sum(reshape(simSpiketrain,1/dtStim,[],SimRep));
    % a = squeeze(a);
    % meanFR = mean(a,2);
end
