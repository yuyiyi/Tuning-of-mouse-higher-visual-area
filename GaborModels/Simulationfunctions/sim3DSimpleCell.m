function [resp, spsum, meanFR] = sim3DSimpleCell(stimtmp, SimpleUnit0, simpara)
dtStim = simpara.dtStim;
Simdt = simpara.Simdt;
binWidth = simpara.binWidth;
filtsize = simpara.filtsize;
t = simpara.t_intgrat;
if ~isfield(simpara, 'nlfun')
    simpara.nlfun = @exp;
end
nlfun = simpara.nlfun;

r = bsxfun(@times, single(stimtmp),  single(SimpleUnit0));
r = squeeze(nansum(nansum(nansum(r,1),2),3));
norm0 = sqrt(sum(sum(sum(SimpleUnit0.^2),2),3));
r = r/norm0;
% r = r*dtStim;
resp = r;
r = repmat(r', round(dtStim/Simdt), 1);
if isfield(simpara, 'bias')
    b = simpara.bias;
else
    b = max(r(:))*0.6;
end
r = r(:)-b;

if nargout>1
    SimRep = simpara.SimRep;
    [spiketime, simSpiketrain] = simuByinputV2(r, Simdt, SimRep,nlfun);
    % bin simulated spike train at 100 ms
    simSpiketrain = squeeze(sum(reshape(simSpiketrain', round(dtStim/Simdt)*(binWidth/dtStim), [],SimRep)));
    spsum = sum(simSpiketrain(:))/SimRep;
    meanFR = mean(simSpiketrain,2);
    % a = sum(reshape(simSpiketrain,1/dtStim,[],SimRep));
    % a = squeeze(a);
    % meanFR = mean(a,2);
end

