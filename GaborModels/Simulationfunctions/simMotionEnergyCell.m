function [resp, spsum, meanFR] = simMotionEnergyCell(stimtmp, SimpleUnit0, SimpleUnit90, simpara)
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

r0 = bsxfun(@times, single(stimtmp),  single(SimpleUnit0));
r0 = squeeze(nansum(nansum(nansum(r0,1),2),3));
norm0 = sqrt(sum(sum(sum(SimpleUnit0.^2),2),3));
r0 = r0/norm0;
r90 = bsxfun(@times, single(stimtmp),  single(SimpleUnit90));
r90 = squeeze(nansum(nansum(nansum(r90,1),2),3));
norm90 = sqrt(sum(sum(sum(SimpleUnit90.^2),2),3));
r90 = r90/norm90;
r = sqrt(max(0,r0).^2+max(0,r90).^2)/2; %motion energy model
% r = r*dtStim;
resp = r;
r = repmat(r', round(dtStim/Simdt), 1);
if isfield(simpara, 'bias')
    b = simpara.bias;
else
    b = max(r(:))*0.6;
end
% b = max(r(:))*0.85;
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

