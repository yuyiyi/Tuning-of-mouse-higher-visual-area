% GLM simulation 
function [spiketime, simSpiketrain] = simNeuCpl(resp, dt, rep, nlfun, ih, dc, dc0)
if nargin < 6
    dc = max(resp)*0.85;
end
if nargin<7
    dc0 = 0;
end
resp = bsxfun(@minus,resp,dc);
[rlen, ncells] = size(resp);
simSpiketrain = zeros(rep, rlen, ncells);
clear spiketime;
nbinsPerEval = 100;
hlen = size(ih,1);
% hlen = 0;
for perm = 1:rep
tsp(1,1:ncells) = {zeros(round(rlen/500),1)};  % allocate space for spike times
nsp = zeros(1,ncells);
jbin = 1;  % counter ?
tspnext = exprnd(1,1,ncells); % time of next spike (in rescaled time) 
rprev = zeros(1,ncells); % Integrated rescaled time up to current point
Itot = resp;
% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen); % Bins to update in this iteration
    nii = length(iinxt);  % Number of bins
    rrnxt = (nlfun(Itot(iinxt,:))-dc0)*dt; % Cond Intensity
    rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,ncells)],1);  % Cumulative intensity
    if all(tspnext >= rrcum(end,:)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end,:);
    else   % Spike!
        [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
        spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
        ispk = iinxt(min(ispks)); % time bin of spike(s)
        rprev = rrcum(min(ispks),:); % grab accumulated history to here

        % Record this spike
        mxi = min(rlen, ispk+hlen); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi;
        for ic = 1:length(spcells) % neurons that fire at time ispk
            icell = spcells(ic);
            nsp(icell) = nsp(icell)+1; % number of spikes
            tsp{icell}(nsp(icell),1) = ispk*dt; % spike timing
            if ~isempty(iiPostSpk)
                Itot(iiPostSpk,:) = Itot(iiPostSpk,:)+ih(1:mxi-ispk,:,icell);
            end
            rprev(icell) = 0;  % reset this cell's integral
            tspnext(icell) = exprnd(1); % draw RV for next spike in this cell
        end
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/(sum(nsp));
        nbinsPerEval = max(100, round(1.5*muISI)); 
    end
end

% Remove any extra bins from tsp and compute binned spike train 'sps'
sps = zeros(rlen,ncells);
for jj = 1:ncells
    tsp{jj} = tsp{jj}(1:nsp(jj));
    if ~isempty(tsp{jj})
        sps(round(tsp{jj}/dt),jj) = 1;
    end
end
simSpiketrain(perm,:,:) = sps;
spiketime(perm,:) = tsp;
end
% [simRate, simtime] = myPSTH(spiketime*1000, [0 rlen], binWidth, smoothSD);