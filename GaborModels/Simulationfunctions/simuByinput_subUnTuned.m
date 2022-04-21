% GLM simulation 
function [spiketime, simSpiketrain, simRate] = simuByinput_subUnTuned(resp, dt, rep, nlfun, ih, resp_untuned, dc)
if nargin < 6
    resp_untuned = zeros(length(resp),1);
    dc = 1.5;
end
resp_untuned = max(0,resp_untuned);

rlen = length(resp);
simSpiketrain = zeros(rep, rlen);
spiketime = [];
nbinsPerEval = 100;
hlen = length(ih);
% hlen = 0;
for perm = 1:rep
    Itot = resp-resp_untuned;
    tspnext = exprnd(1);  % time of next spike (in rescaled time)
    rprev = 0;  % Integrated rescaled time up to current point
    jbin = 1;
    sps = zeros(1, rlen);
    nsp = 0;
% --------------- Run dynamics ---------------------------------------
while jbin < rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = (nlfun(Itot(iinxt))-dc)*dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if ~isempty(rrcum)
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
%             Itot(iiPostSpk) = Itot(iiPostSpk) - 10;
            Itot(iiPostSpk) = Itot(iiPostSpk)+ih(1:mxi-ispk);
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
    end
end
tsp = find(sps>0)*dt;
simSpiketrain(perm, 1:length(sps)) = sps;
spiketime = vectCat(spiketime, tsp');
end

binWidth = 50;
smoothSD = 2;
if nargout==3
    [simRate, ~] = myPSTH(spiketime*1000, [0 rlen*dt*1000], binWidth, smoothSD);
end