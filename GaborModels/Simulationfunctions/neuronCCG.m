function [ccgRaw, ccgshiftpred, ccgRawstd, ccgpair] =...
    neuronCCG(spiketrainCCG1,spiketrainCCG2, dt, tau, cellforccg1,cellforccg2,figureflag, correctionflag)
% compute cross-correlagram as Bair, 2001

N1 = length(spiketrainCCG1);
N2 = length(spiketrainCCG2);
S = length(spiketrainCCG1(1).SP);
lagbin = (length(tau)-1)/2;
baseRegion = find(abs(tau)>900);
peakRegion = find(abs(tau<900));

clear ccgRaw ccgshiftpred ccgRawstd ccgpair
ind = 1;
for n1 = 1:N1
    for n2 = n1+1:N2;
        for s = 1:S
            sc1_1 = spiketrainCCG1(n1).SP{s};
            sc1_2 = spiketrainCCG2(n2).SP{s};
            [R,T] = size(sc1_1);
            phi = (T-abs(tau))/1000; % unit sec

            sc1_1 = reshape(sc1_1,R, dt, T/dt);
            y1 = squeeze(sum(sc1_1, 2));
            y1(y1>1) = 1;
            trialout1 = find(sum(y1,2)<1);
            
            sc1_2 = reshape(sc1_2,R, dt, T/dt);
            y2 = squeeze(sum(sc1_2, 2));
            y2(y2>1) = 1;
            trialout2 = find(sum(y2,2)<1);
            
            trialout = unique([trialout1;trialout2]);
            trialkeep = setdiff(1:R,trialout);
            if length(trialkeep)>20
                y1 = y1(trialkeep,:); y2 = y2(trialkeep,:);
                y1_psth = mean(y1,1); y2_psth = mean(y2,1);
                meanfr1 = sum(y1_psth)/(T/1000); % unit spikes/sec 
                meanfr2 = sum(y2_psth)/(T/1000);
                ccg1perm = [];
                parfor perm = 1:100
                    permID = ceil(rand(round(length(trialkeep)*0.8),1)*length(trialkeep));
                    y1tmp = y1(permID,:); y2tmp = y2(permID,:);
                    ccg1 = [];
                    for r = 1:length(permID)                          
                        ccg1(r,:) = xcorr(y1tmp(r,:)',y2tmp(r,:)',lagbin);
                    end
                    ccg1perm(perm,:) = nanmean(ccg1)/sqrt(meanfr1*meanfr2)./phi; % unit coincidence/spike
                end
                ccgRaw(ind,:) = nanmean(ccg1perm,1);
                ccgRawstd(ind,:) = nanstd(ccg1perm,[],1);
                
                % correction shift predictor
                if correctionflag == 2
                    ccg2 = [];
                    for r = 1:length(trialkeep)-1
                        ccg2(r,:) = xcorr(y1(r,:)',y2(r+1,:)',lagbin);
                    end
                    ccgshiftpred(ind,:) = nanmean(ccg2)/sqrt(meanfr1*meanfr2)./phi; % unit coincidence/spike
                end
                % correction by all-way shuffle correction
                if correctionflag == 1
                    ccg2 = [];
                    for r = 1:length(trialkeep)
                        ccg2(r,:) = xcorr(y1(r,:)',y2(r,:)',lagbin);
                    end
                    ccg2 = nanmean(ccg2)/sqrt(meanfr1*meanfr2)./phi; % unit coincidence/spike
                    ccg1shift = (xcorr(y1_psth', y2_psth',lagbin )'*length(trialkeep)...
                        - nanmean(ccg2))/(length(trialkeep)-1);   
                    ccgshiftpred(ind,:) = ccg1shift/sqrt(meanfr1*meanfr2)./phi;
                end
                ccgpair(ind,:) = [cellforccg1(n1),cellforccg2(n2),s];
                ind = ind+1;
            end
        end
    end
end

if ind==1
    ccgRaw = [];
    ccgshiftpred = [];
    ccgRawstd = [];
    ccgpair = [];
end   
    
if figureflag==1 && ind>1
for c = 1:ceil(size(ccgRaw,1)/20)
    ii = (c-1)*20+1:min(c*20,size(ccgRaw,1));
    for i = 1:length(ii)
        j = ii(i);
        figure(1), subplot(4,5,i); errorbar(tau, smooth(ccgRaw(j,:)), ccgRawstd(j,:)/sqrt(100))
        hold on, plot(tau, smooth(ccgshiftpred(j,:))); title(num2str(ccgpair(j,3)))
        %axis([-1000 1000 -0.01 0.1])
        figure(2), subplot(4,5,i);
        errorbar(tau, smooth(ccgRaw(j,:)-ccgshiftpred(j,:)),ccgRawstd(j,:)/sqrt(50),'k')
        %axis([-1000 1000 -0.01 0.05])
    end
    waitforbuttonpress
    figure(1), clf('reset')
    figure(2), clf('reset')
end
end
