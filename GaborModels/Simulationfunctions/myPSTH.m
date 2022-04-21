function [binRates, binTimes] = myPSTH(spikes, timeEnds, binWidth, smoothSD)
%INPUTS - spikes- an n x m matrix, with m = number of trials and n =
% max number of spikes over all trials. Spike times (in ms) are listed
% along the columns for each trial, followed by zeros
% timeEnds- [timeStart timeEnd] (in ms) along which to compute
% binWidth- time binning (in ms) of output psth
% smoothSD- standard deviation of gaussian (in ms) used to smooth psth
%   something like 1-2 ms is standard
% OUTPUTS- binRates- output PSTH
% binTimes - time sampling points of binRates

numTrials = size(spikes,2);
startTime = timeEnds(1);
endTime = timeEnds(2);

% binWidth = 1;

binTimes = startTime:binWidth:endTime;

binCounts = histc(spikes,binTimes,1);
binSums = sum(binCounts,2);
binSums(1) = 0;
binProbs = binSums/numTrials;
binRates = binProbs/((.001)*binWidth);

binTimes = (binTimes(1:end-1)+binTimes(2:end))/2;
binRates = binRates(1:end-1);

if nargin > 3
    smoothTpts = -5:1:5;
    smoothVals = inv(2*pi*smoothSD^2)*exp((-(smoothTpts).^2)/(2*smoothSD^2));
    smoothVals = smoothVals/sum(smoothVals);
    
    binRatesSmooth = conv(binRates, smoothVals', 'same');
    binRates = binRatesSmooth;
end




