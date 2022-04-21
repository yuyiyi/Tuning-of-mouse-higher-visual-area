function [stim, stim_mask, stimtrain] = makeStimStructure(video, loc, filtsize, t_intgrat, ZeroStim)

if filtsize > min(size(video(:,:,1)))
    stim_padding = bsxfun(@times, ones(filtsize,filtsize,size(video,3)), mean(mean(video,1),2));    
    stim_mask = stim_padding;
    i1 = floor((filtsize - size(video, 1))/2);
    i2 = floor((filtsize - size(video, 2))/2);
    stim_padding(i1+1:i1+size(video, 1), i2+1:i2+size(video, 2), :) = video;
    stim_mask(i1+1:i1+size(video, 1), i2+1:i2+size(video, 2), :) = 1;
    stim = stim_padding;
elseif ~isempty(loc) && filtsize <= min(size(video(:,:,1)))
    r1 = loc(1);
    c1 = loc(2);
    stim = single(video(r1-filtsize/2:(r1+filtsize/2-1), c1-filtsize/2:(c1+filtsize/2-1),:));
end

if nargin<5
    base = mean(stim(:));
else
    base = mean(ZeroStim(:));
end
stim = stim-base;
stim = stim/std(stim(:));

if nargout==3
    stimtrain = []; 
    for i = 1:t_intgrat
        stimtrain(:,:,i, i:size(stim,3)) = stim(:,:,1:end-(i-1));
    end
end
