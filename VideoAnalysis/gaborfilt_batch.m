function [data_gabor_HV, data_gabor_diag, data_filt_norm] = gaborfilt_batch(data, gabor_wavelength, gabor_orientation, useGPU)
%%
gabor_filters = gabor(gabor_wavelength, gabor_orientation);
filtersize = size(gabor_filters(1).SpatialKernel,1);
[d1,d2,nFrames] = size(data);
NFFT = 2.^nextpow2(max([d1,d2]));
d10 = floor((NFFT-d1)/2);
d20 = floor((NFFT-d2)/2);

% h = zeros(NFFT, NFFT, length(gabor_filters)); 
h_fft = []; gaussfilt = [];
% dh = floor((NFFT-filtersize)/2);
for i = 1:length(gabor_filters)
    h_fft(:,:,i) = makeFrequencyDomainTransferFunction(gabor_filters(i),[NFFT, NFFT],'double');
%     h(dh+1:dh+filtersize,dh+1:dh+filtersize,i) = real(gabor_filters(i).SpatialKernel);
%     h_fft(:,:,i) = fft2(h(:,:,i));
end

if useGPU
    g = gpuDevice(1);
    reset(g)    
    batchSize = 2^22/2^ceil(log2(NFFT^2)); % works well on GTX 970
else
    batchSize = 20;
end

%%%% filtering with gabor filter
data_gabor_HV = zeros(d1, d2, nFrames, 'single');
data_gabor_diag = zeros(d1, d2, nFrames, 'single');
nBatches = ceil(nFrames/batchSize);

for bi = 1:nBatches
  fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
  batchData = zeros(NFFT,NFFT,length(fi))+mean(data(:));
  data_filt1 = zeros(NFFT,NFFT,length(fi));
  data1 = data(:,:,fi);
  data_gaborfilt = zeros(d1, d2, length(fi), length(gabor_filters), 'single');
  batchData(d10+1:d10+d1,d20+1:d20+d2,:) = single(data1);
  batchData(1:d10,d20+1:d20+d2,:) = repmat(data1(1,:,:),d10,1,1);
  batchData(d10+d1+1:NFFT,d20+1:d20+d2,:) = repmat(data1(d1,:,:),NFFT-(d10+d1),1,1);
  batchData(d10+1:d10+d1,1:d20,:) = repmat(data1(:,1,:),1,d20,1);
  batchData(d10+1:d10+d1,d20+d2+1:NFFT,:) = repmat(data1(:,d2,:),1, NFFT-(d20+d2),1);   
  for k = 1:length(gabor_filters)
    if useGPU
        hh = gpuArray(h_fft(:,:,k));
        batchData = gpuArray(batchData);
    else
        hh = h_fft(:,:,k);
    end
  data_fft = fft2(batchData);
  data_filt = abs(ifft2(bsxfun(@times, data_fft, ifftshift(hh))));
  datatmp = gather(data_filt);
%   data_filt1(1:NFFT/2,1:NFFT/2,:) = datatmp(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:);
%   data_filt1(1:NFFT/2,1+NFFT/2:NFFT,:) = datatmp(1+NFFT/2:NFFT,1:NFFT/2,:);
%   data_filt1(1+NFFT/2:NFFT,1:NFFT/2,:) = datatmp(1:NFFT/2,1+NFFT/2:NFFT,:);
%   data_filt1(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:) = datatmp(1:NFFT/2,1:NFFT/2,:);  
  data_gaborfilt(:,:,:,k) = datatmp(d10+1:d10+d1,d20+1:d20+d2,:);
  end
  data_filt_energy = data_gaborfilt.^2;
  %%%% contrast invariant
  data_filt_norm = bsxfun(@rdivide, data_filt_energy, sum(data_filt_energy,4)); % contrast invariance
%   tmp = reshape(data_filt_norm,[],length(gabor_filters));
%   [vpref, oripref] = max(tmp,[],2);
%   orinull = mod(oripref+length(gabor_filters)/2,length(gabor_filters));
%   orinull(orinull==0) = length(gabor_filters);
%   aa1 = (orinull-1)*length(orinull);
%   aa2 = 1:length(orinull);
%   aa3 = aa1+aa2';
%   vnull = tmp(aa3);
%   osi = (vpref-vnull)./(vpref+vnull);
%   osi = reshape(osi, d1,d2,length(fi));
%   oripref = cos(deg2rad(gabor_orientation(oripref)));
%   oripref = reshape(oripref,d1,d2,length(fi));
  data_gabor_HV(:,:,fi) = data_filt_norm(:,:,:,1) - data_filt_norm(:,:,:,3);  
  data_gabor_diag(:,:,fi) = data_filt_norm(:,:,:,2) - data_filt_norm(:,:,:,4);
%   data_gabor(:,:,fi) = data_filt_norm1 + data_filt_norm2;
end

if useGPU
    reset(g); 
end
