function data_DOG = DOG_batch(data, filtersize, useGPU, filter_weights)
%%
if nargin < 4
    filter_weights = [1,1];
end
[d1,d2,nFrames] = size(data);
NFFT = 2.^nextpow2(max([d1,d2]));
h1 = fspecial('gaussian', NFFT, filtersize);
h2 = fspecial('gaussian', NFFT, filtersize*2);
d10 = floor((NFFT-d1)/2);
d20 = floor((NFFT-d2)/2);
if useGPU
    g = gpuDevice(1);
    reset(g)    
    batchSize = 2^20/2^ceil(log2(NFFT^2)); % works well on GTX 970
    h1_f = gpuArray(fft2(h1));
    h2_f = gpuArray(fft2(h2));
    hh = filter_weights(1)*h1_f-filter_weights(2)*h2_f;
else
    batchSize = 20;
    h1_f = fft2(h1);
    h2_f = fft2(h2);
    hh = filter_weights(1)*h1_f-filter_weights(2)*h2_f;
end

%%%% filtering with difference of gaussian filter
data_DOG = zeros(d1, d2, nFrames, 'single');
nBatches = ceil(nFrames/batchSize);

for bi = 1:nBatches
  fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
  batchData = zeros(NFFT,NFFT,length(fi))+mean(data(:));
  data_filt1 = zeros(NFFT,NFFT,length(fi));
  data1 = data(:,:,fi);
  if useGPU
      batchData(d10+1:d10+d1,d20+1:d20+d2,:) = single(data1);
      batchData(1:d10,d20+1:d20+d2,:) = repmat(data1(1,:,:),d10,1,1);
      batchData(d10+d1+1:NFFT,d20+1:d20+d2,:) = repmat(data1(d1,:,:),NFFT-(d10+d1),1,1);
      batchData(d10+1:d10+d1,1:d20,:) = repmat(data1(:,1,:),1,d20,1);
      batchData(d10+1:d10+d1,d20+d2+1:NFFT,:) = repmat(data1(:,d2,:),1, NFFT-(d20+d2),1);      
      batchData = gpuArray(batchData);
  else
      batchData(d10+1:d10+d1,d20+1:d20+d2,:) = single(data1);
      batchData(1:d10,d20+1:d20+d2,:) = repmat(data1(1,:,:),d10,1,1);
      batchData(d10+d1+1:NFFT,d20+1:d20+d2,:) = repmat(data1(d1,:,:),NFFT-(d10+d1),1,1);
      batchData(d10+1:d10+d1,1:d20,:) = repmat(data1(:,1,:),1,d20,1);
      batchData(d10+1:d10+d1,d20+d2+1:NFFT,:) = repmat(data1(:,d2,:),1, NFFT-(d20+d2),1);      
  end
  data_fft = fft2(batchData);
  data_filt = real(ifft2(bsxfun(@times, data_fft, hh)));
  datatmp = gather(ifftshift(data_filt));
%   data_filt1(1:NFFT/2,1:NFFT/2,:) = datatmp(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:);
%   data_filt1(1:NFFT/2,1+NFFT/2:NFFT,:) = datatmp(1+NFFT/2:NFFT,1:NFFT/2,:);
%   data_filt1(1+NFFT/2:NFFT,1:NFFT/2,:) = datatmp(1:NFFT/2,1+NFFT/2:NFFT,:);
%   data_filt1(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:) = datatmp(1:NFFT/2,1:NFFT/2,:);  
  data_filt1 = datatmp;  
  data_DOG(:,:,fi) = data_filt1(d10+1:d10+d1,d20+1:d20+d2,:);
end

if useGPU
    reset(g); 
end