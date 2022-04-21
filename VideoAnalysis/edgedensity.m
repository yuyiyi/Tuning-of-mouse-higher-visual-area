function [ed, Imgfilt] = edgedensity(Img, useGPU, filtersize, edge_thresh, edge_sigma)

if nargin==2
    filtersize = [];
end


[d1,d2, nFrames] = size(Img);
NFFT = 2.^nextpow2(max([d1,d2]));

if ~isempty(filtersize)
    h = fspecial('gaussian', NFFT, filtersize);
    d10 = floor((NFFT-d1)/2);
    d20 = floor((NFFT-d2)/2);
    if useGPU
        g = gpuDevice(1);
        reset(g)    
        batchSize = 2^20/2^ceil(log2(NFFT^2)); % works well on GTX 970
        hh = gpuArray(fft2(h));
    else
        batchSize = 20;
        hh = fft2(h);
    end

    %%%% filtering with difference of gaussian filter
    Imgfilt = zeros(d1, d2, nFrames, 'single');
    nBatches = ceil(nFrames/batchSize);

    for bi = 1:nBatches
        fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
        batchData = zeros(NFFT,NFFT,length(fi))+mean(Img(:));
        data_filt1 = zeros(NFFT,NFFT,length(fi));
        data1 = Img(:,:,fi);
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
        datatmp = gather(data_filt);
        data_filt1(1:NFFT/2,1:NFFT/2,:) = datatmp(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:);
        data_filt1(1:NFFT/2,1+NFFT/2:NFFT,:) = datatmp(1+NFFT/2:NFFT,1:NFFT/2,:);
        data_filt1(1+NFFT/2:NFFT,1:NFFT/2,:) = datatmp(1:NFFT/2,1+NFFT/2:NFFT,:);
        data_filt1(1+NFFT/2:NFFT,1+NFFT/2:NFFT,:) = datatmp(1:NFFT/2,1:NFFT/2,:);  
        Imgfilt(:,:,fi) = data_filt1(d10+1:d10+d1,d20+1:d20+d2,:);
    end
    if useGPU
        reset(g); 
    end    
else
    Imgfilt = Img;
end

ed = [];
if nargin<=3
    parfor i = 1:nFrames
        BW = edge(Imgfilt(:,:,i), 'Canny');         
        ed(i, 1) = sum(BW(:));
    end
else
    parfor i = 1:nFrames
        BW = edge(Imgfilt(:,:,i), 'Canny',edge_thresh, edge_sigma);
        ed(i, 1) = sum(BW(:));
    end    
end
ed = ed/(d1*d2);


