function video = readInStimulus(folder, filename)
[~, ~, extension] = fileparts(filename);

if strcmp(extension,'.tif')
    imglist = imfinfo(fullfile(folder, filename));
    video = [];
    for i = 1:length(imglist)
        tmp = imread(fullfile(folder, filename), i);
        video(:,:,i) = tmp(:,:,1);
    end
elseif strcmp(extension,'.avi')
    v = VideoReader(fullfile(folder, filename));
    video = []; k = 0;
    while hasFrame(v)
        tmp = readFrame(v);
        k=k+1;
        video(:,:,k) = tmp(:,:,1);
    end
end