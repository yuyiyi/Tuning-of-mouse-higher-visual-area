function [Vx, Vy]= movie_opticalflow(Ifiltered)
opticFlow = opticalFlowHS;
clear movfilt_OF
for i = 1:size(Ifiltered,3)
    movfilt_OF(i) = estimateFlow(opticFlow,Ifiltered(:,:,i));
end
parfor i = 1:size(Ifiltered,3)
    Vx(:,:,i) = movfilt_OF(i).Vx;
    Vy(:,:,i) = movfilt_OF(i).Vy;
end
% 
% Vxmean = []; Vymean = []; vmean = [];
% parfor i = 1:size(Vx,3)
%     vabs(:,:,i) = sqrt(abs(Vx(:,:,i)).^2+abs(Vy(:,:,i)).^2);
%     v1 = vabs(:,:,i);
%     vmean(i) = mean(v1(v1<0.1));
%     s = nanmean(nanmean(Vx(:,:,i)))./vmean(i);
%     c = nanmean(nanmean(Vy(:,:,i)))./vmean(i);
%     if s<0
%         Vxmean(i) = -mean(mean(sqrt(abs(Vx(:,:,i)).^2)));
%     else
%         Vxmean(i) = mean(mean(sqrt(abs(Vx(:,:,i)).^2)));
%     end
%     if c<0
%         Vymean(i) = -mean(mean(sqrt(abs(Vy(:,:,i)).^2)));
%     else
%         Vymean(i) = mean(mean(sqrt(abs(Vy(:,:,i)).^2)));
%     end   
% end