function plot3DGabor(SimpleUnit0, SimpleUnit90, filtsize)
figure(1),
subplot(1,2,1), imagesc(SimpleUnit0(:,:,1))
subplot(1,2,2), imagesc(SimpleUnit90(:,:,1))
figure(2), 
subplot(1,4,1), imagesc(squeeze(SimpleUnit0(round(filtsize/2),:,:)))
subplot(1,4,2), imagesc(squeeze(SimpleUnit0(:,round(filtsize/2),:)))
subplot(1,4,3), imagesc(squeeze(SimpleUnit90(round(filtsize/2),:,:)))
subplot(1,4,4), imagesc(squeeze(SimpleUnit90(:,round(filtsize/2),:)))
