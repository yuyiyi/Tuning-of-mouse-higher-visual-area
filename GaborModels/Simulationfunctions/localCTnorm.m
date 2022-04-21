function imn = localCTnorm(imw)

%%%% local contrast normalization
D = 16;
[x y] = meshgrid(-D/2:D/2-1);
G = exp(-0.5*((x.^2+y.^2)/(D/2)^2));
G = G/sum(G(:));
imv = conv2(imw.^2,G,'same');
imn = imw./sqrt(imv);