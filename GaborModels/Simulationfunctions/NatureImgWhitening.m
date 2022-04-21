function imw = NatureImgWhitening(I)
I = double(I);
%%%% 2D FFT of square image
NFFT = 2.^nextpow2(size(I));
imf = fftshift(fft2(I,NFFT(1), NFFT(2))); % FT of image I
% amplitude and phase
impha = angle(imf); % phase of frequency domain 
imamp = abs(imf); % amplitude of frequency domain 
% power spectrum
impf = abs(imf).^2;
f = -NFFT/2:NFFT/2-1;

Pf=rotavg(impf); % angular average of power spectrum
f1 = 0:NFFT(1)/2;

%%%% whitening
[fx fy] = meshgrid(f);
[theta rho]=cart2pol(fx,fy); % polar coordinates of frequencies
filtf = rho.*exp(-0.5*(rho/(0.7*NFFT(1)/2)).^2); % 

imwf = filtf.*imf; % normalize by average power spectrum
imw = ifft2(fftshift(imwf));
