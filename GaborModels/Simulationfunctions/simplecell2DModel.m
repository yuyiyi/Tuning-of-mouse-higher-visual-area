function simSimpleUnit = simplecell2DModel(Gaussbandwidth, orientation, gaborphase, gaborwavelen, filtsize, pixelPerDeg)

% filtsize = 48;                           % image size: n X n pixel
% gaborwavelen = 5;                        % wavelength (number of deg per cycle)
% orientation = 0;                              % grating orientation
% Gaussbandwidth = 5;                             % gaussian standard deviation in deg
% gaborphase = deg2rad(90);                            % phase (0 -> 1)
trim = .005;                             % trim off gaussian values smaller than this
% pixelPerDeg = 2;

cycPerpixel = (1/gaborwavelen)/pixelPerDeg; % cycle per pixel 
X = 1:filtsize;                           % X is a vector from 1 to imageSize
X0 = ((X/filtsize)-0.5)*filtsize;         
[Xm, Ym] = meshgrid(X0, X0);             % 2D matrices

thetaRad = deg2rad(orientation);        % convert theta (orientation) to radians
X1 = Xm * cos(thetaRad)*2*pi;                % compute proportion of Xm for given orientation
Y1 = Ym * sin(thetaRad)*2*pi;                % compute proportion of Ym for given orientation
grating = sin(X1*cycPerpixel+Y1*cycPerpixel + deg2rad(gaborphase)); 

s = Gaussbandwidth*pixelPerDeg;                     % gaussian width as fraction of imageSize
gaussfilt = 1/(s*sqrt(2*pi))*exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
gaussfilt(gaussfilt < trim) = 0;                 % trim around edges (for 8-bit colour displays)
simSimpleUnit = grating .* gaussfilt;                % use .* dot-product

