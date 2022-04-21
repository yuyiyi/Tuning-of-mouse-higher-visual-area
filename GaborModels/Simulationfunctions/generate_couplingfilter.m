function [Ih_cpl, ihbasis, Ih_para] = generate_couplingfilter(NeuNum, dt, ch_para, display)
%%%%
% NeuNum: number of neurons 
% hyperparameter for coupling filter, default = 0.1; a NeuNum * NeuNum matrix
% or a scaler. regulate correlation value
% dt: simulation time 
% display coupling filter
% --- Make basis for post-spike (h) filter ------
if nargout<3
    ch_para = 0;
end

ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [dt, dt*30];  % Peak location for first and last vectors
ihbasprs.b = dt*ihbasprs.ncols;  % Determines how nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional)
[iht,~,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);

%%%% generate coupling filter and post spike filter using basis
% coupling filter parameter vs. coupling filter shape
Ih_para = bsxfun(@plus, bsxfun(@times,randn(NeuNum, NeuNum, ihbasprs.ncols),shiftdim([3, 1, 0.2, 0.1, 0.01], -1)), ch_para);
Ih_para2 = bsxfun(@times, shiftdim(Ih_para,2), shiftdim(log(10.^rand(NeuNum, NeuNum)),-1));
Ih_cpl = [];
for k = 1:size(ihbasis,1)
    Ih_cpl(k,:,:) = sum(bsxfun(@times, ihbasis(k,:)', Ih_para2),1);
    % 3rd dimension is source
end
for n = 1:NeuNum
    Ih_cpl(:,n,n) = ihbasis*[-10 -4 -2 .25 -1]';
end

if nargout < 4    
    display = false;
end

if display
    figure,
    subplot(121), plot(Ih_cpl(:,:,1))
    subplot(122), plot(Ih_cpl(:,:,2))
end
