function [P, D, h] = PartialCorr_SparseLatentNetwork(extras, NetworkDegree_bin, PartialCorr_threshold)

C = inv(extras.S - extras.H*extras.H');
sigma1 = diag(diag(inv(C)));
% partial correlation 
P = -sqrt(sigma1)\inv(C)/sqrt(sigma1); % parital correlation matrix

NeuNum = size(P,1);

if nargin < 2
    x = 0:0.1:1;
else
    x = NetworkDegree_bin;
end

if nargin < 3
    PartialCorr_threshold = 0.01;
end

% ensure partial correlation matrix is symmetric
P1 = sqrt(P.*P'); 
A = ones(size(P));
A(P<0) = -1;
P1 = P1.*A;
P1(abs(P1)<PartialCorr_threshold) = 0;
B = 1-diag(ones(size(P,1),1));
P1 = P1.*B;
G = graph(P1);

%%%% degree
D = degree(G);
Dtmp = D/NeuNum;
h = hist(Dtmp, x)/NeuNum;

