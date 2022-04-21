function e = Feature_Entropy(data, x1, f, localmask)
if nargin <4
    [d, l] = size(data);
    e = [];    
    parfor i = 1:l
        tmp = data(:,i);
        tmp(isnan(tmp)) = [];
        h = hist(tmp, x1)/length(tmp);
        e(i,1) = -nansum(h.*log(h+eps));
    end
else
    e = [];    
    [d, l] = size(data);
    for k = 1:max(localmask(:))+1
        datatmp = data(localmask==k-1,:);
        datatmp = reshape(datatmp,[], l);
        parfor i = 1:l
            h = hist(datatmp(:,i), x1);
            e(k,i) = -nansum(h.*log(h+eps));
        end
    end
end