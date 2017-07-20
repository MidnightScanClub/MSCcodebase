function [mean_PC_pctiles, mean_PCs, all_PCs] = PC_calc(all_parcel_corrmat,parcel_distances_tooclose,communities,thresholdarray,highdegree,borderverts)

all_parcel_corrmat(logical(parcel_distances_tooclose)) = 0;
parcel_correlation_vec = all_parcel_corrmat(triu(true(size(all_parcel_corrmat)),1));
sorted_vals = sort(parcel_correlation_vec,'descend');

all_PCs = zeros(size(all_parcel_corrmat,1),length(thresholdarray));

if size(communities,2)==1
    communities = repmat(communities,1,length(thresholdarray));
end

%PCvals = zeros(size(communities));

for threshnum = 1:length(thresholdarray)
    
    thresh_communities = unique(communities(:,threshnum));
    thresh_communities(thresh_communities<1) = [];
    kdenthresh = thresholdarray(threshnum);
    rthresh = sorted_vals(ceil(length(sorted_vals) .* kdenthresh));
    corrmat_weighted = all_parcel_corrmat .* (all_parcel_corrmat > rthresh);
    
    parcel_degree = sum(corrmat_weighted,2);
    moduleconnectionssum = zeros(size(corrmat_weighted,1),1);
    
    for communitynum = 1:length(thresh_communities)
        communityID = thresh_communities(communitynum);
        communityindices = communities(:,threshnum)==communityID;
        parcel_community_degree = sum(corrmat_weighted(:,communityindices),2);
        parcel_community_ratio = (parcel_community_degree ./ parcel_degree);
        moduleconnectionssum = moduleconnectionssum + (parcel_community_ratio.^2);
    end
    
    all_PCs(:,threshnum) = 1-moduleconnectionssum;
end
all_PCs(~highdegree,:) = 0;
if exist('borderverts')
    if size(borderverts,2)==1
        all_PCs(borderverts,:) = 0;
    else
        all_PCs(borderverts) = 0;
    end
end
mean_PCs = nanmean(all_PCs,2);

all_PC_pcts = zeros(size(all_PCs));
for i = 1:size(all_PCs,2);
    nonanvec = all_PCs(:,i); nonanvec(isnan(nonanvec)) = 0;
    all_PC_pcts(:,i) = calc_percentiles(nonanvec);
end
all_PC_pcts(isnan(all_PCs)) = NaN;

mean_PCs_temp = nanmean(all_PC_pcts,2);

mean_PCs_temp(isnan(mean_PCs_temp)) = 0;

mean_PC_pctiles = calc_percentiles(mean_PCs_temp)';

end