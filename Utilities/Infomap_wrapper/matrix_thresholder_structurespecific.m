function i = matrix_thresholder_structurespecific(matrix,threshold,structinds)
%i = matrix_thresholder_simple(matrix,threshold,structinds)
%
% This script takes a 2D matrix with defined structural divisions along
% each dimension (e.g., cortical vs subcortical) and, for each pairwise
% combination of structural divisions in the matrix (e.g.
% cortical-to-cortical, cortical-to-subcortical, etc.), thresholds it at a
% given edge density. The script then returns the indices of all cells in
% the upper triangle of the whole matrix that survive the
% structure-specific threshold, in descending order by their original
% matrix value. 
%
% NOTE:
% This script presumes an undirected (symmetric) matrix !!!!
% 
% It also ignores values on the diagonal (should be zero anyways) !!!!
%
% This can be very RAM-intensive for a large matrix....
%
% EMG 08/18/15

d=size(matrix);
if size(d,1)>2
    error('Matrix_thresholder only works on 2D matrices (with reason!)');
end

numpossibleedges=d(1)*(d(1)-1)/2;

inds = triu(true(size(matrix)),1);

structs = unique(structinds);
for structnum1 = 1:length(structs)
    for structnum2 = structnum1:length(structs)
        
        smallinds = inds; smallinds(structinds~=structs(structnum1),:) = 0; smallinds(:,structinds~=structs(structnum2)) = 0;
        
        smallvals = matrix(smallinds);
        
        [~, smallorder]=sort(smallvals,'descend');
                
        smallvals(smallorder) = [numel(smallvals):-1:1] / numel(smallvals);
        clear smallorder
        
        matrix(smallinds) = smallvals;
                
    end
end

matrix=triu(matrix,1);
[v i]=sort(matrix(:),'descend');
i = i(1:(round(threshold(end)*numpossibleedges)));

