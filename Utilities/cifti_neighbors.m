function neighbors = cifti_neighbors(ciftifile,volneighbor_method,separate_structures)
%neighbors = cifti_neighbors(ciftifile,[volneighbor_method],[separate_structures])
%
%Determines adjacencies of data points in a cifti file
%
% 'ciftifile' is the file to calculate adjacencies from.
%
% 'volneighbor_method' determines what sort of adjacencies to allow in
% the volumetric elements of the file. Options are 'faces' (only voxels
% sharing a face count as neighbors), 'edges' (voxels sharing an edge count
% as neighbors), and 'corners' (voxels sharing a corner count as
% neighbors). Omit to default to 'edges', which has proved most useful in
% volumetric parcellations.
%
% 'separate_structures' is a binary argument determining whether neighbors
% are allowed across different but adjacent volumetric structures (e.g. do
% hippocampal voxels have amygdala neighbors). Set to 1 or omit to prevent
% neighbors from existing across different structures.
%
% EMG 06/24/15


if ~exist('volneighbor_method')
    volneighbor_method = 'edges';
end

if ~exist('separate_structures')
    separate_structures = 0;
end

if isstruct(ciftifile)
    cifti = ciftifile; clear ciftifile;
else
    cifti = ft_read_cifti_mod(ciftifile); 
end
    cifti.data = [];

bufsize=16384;
surfneighborfile = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/node_neighbors.txt';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[surfneighbors(:,1) surfneighbors(:,2) surfneighbors(:,3) surfneighbors(:,4)...
    surfneighbors(:,5) surfneighbors(:,6) surfneighbors(:,7)] = ...
    textread(surfneighborfile,'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
surfneighbors = surfneighbors+1;

nsurfverts = size(surfneighbors,1);

for hem = 1:2
    
    newordering = ones(nsurfverts,1) * NaN;
    verts_indata = (cifti.brainstructure((1:nsurfverts)+(nsurfverts*(hem-1)))==hem);
    verts_notindata = (cifti.brainstructure((1:nsurfverts)+(nsurfverts*(hem-1)))==-1);
    
    newordering(verts_indata) = [1:nnz(verts_indata)] + (nnz(cifti.brainstructure==1)*(hem-1));
    
    naninds = find(isnan(newordering));
    surfneighbors_mod = surfneighbors;
    surfneighbors_mod(isnan(surfneighbors_mod)) = naninds(1);
    
    hem_surfneighbors{hem} = newordering(surfneighbors_mod);
    
    hem_surfneighbors{hem}(verts_notindata,:) = [];
    
    hem_surfneighbors{hem}(:,end+1:27) = NaN;
    
end

neighbors = [hem_surfneighbors{1} ; hem_surfneighbors{2}];


if isfield(cifti,'transform')
    
    dims = [cifti.transform(1,1) cifti.transform(2,2) cifti.transform(3,3)];
    switch volneighbor_method
        case 'faces'
            toofar_tobeneighbors_dist = min([sqrt(dims(1)^2 + dims(2)^2) sqrt(dims(1)^2 + dims(3)^2) sqrt(dims(2)^2 + dims(3)^2)]);
        case 'edges'
            toofar_tobeneighbors_dist = sqrt(dims(1)^2 + dims(2)^2 + dims(3)^2);
        case 'corners'
            toofar_tobeneighbors_dist = sqrt(dims(1)^2 + dims(2)^2 + dims(3)^2) + .01;
    end
    volpos = cifti.pos(cifti.brainstructure>2,:);
    volstructs = cifti.brainstructure(cifti.brainstructure>2,:);
    
    [distances, volneighbors] = pdist2(volpos,volpos,'Euclidean','Smallest',27); 
    volneighbors = volneighbors';
    distances = distances';
    
    
    if separate_structures
        for vox = 1:size(volneighbors,1)
            thisvox_neighbors = volneighbors(vox,:);
            volneighbors(vox,volstructs(thisvox_neighbors)~=volstructs(vox)) = NaN;
        end
    end
    
    volneighbors(distances >= toofar_tobeneighbors_dist) = NaN;
    volneighbors = volneighbors + nnz(cifti.brainstructure==1) + nnz(cifti.brainstructure==2);
    
    neighbors = [neighbors; volneighbors];
else
    neighbors = neighbors(:,1:7);
end

