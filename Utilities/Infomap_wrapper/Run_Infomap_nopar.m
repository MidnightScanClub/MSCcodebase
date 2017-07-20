function Run_Infomap_nopar(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, structure_indices)
%Run_Infomap(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, [structure_indices])
%
% Run infomap on a matrix with a given distance exclusion, at various
% density thresholds, and write the results from all thresholds into a
% single text file named "rawassn.txt". This can take a long time for large
% matrices. It will run up to eight infomaps simultaneously if the Parallel
% Computing Toolbox is installed.
%
% Inputs:
%
% rmat - a correlation matrix to be infomapped. Can be a numeric matrix or
%  a cifti file that will be loaded
% dmatname - a .mat file containing a node-to-node distance matrix
% xdistance - the distance exclusion to apply, in mm (i.e., nodes closer
%  than xdistance are not allowed to be connected)
% thresholdarray - a vector of thresholds to apply to the matrix. Infomap
%  will be run separately for each threshold.
% makebinary - whether the matrix is binarized after thresholding. 1 = make
%  it binary; 0 = leave it weighted.
% outdir - the folder results will be written to. Will be created if it
%  doesn't exist.
%
% structure_indices - an OPTIONAL vector input of the same length as one
%  dimension of rmat that defines known divisions of interest in the matrix.
%  For example, one might provide a vector defining cortical nodes as "1"
%  and subcortical nodes as "2". If this input is provided, connection
%  values in the matrix will be thresholded separately for each combination
%  of defined structures. This is useful if one wants to e.g. allow
%  subcortical nodes to participate in communities even though they 
%  have very weak connection strengths.
%
%
% Requires the Resources scripts to be in your path (e.g.,
% /home/data/scripts/Resources/ and subfolders)
%
%EMG 06/25/15

dlmwrite([outdir '/thresholds.txt'],thresholdarray,'delimiter',' ')

[~,~] = system(['rm ' outdir '/pajek*']);

prevstring = [];

string = ['loading correlations...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

if ischar(rmat)
    rmat = ft_read_cifti_mod(rmat);
    rmat = rmat.data;
end
rmat = single(rmat);

warning off

string = ['applying distance exclusion...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;
% tic

if ischar(dmatname)
    dmat = smartload(dmatname);
elseif isnumeric(dmatname)
    dmat = dmatname;
    clear dmatname
end

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        rmat(dmat < xdistance) = 0;
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

clear dmat
%toc


numanalyses = length(thresholdarray);
numnodes = size(rmat,1);

if ~exist(outdir)
    mkdir(outdir);
end


string = ['finding r thresholds...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

%Remove diagonals
for i = 1:numnodes
    rmat(i,i) = 0;
end

%tic
if exist('structure_indices','var') && (numel(unique(structure_indices)) > 1)
    ind = matrix_thresholder_structurespecific(rmat,thresholdarray(end),structure_indices);
else
    ind = matrix_thresholder_simple(rmat,thresholdarray(end));
end
%toc

if makebinary
    rmat=ceil(rmat);
end

string = ['saving pajek files...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [outdir '/pajek_col' num2str(length(thresholdarray)) '.net'];
mat2pajek_byindex(rmat,ind,pajekfileorig)
numpossibleedges=(numnodes*(numnodes-1))/2;

clear rmat

% Write out other thresholds
for i = 1:numanalyses
    if i<numanalyses
        pajekfile = [outdir '/pajek_col' num2str(i) '.net'];
        edgesleft=round(thresholdarray(i)*numpossibleedges);
        numuse = edgesleft + numnodes + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
        evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]);
    end
end
    
string = ['running infomap'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
disp(' ')

% do analyses at each threshold
for i=1:numanalyses
    %fprintf('Thr/box %d, pass %d\n',i);
    
    %tic
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    rawclrs = run_infomap_on_pajekfile(pajekfile,100);
    dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
    %toc
end


for i = 1:numanalyses
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end



for i = 1:numanalyses
    rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
end
% write the raw assignments as .txt
dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
delete([outdir '/rawassn_col*.txt'])
    
end


