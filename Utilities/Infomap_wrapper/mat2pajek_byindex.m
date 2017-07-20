function mat2pajek_byindex(mat,ind,outputname)
%mat2pajek_byindex(mat,ind,outputname)
%
% This takes a matrix, a list of indices within the matrix that survive
% some threshold, and a name for the output, and writes the graph in the
% pajek format, which is used in Pajek (for windows only), and also for
% some compiled versions of things, like Infomap. 
% 
% NOTE: pajek files should have .net extensions
%
%TOL 06/26/15; modified from JDP 10/10/10


% only the upper triangle
mat=triu(mat,1);

% get edges and values
[x y] = ind2sub(size(mat),ind);
z = mat(ind);

towrite = [x y z];
%%% make the input file %%%
nodes = size(mat,1);
nodenum = 1:nodes;

c=clock; 
fprintf('\t%2.0f:%2.0f:%2.0f: mat2pajek: writing .net file, with %d vertices and %d edges\n',c(4),c(5),c(6),nodes,length(x));

fid=fopen(outputname,'W');
fprintf(fid,'*Vertices %d\n',nodes);
fprintf(fid,'%d "%d"\n',[nodenum; nodenum]);
fprintf(fid,'*Edges %d\n',length(x));

fprintf(fid,'%d %d %f\n',towrite');

fclose(fid);