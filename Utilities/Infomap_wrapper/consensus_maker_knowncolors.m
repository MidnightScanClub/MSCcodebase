function consensus_maker_knowncolors(regularized_ciftifile,manualset,groupnetworksfile,dostripes,mincol,minsize,orig_parcelsfile)
%consensus_maker_knowncolors(regularized_ciftifile,[manualset],[groupnetworksfile],[dostripes],[mincol],[minsize],[orig_parcelsfile])

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('dostripes','var') || isempty(dostripes)
    dostripes = 0;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end

if ~exist('groupnetworksfile','var') || isempty(groupnetworksfile)
    groupnetworksfile = '/data/nil-bluearc/GMT/Evan/MSC/Analysis_V2/Networks_template.dscalar.nii';
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:100];

cifti_data = ft_read_cifti_mod(regularized_ciftifile); assigns = cifti_data.data;
assigns(assigns<0) = 0;
assigns(isnan(assigns)) = 0;



groupfile = ft_read_cifti_mod(groupnetworksfile);
groupdata = groupfile.data;
ncortverts = nnz(groupfile.brainstructure==1) + nnz(groupfile.brainstructure==2);
groupdata = groupdata(1:ncortverts,1);

potential_colors = [1 2 10 9 3 5 6 11 16 15 7 8 17 12 4 14 13];
newcolors = setdiff(all_color_values,potential_colors);

unassigned_networks = cell(1,size(assigns,2));

all_recolored = zeros(size(assigns));
for c = 1:size(all_recolored,2)
    col_consensusmap = assigns(:,c);
    
    
    unassigned = find(col_consensusmap<1);
    for unassignedindex = unassigned'
        thisassignments = assigns(unassignedindex,mincol:end);
        thisassignments(thisassignments<1) = [];
        if ~isempty(thisassignments)
            col_consensusmap(unassignedindex) = thisassignments(1);
        end
    end
    networks = unique(col_consensusmap); networks(networks<=0) = [];
    new_networks = networks;
    assigning_networks = networks;
    
    col_out = zeros(size(col_consensusmap));
    
    if exist('manualset')
        if isempty(manualset)
            clear manualset
        else
            for i = 1:size(manualset,1)
                col_out(col_consensusmap==manualset(i,1)) = manualset(i,2);
                new_networks(new_networks==manualset(i,1)) = [];
                if all(manualset(i,2)~=potential_colors); new_networks = [new_networks; manualset(i,2)]; end
                assigning_networks(assigning_networks==manualset(i,1)) = [];
            end
        end
    end
    
    col_consensusmap_nosubcort = col_consensusmap(1:size(groupdata,1));
    
    for i = 1:length(potential_colors)
        
        if exist('manualset') && any(manualset(:,2)==potential_colors(i))
            
            
            
        else
            
            
            if ~isempty(assigning_networks)
                groupnetwork_comp = groupdata==potential_colors(i);
                D = zeros(length(assigning_networks),1);
                P = zeros(length(assigning_networks),1);
                for j = 1:length(assigning_networks)
                    
                    network_comp = col_consensusmap_nosubcort==assigning_networks(j);
                    P(j) = nnz(groupnetwork_comp & network_comp);
                    D(j) = P(j)/nnz(groupnetwork_comp | network_comp);
                    
                    if c>2 && any(any(all_recolored(1:ncortverts,1:(c-1))==potential_colors(i)))
                        prevnetwork_comp = any(all_recolored(1:ncortverts,1:(c-1))==potential_colors(i),2);
                        D_withprev = nnz(prevnetwork_comp & network_comp) ./ nnz(network_comp | prevnetwork_comp);
                        if D_withprev < .1
                            D(j) = 0;
                        end
                    end
                end
                [maxval, maxind(i)] = max(D);
                
                if maxval > .1
                    col_out(col_consensusmap==assigning_networks(maxind(i))) = potential_colors(i);
                    new_networks(new_networks==assigning_networks(maxind(i))) = [];
                    assigning_networks(assigning_networks==assigning_networks(maxind(i))) = [];
                end
            end
        end
        
    end
    clear maxind D P
    for j = 1:length(new_networks)
        col_out(col_consensusmap==new_networks(j)) = newcolors(j);
    end
    all_recolored(:,c) = col_out;
    unassigned_networks{c} = assigning_networks;
end


all_recolored(assigns<=0) = 0;
cifti_data.data = all_recolored;
if ~exist('cifti_data.mapname')
    for col = 1:size(all_recolored,2)
        cifti_data.mapname{col} = ['Column number ' num2str(col)];
    end
    cifti_data.dimord = 'scalar_pos';
end

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_allcolumns_recolored'];
ft_write_cifti_mod(outname,cifti_data);
set_cifti_powercolors([outname '.dscalar.nii'])




out = all_recolored(:,mincol);

uniquevals = unique(out); %uniquevals(uniquevals<1) = [];
colors_tofix = setdiff(uniquevals,potential_colors);
verts_tofix = [];
for colornum = 1:length(colors_tofix)
    verts_thiscolor = find(out==colors_tofix(colornum));
    verts_tofix = [verts_tofix ; verts_thiscolor(:)];
end
for vertnum = verts_tofix'
    for col = (mincol+1):size(all_recolored,2)
        if any(all_recolored(vertnum,col)==potential_colors)
            out(vertnum) = all_recolored(vertnum,col);
            break
        end
    end
end


temp_out = out;

if logical(minsize)
    
    
    % Clean up tiny pieces
    
    
    neighbors = cifti_neighbors(regularized_ciftifile);
    allcolors = unique(out); allcolors(allcolors<=0) = [];
    
    
    for color = allcolors(:)'
        clusteredmetric = zeros(size(temp_out));
        thiscolorverts = find(temp_out==color);
        for vertex = thiscolorverts'
            
            %find the neighbors of this vertex
            vertexneighbors = neighbors(vertex,:);
            vertexneighbors(isnan(vertexneighbors)) = [];
            
            %find which of those neighbors also pass the thresholds
            vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
            
            %find if those neighbors have already been assigned different cluster values
            uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
            uniqueneighborvals(uniqueneighborvals==0) = [];
            
            %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
            if isempty(uniqueneighborvals)
                clusteredmetric(vertexneighbors_thiscolor) = vertex;
                %if there is only one previous cluster identifier present, make all the neighbors that value
            elseif length(uniqueneighborvals)==1
                clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
                %if there are multiple cluster identifier values in the neighborhood, merge them into one
            else
                for valuenum = 2:length(uniqueneighborvals)
                    clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                end
            end
        end
        uniqueclustervals = unique(clusteredmetric);
        uniqueclustervals(uniqueclustervals==0) = [];
        
        for clusternum = uniqueclustervals'
            if nnz(clusteredmetric==clusternum) < minsize
                neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                neighborverts(isnan(neighborverts)) = [];
                borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                borderverts(temp_out(borderverts)<1) = [];
                mode_neighborval = mode(temp_out(borderverts));
                temp_out(clusteredmetric==clusternum) = mode_neighborval;
            end
        end
    end
    
    
    out = temp_out;
    
end


cifti_data.data = out;
if ~exist('cifti_data.mapname')
    cifti_data.mapname = {'Column number ' num2str(mincol)};
    cifti_data.dimord = 'scalar_pos';
else
    cifti_data.mapname = cifti_data.mapname(mincol);
end

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_recolored'];
ft_write_cifti_mod(outname,cifti_data);

set_cifti_powercolors([outname '.dscalar.nii'])

if exist('orig_parcelsfile') && ~isempty(orig_parcelsfile)
    parcels = ft_read_cifti_mod(orig_parcelsfile);
    parcels = parcels.data;
    parcels((length(out)+1):end) = [];
    IDs = unique(parcels); IDs(IDs<1) = [];
    outtext = zeros(length(IDs),1);
    outtext_bycol = zeros(length(IDs),size(all_recolored,2));
    for IDnum = 1:length(IDs)
        outtext(IDnum) = mode(out(parcels==IDs(IDnum)));
        outtext_bycol(IDnum,:) = mean(all_recolored(parcels==IDs(IDnum),:),1);
    end
    dlmwrite([outname '.txt'],outtext,'delimiter',' ')
    dlmwrite([basename '_allcolumns_recolored.txt'],outtext_bycol,'delimiter',' ')
end


%Stripes

if dostripes
    
    all_recolored = all_recolored(1:ncortverts,:);
    allcolors = unique(all_recolored); allcolors(allcolors==0) = [];
    unknown_colors = setdiff(allcolors,potential_colors);
    for color = unknown_colors(:)'
        all_recolored(all_recolored==color) = 0;
    end
    
    to_be_striped = out;
    change = logical(diff(all_recolored,1,2));
    change = change .* (all_recolored(:,1:end-1)>0) .* (all_recolored(:,2:end)>0);
    for col = 1:size(change,2)
        colvals = unique(all_recolored(:,col+1)); colvals(colvals==0) = [];
        for val_totest = colvals(:)'
            if nnz((all_recolored(:,col)==val_totest) & (all_recolored(:,col+1)==val_totest)) / nnz((all_recolored(:,col)==val_totest) | (all_recolored(:,col+1)==val_totest)) > .5
                thiscolor_changed_inds = logical(change(:,col) .* (all_recolored(:,col)==val_totest));
                to_be_striped(thiscolor_changed_inds,2+((col-1)*2)) = all_recolored(thiscolor_changed_inds,col+1);
                to_be_striped(thiscolor_changed_inds,3+((col-1)*2)) = val_totest;
            end
        end
    end
    
    
    to_be_striped_final = zeros(size(to_be_striped));
    for vert = 1:size(to_be_striped_final,1)
        uniquevals = unique([out(vert) to_be_striped(vert,:)]); uniquevals(uniquevals==0) = [];
        to_be_striped_final(vert,1:length(uniquevals)) = uniquevals;
    end
    make_striped_cifti(to_be_striped_final,0,[outname '_striped_164'],1/30)
    set_cifti_powercolors([outname '_striped_164.dtseries.nii'])
    
end
