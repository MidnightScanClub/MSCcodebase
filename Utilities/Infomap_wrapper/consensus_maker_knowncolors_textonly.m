function consensus_maker_knowncolors_textonly(assigns_file,groupnetworksfile,basename)
%consensus_maker_knowncolors(regularized_ciftifile,groupnetworksfile,dostripes,mincol,minsize,orig_parcelsfile)

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('dostripes','var') || isempty(dostripes)
    dostripes = 0;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end



% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:100];

if ischar(assigns_file)
    assigns = load(assigns_file);
else
    assigns = assigns_file;
end
assigns(assigns<0) = 0;
assigns(isnan(assigns)) = 0;

if ischar(groupnetworksfile)
    groupdata = load(groupnetworksfile);
else
    groupdata = groupnetworksfile;
end
groupdata(groupdata<0) = 0;
groupdata(isnan(groupdata)) = 0;


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
                    
                    if c>2 && any(any(all_recolored(:,1:(c-1))==potential_colors(i)))
                        prevnetwork_comp = any(all_recolored(:,1:(c-1))==potential_colors(i),2);
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

outname = [basename '_allcolumns_recolored'];

dlmwrite([outname '.txt'],all_recolored,'delimiter',' ');


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




outname = [basename '_recolored'];
dlmwrite([outname '.txt'],out,'delimiter',' ');


