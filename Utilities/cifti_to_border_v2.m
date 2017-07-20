function cifti_to_border_v2(ciftifile,filled,allcolors_inonemap,bordercolors)
%cifti_to_border_v2(ciftiname,filled,allcolors_inonemap,bordercolors)
%
%Make left and right surface border files from a cifti file
%
%ciftiname - the path to the cifti
%filled - the script can handle both cifti files that already look like
% borders (one-vertex-wide lines) as well as files that have filled blobs
% (the script will draw borders around the blobs). Set filled to zero for
% the former; set it to 1 for the latter.
%allcolors_inonemap - the script can handle both a multi-map cifti for which
% borders should be drawn separately for each map, as well as a single-map
% cifti with multiple nonzero values in the one map (i.e. a typical network map).
% Set allcolors_inonemap to 0 for the former; set it to 1 for the latter.
%bordercolors - by default (i.e. if bordercolors is not specified), all
% borders will be made black. If you don't want black borders, you can a)
% set bordercolors to be an rgb triple; b) set bordercolors to be an N by 3
% matrix of N rgb triples representing appropriate colors for the N maps /
% N discrete values in the cifti; or c) set bordercolors to 'default',
% which attempts to map colors to the standard values for the power_surf
% pallete. With the 'default' option, any values that are not integers
% between 1 and 17 will be set to "Unknown" and colored gray.
%
%E. Gordon 07/10/16



data_template = ft_read_cifti_mod(ciftifile);
data = data_template.data;

num_wholesurfverts = nnz(data_template.brainstructure <= 2);
data_wholesurf = zeros(num_wholesurfverts,size(data,2));
wholesurf_indswithdata = data_template.brainstructure(1:num_wholesurfverts) > 0;
data_wholesurf(wholesurf_indswithdata,:) = data(1:nnz(wholesurf_indswithdata),:);

data_bothhems = data_wholesurf;
IDs = unique(data_bothhems); IDs(IDs==0) = [];
    
if ~exist('allcolors_inonemap'); allcolors_inonemap = 0; end
if logical(allcolors_inonemap)
    temp = data_bothhems;
    data_bothhems = zeros(size(temp,1),length(IDs));
    for IDnum = 1:length(IDs)
        data_bothhems(:,IDnum) = temp==IDs(IDnum);
    end
    clear temp
end

if isfield(data_template,'mapname') && ~logical(allcolors_inonemap);
    colnames_bothhem = data_template.mapname;
else
    colnames_bothhem = cellstr(num2str(IDs));
end


if ~exist('bordercolors')
    bordercolors = [0 0 0];
end

if ischar(bordercolors) && strcmp(bordercolors,'default') 
    if any(IDs < 1) || any(IDs > 18) || any(mod(IDs,1))
        disp('Warning: Use of (default) power colors requires values to be integers between 1 and 18. Non-compliant values will be set to "Unknown".')
        data_bothhems(:,IDs<1) = 18;%[];
        IDs(IDs<1) = 18;%[];
        data_bothhems(:,IDs>18) = 18;%[];
        IDs(IDs>18) = 18;%[];
        data_bothhems(:,logical(mod(IDs,1))) = 18;%[];
        IDs(logical(mod(IDs,1))) = 18;%[];
    end
    defaultcolors = [1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.85 .85 .85;0 .5 0;.2 .2 .2;.5 .75 .2];
    bordercolors = defaultcolors(IDs,:);
    colnames_bothhem = {'Default','Visual','Fronto-Parietal','PrimaryVisual','DorsalAttention','PreMotor','VentralAttentionLanguage','Salience','Cingulo-Opercular','MotorHand','MotorMouth','Auditory','AntMedialTemporal','PostMedialTemporal','Cingulo-Parietal','Parieto-Occipital','MotorFoot','Unknown'};
    colnames_bothhem = colnames_bothhem(IDs);
end

if size(bordercolors,1) == 1;
    bordercolors = repmat(bordercolors,size(data_bothhems,2),1);
end

if size(bordercolors,1) ~= size(data_bothhems,2)
    error('Border colors variable has the wrong number of colors!')
end




clear neighbors

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


hemcaps = {'LEFT', 'RIGHT'};

MNI{1} = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii']); MNI{1} = MNI{1}.vertices;
MNI{2} = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii']); MNI{2} = MNI{2}.vertices;


for hem = 1:2
    
    data = data_bothhems([1:size(MNI{hem},1)] + (size(MNI{1},1) * (hem-1)),:);
    
    cols_withdata = logical(sum(data,1));
    data = data(:,cols_withdata);
    colnames = colnames_bothhem(cols_withdata);
    thishem_bordercolors = bordercolors(cols_withdata,:);
    
    
    if exist('filled') && logical(filled)
        data_temp = zeros(size(data));
        [verts,cols] = find(data);
        for vertnum = 1:size(verts)
            vertneighbors = neighbors(verts(vertnum),:); vertneighbors(isnan(vertneighbors)) = [];
            if any(data(vertneighbors,cols(vertnum))==0)
                data_temp(verts(vertnum),cols(vertnum)) = 1;
            end
        end
        data = data_temp;
    end
    
    
    extensionloc = strfind(ciftifile,'.dtseries.nii');
    if isempty(extensionloc)
        extensionloc = strfind(ciftifile,'.dscalar.nii');
    end
    
    outname = [ciftifile(1:(extensionloc-1)) '_' hemcaps{hem}(1) '.border'];
    
    warning off
    delete(outname)
    fid = fopen(outname,'at'); %open the output file for writing
    fclose(fid);
    
    header = {'<?xml version="1.0" encoding="UTF-8"?>',...
        '<BorderFile Version="1">',...
        '   <MetaData>',...
        '      <MD>',...
        '         <Name><![CDATA[UniqueID]]></Name>',...
        '         <Value><![CDATA[{f0f2bfce-760e-4474-b8e4-2a3b1e963eff}]]></Value>',...
        '      </MD>',...
        '   </MetaData>',...
        '   <BorderClassColorTable>',...
        '      <LabelTable>',...
        ['         <Label Key="1" Red="0" Green="0" Blue="0" Alpha="1"><![CDATA[' ciftifile(1:(extensionloc-1)) ']]></Label>'],...
        '      </LabelTable>',...
        '   </BorderClassColorTable>'};
    for i = 1:length(header)
        dlmwrite(outname,header{i},'-append','delimiter','');
    end
    
    
    header = {'   <BorderNameColorTable>',...
        '      <LabelTable>'};
    for col = 1:size(data,2)
        header{end+1} = ['         <Label Key="' num2str(col) '" Red="' num2str(thishem_bordercolors(col,1)) '" Green="' num2str(thishem_bordercolors(col,2)) '" Blue="' num2str(thishem_bordercolors(col,3)) '" Alpha="1"><![CDATA[' colnames{col} ']]></Label>'];
    end
    header(end+1 : end+2) = {'      </LabelTable>',...
        '   </BorderNameColorTable>'};
    
    for i = 1:length(header)
        dlmwrite(outname,header{i},'-append','delimiter','');
    end
    
    
    for col = 1:size(data,2)
        stufftowritehere = {'   <Border>',...
            ['      <Name>' colnames{col} '</Name>'],...
            ['      <ClassName>' ciftifile(1:(extensionloc-1)) '</ClassName>'],...
            '      <ColorName>CLASS</ColorName>'};
        for i = 1:length(stufftowritehere)
            dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
        end
        
        
        
        verts = find(data(:,col));
        for vert = verts(:)'
            
            thirdneighbors = intersect(neighbors(vert,3:end),neighbors(neighbors(vert,2),2:end));
            
            stufftowritehere = {'       <SurfaceProjectedItem>',...
                ['          <Structure>CORTEX_' hemcaps{hem} '</Structure>'],...
                ['          <StereotaxicXYZ>' num2str(MNI{hem}(vert,:)) '</StereotaxicXYZ>'],...
                '           <ProjectionBarycentric>',...
                '              <TriangleAreas>1 0 0</TriangleAreas>',...
                ['              <TriangleNodes>' num2str([vert-1 neighbors(vert,2)-1 thirdneighbors(1)-1]) '</TriangleNodes>'],...
                '              <SignedDistanceAboveSurface>0</SignedDistanceAboveSurface>',...
                '           </ProjectionBarycentric>',...
                '       </SurfaceProjectedItem>'};
            
            for i = 1:length(stufftowritehere)
                dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
            end
        end
        dlmwrite(outname,'   </Border>','-append','delimiter','');
    end
    dlmwrite(outname,'</BorderFile>','-append','delimiter','');
    
end