function spring_embedding_func_easy_crossthresh(corrmat,assignments,Kc,L0,distances,xdist,thresholds,outname)
% Spring-embedding algorithm
%[h nodepos] = spring_embedding_func_easy_crossthresh(corrmat,assignments,Kc,L0,distances,xdist,kden_thresholds,output_figurename)
%Kc: try 1
%L0:try 25




if ischar(corrmat)
    corrmat = smartload(corrmat);
end
if ischar(assignments)
    assignments = load(assignments);
end

if exist('distances')
    if ischar(distances)
        distances = smartload(distances);
    end
    corrmat(distances < xdist) = 0;
end

power_surf_colormap = [.67 .67 .67;.67 .67 .67;1 0 0;0 0 .6;1 1 0;1 .7 .4;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .2 .4;.6 1 .6;0 0 1;1 1 .8;0 .4 0;repmat([.25 .25 .25],100,1)];


colors = [];
assignments = assignments + 2;
for i = 1:length(assignments)
    if isinteger(assignments(i))
        colors(i,:) = power_surf_colormap(assignments(i),:);
    else
        colors(i,:) = mean([power_surf_colormap(floor(assignments(i)),:);power_surf_colormap(ceil(assignments(i)),:)],1);
    end
end

%graph = abs(graph);
%mods = {1:28;29:58;59:63;64:77;78:90;91:148;149:153;154:184;185:209;210:227;228:240;241:249;250:260;261:264};  % Consensus foci

%mods = {1:28;29:58;59:63;64:81;82:94;95:152;153:157;158:188;189:213;214:231;232:244;245:253;254:264};  % Consensus foci
%mods = {1:17;18:49;50:57;58:68;69:97;98:122;123:126;127:137;138:144;145:203;204:208;209:228;229:264};
nodes = size(corrmat,1);



%Define circle
%r = 25;
r = L0;
theta = 0:2*pi/nodes:2*pi-2*pi/nodes;
nodepos = [r*sin(theta); r*cos(theta)]';

% Length of display area
%L0 = 10;


%
% for m = 1:nodes
%     xm = repmat(nodepos(m,1),nodes-1,1);
%     ym = repmat(nodepos(m,2),nodes-1,1);
%     xother = nodepos(:,1);
%     xother(m) = [];
%     yother = nodepos(:,2);
%     yother(m) = [];
%     Im = I(m,:)';
%     Im(m) = [];
%     Km = K(m,:)';
%     Km(m) = [];
%
%     dEdxm = sum(Km.*(xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5));
%     dEdym = sum(Km.*(ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5));
%     deltam(m) = sqrt((dEdxm).^2+(dEdym).^2);
% end


xm = repmat(nodepos(:,1),1,nodes);
xm(logical(eye(nodes))) = 0;
ym = repmat(nodepos(:,2),1,nodes);
ym(logical(eye(nodes))) = 0;
xother = repmat(nodepos(:,1),1,nodes)';
xother(logical(eye(nodes))) = 0;
yother = repmat(nodepos(:,2),1,nodes)';
yother(logical(eye(nodes))) = 0;


for kden = thresholds
    
    [graph,~,~] = matrix_thresholder(corrmat,kden,'kden');
    degree = zeros(nodes,1);
    for n = 1:nodes
        row = zeros(1,nodes);
        row(n) = 1;
        degree(n) = sum(graph(n,~row));
    end
    
    [~,connectedcomps] = conncomp(graph);
    compIDs = unique(connectedcomps);
    biggestcomp = [];
    biggestcompsize = 0;
    for i = 1:length(compIDs)
        if nnz(connectedcomps==compIDs(i)) > biggestcompsize
            biggestcompsize = nnz(connectedcomps==compIDs(i));
            biggestcomp = compIDs(i);
        end
    end
    inbiggestcomp = (connectedcomps==biggestcomp)';
    
    outnodevec = (degree==0);
    outnodeind = find(degree==0);
    % Find shortest path using dijkstra algorithm
    D = dijkstra(graph,graph);
    D(isinf(D))=0;
    %D = graph;
    
    % Find desirable length between nodes
    maxd = max(D(:));
    L = L0/maxd;
    I = L.*D;
    
    % Strength of springs between nodes
    %Kc = 1;
    K = Kc./(D.^2);
    K(isinf(K))=0;
    Im = I;
    Km = K;
    
    dEdxmstep = Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5));
    dEdxmstep(isnan(dEdxmstep))=0;
    dEdxm = sum(dEdxmstep,2);
    dEdymstep = Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5));
    dEdymstep(isnan(dEdymstep))=0;
    dEdym = sum(dEdymstep,2);
    deltam = sqrt((dEdxm).^2+(dEdym).^2)';
    
    count = 1;
    while max(deltam)>.05
        string{count} = ['Iteration ' num2str(count)];
        if count==1; fprintf('%s',string{count}); else fprintf([repmat('\b',1,length(string{count-1})) '%s'],string{count}); end
        if count>20000
            break
        end
        [maxval maxi] = max(deltam);
        xm = repmat(nodepos(maxi,1),nodes-1,1);
        ym = repmat(nodepos(maxi,2),nodes-1,1);
        xother = nodepos(:,1);
        xother(maxi) = [];
        yother = nodepos(:,2);
        yother(maxi) = [];
        Im = I(maxi,:)';
        Im(maxi) = [];
        Km = K(maxi,:)';
        Km(maxi) = [];
        dEdxm = sum(Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
        dEdym = sum(Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
        
        d2Edxm2 = sum(Km.*(1-((Im.*(ym-yother).^2)./((xm-xother).^2+(ym-yother).^2).^(3/2))));
        d2Edym2 = sum(Km.*(1-((Im.*(xm-xother).^2)./((xm-xother).^2+(ym-yother).^2).^(3/2))));
        
        d2Edxmdym = sum(Km.*((Im.*(xm-xother).*(ym-yother))./((xm-xother).^2+(ym-yother).^2).^(3/2)));
        d2Edymdxm = sum(Km.*((Im.*(xm-xother).*(ym-yother))./((xm-xother).^2+(ym-yother).^2).^(3/2)));
        
        A = [d2Edxm2 d2Edxmdym; d2Edymdxm d2Edym2];
        B = [-dEdxm;-dEdym];
        
        X = linsolve(A,B);
        deltax = X(1);
        deltay = X(2);
        nodepos(maxi,1) = nodepos(maxi,1)+deltax;
        nodepos(maxi,2) = nodepos(maxi,2)+deltay;
        
        %     for m = 1:nodes
        %         xm = repmat(nodepos(m,1),nodes-1,1);
        %         ym = repmat(nodepos(m,2),nodes-1,1);
        %         xother = nodepos(:,1);
        %         xother(m) = [];
        %         yother = nodepos(:,2);
        %         yother(m) = [];
        %         Im = I(m,:)';
        %         Im(m) = [];
        %         Km = K(m,:)';
        %         Km(m) = [];
        %
        %         dEdxm = sum(Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
        %         dEdym = sum(Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
        %         deltam(m) = sqrt((dEdxm).^2+(dEdym).^2);
        %     end
        
        xm = repmat(nodepos(:,1),1,nodes);
        xm(logical(eye(nodes))) = 0;
        ym = repmat(nodepos(:,2),1,nodes);
        ym(logical(eye(nodes))) = 0;
        xother = repmat(nodepos(:,1),1,nodes)';
        xother(logical(eye(nodes))) = 0;
        yother = repmat(nodepos(:,2),1,nodes)';
        yother(logical(eye(nodes))) = 0;
        
        Im = I;
        Km = K;
        
        dEdxmstep = Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5));
        dEdxmstep(isnan(dEdxmstep))=0;
        dEdxm = sum(dEdxmstep,2);
        dEdymstep = Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5));
        dEdymstep(isnan(dEdymstep))=0;
        dEdym = sum(dEdymstep,2);
        deltam = sqrt((dEdxm).^2+(dEdym).^2)';
        
        %     h=figure('position',[1500 500 700 500]);
        %     plot(nodepos(:,1),nodepos(:,2),'r.','MarkerSize',30)
        %
        %     plot(nodepos(:,1),nodepos(:,2),'b.','MarkerSize',30)
        %
        %     for i = 1:size(graph,1)
        %         initcoordmat = repmat(nodepos(i,:),nnz(graph(i,:)),1);
        %         line([initcoordmat(:,1) nodepos(logical(graph(i,:)),1)]',[initcoordmat(:,2) nodepos(logical(graph(i,:)),2)]','Color','r');
        %     end
        count = count+1;
    end
    disp(' ')
    
    nodepos_forfig = nodepos;
    
    h=figure('position',[1500 500 1000 750],'Color','white');
    hold
    axis([-r r -r r])
    
    %Paste nodes not in network at bottom
    outnodepos = -r:(2*r)/(length(outnodeind)-1):r;
    for j = 1:length(outnodeind)
        nodepos_forfig(outnodeind(j),:) = [outnodepos(j) -100];%[outnodepos(j) -r];
    end
    
    %Draw lines between nodes
    for i = 1:size(graph,1)
        initcoordmat = repmat(nodepos_forfig(i,:),nnz(graph(i,:)),1);
        line([initcoordmat(:,1) nodepos_forfig(logical(graph(i,:)),1)]',[initcoordmat(:,2) nodepos_forfig(logical(graph(i,:)),2)]','Color',[155/255 155/255 155/255]);
    end
    
    %Mark and color node positions
    int=1;
    for i = 1:nodes
        plot(nodepos_forfig(int,1),nodepos_forfig(int,2),'MarkerFaceColor',[colors(i,1) colors(i,2) colors(i,3)],'Marker','o','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1);
        int = int+1;
    end
    
    highdegree_maincomponentnodes = logical((degree >= 3) .* inbiggestcomp);
    
    minx_forzoom = min(nodepos_forfig(highdegree_maincomponentnodes,1)) - std(nodepos_forfig(highdegree_maincomponentnodes,1));
    maxx_forzoom = max(nodepos_forfig(highdegree_maincomponentnodes,1)) + std(nodepos_forfig(highdegree_maincomponentnodes,1));
    miny_forzoom = min(nodepos_forfig(highdegree_maincomponentnodes,2)) - std(nodepos_forfig(highdegree_maincomponentnodes,2));
    maxy_forzoom = max(nodepos_forfig(highdegree_maincomponentnodes,2)) + std(nodepos_forfig(highdegree_maincomponentnodes,2));
    axis(gca,[minx_forzoom maxx_forzoom miny_forzoom maxy_forzoom])
    axis off
    if exist('outname')
        saveas(gcf,[outname '_' num2str(kden) '.fig'])
    end
    
    
%     % Simple fig
%     
%     h=figure('position',[1500 500 1000 750],'Color','white');
%     hold
%     axis([-r r -r r])
%     
%     assignments_forfig = assignments;
%     
%     assignments_forfig(outnodevec) = 0;
%     assignments_forfig(~inbiggestcomp) = 0;
%     assignmentIDs = unique(assignments_forfig); assignmentIDs(assignmentIDs<3) = []; assignmentIDs(assignmentIDs > 19) = [];
%     centroidpos = zeros(length(assignmentIDs),2);
%     
%     for IDnum = 1:length(assignmentIDs);
%         ID = assignmentIDs(IDnum);
%         centroidpos(IDnum,:) = mean(nodepos_forfig(assignments_forfig==ID,:),1);
%     end
%     
%     networkgraph = zeros(length(assignmentIDs));
%     for i = 1:size(networkgraph,1)
%         for j = (i+1):size(networkgraph,1)
%             connectionstrength(i,j) = nnz(graph(assignments_forfig==assignmentIDs(i),assignments_forfig==assignmentIDs(j)));
%         end
%     end
%     
%     for i = 1:size(networkgraph,1)
%         for j = (i+1):size(networkgraph,1)
%             connectionratio = connectionstrength(i,j) ./ (max(max(connectionstrength)));
%             if connectionratio > 0
%                 line([centroidpos(i,1) centroidpos(j,1)],[centroidpos(i,2) centroidpos(j,2)],'Color',[.85-(connectionratio*.35) .85-(connectionratio*.35) .85-(connectionratio*.35)],'LineWidth',(connectionratio*10 + 2));
%                 %line([centroidpos(i,1) centroidpos(j,1)],[centroidpos(i,2) centroidpos(j,2)],'Color',[155/255 155/255 155/255],'LineWidth',(connectionratio*10 + 2));
%             end
%             
%         end
%     end
%     
%     for IDnum = 1:length(assignmentIDs);
%         ID = assignmentIDs(IDnum);
%         if isinteger(ID)
%             networkcolors = power_surf_colormap(ID,:);
%         else
%             networkcolors = mean([power_surf_colormap(floor(ID),:);power_surf_colormap(ceil(ID),:)],1);
%         end
%         networksize = 300 * nnz(assignments_forfig==ID) / nnz(assignments_forfig) + 20;
%         plot(centroidpos(IDnum,1),centroidpos(IDnum,2),'MarkerFaceColor',[networkcolors(1) networkcolors(2) networkcolors(3)],'Marker','o','MarkerSize',networksize,'MarkerEdgeColor','k','LineWidth',3);
%     end
%     
%     axis(gca,[minx_forzoom maxx_forzoom miny_forzoom maxy_forzoom])
%     
%     if exist('outname')
%         saveas(gcf,[outname '_' num2str(kden) '_simple.fig'])
%     end
    
end

%
%     % For module coloring
%      int=1;
%      for m = 1:length(mods)
%          for i = 1:length(mods{m})
%              plot(nodepos(int,1),nodepos(int,2),'MarkerFaceColor',[colors(m,1) colors(m,2) colors(m,3)],'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',2);
%              int = int+1;
%          end
%      end
%     for i = 1:length(nodepos)
%         text(nodepos(i,1),nodepos(i,2),[num2str(seeds(i,1)) '/' num2str(seeds(i,2)) '/' num2str(seeds(i,3))]);
%     end