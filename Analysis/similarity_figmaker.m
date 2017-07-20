function parcel_correlmat_figmaker_alt(corrmat,assignmentsfile,limits,networklabels,titlename)
%parcel_correlmat_figmaker(corrmat,assignmentsfile,limits,titlename)

%networklabels = {'Unassigned','Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip','Unknown','Unknown','Unknown'};
%colors = [1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.85 .85 .85;.5 .5 .3;.8 .35 .5;.5 .75 .2];
colors = repmat([0 0 0],100,1);%[1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.85 .85 .85;.5 .5 .3;.8 .35 .5;.5 .75 .2];

nodes = size(corrmat,1);

yaxis_label_position_params = [1.4 (nodes/100) * -1.2]; %horiz vert
xaxis_label_position_params = [(nodes/100) (nodes/100)*1.75 + .6]; %horiz vert

corrmat(logical(diag(ones(size(corrmat,1),1)))) = 0;

if ischar(assignmentsfile)
    assignments = load(assignmentsfile);
else
    assignments = assignmentsfile;
end
assignments(assignments<0) = 0;
IDs = unique(assignments);
networklabels = networklabels(IDs);
colors = colors(IDs,:);

[communities sorti] = sort(assignments);

sorted_corrmat = corrmat(sorti,sorti);

transitions = find(communities(1:end-1) ~= communities(2:end));
transitions_plusends = [1 transitions(:)' length(communities)];
centers = transitions_plusends(1:end-1) + ((transitions_plusends(2:end) - transitions_plusends(1:end-1))/2);

figure;
imagesc(sorted_corrmat,limits)

% colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
colormap jet
colorbar

set(gca,'XTicklabel','')
set(gca,'YTicklabel','') 

set(gca,'Xtick',centers)
tx= text(centers+ xaxis_label_position_params(1),ones(1,length(centers))*(length(corrmat)) + xaxis_label_position_params(2),networklabels);
set(tx,'HorizontalAlignment','right','VerticalAlignment','bottom','Rotation',45)

for i = 1:length(tx)
     %set(tx_outline(i),'Color','k','FontSize',15,'FontWeight','bold'); 
     set(tx(i),'Color',[colors(i,1) colors(i,2) colors(i,3)],'FontName','Helvetica','FontSize',14,'FontWeight','bold');   
end


set(gca,'Ytick',centers)
ty= text(-1*ones(1,length(centers)) + yaxis_label_position_params(1),centers+yaxis_label_position_params(2),networklabels);
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

for i = 1:length(ty)
    set(ty(i),'Color',[colors(i,1) colors(i,2) colors(i,3)],'FontName','Helvetica','FontSize',14,'FontWeight','bold');   
end



%set(gca,'XtickLabel',networklabels)
%xticklabel_rotate([],60,[],'Fontsize',12);

%set(gca,'YtickLabel',networklabels)
set(gca,'FontSize',20)

set(gcf,'Position',[813 30 1102 805])



title(titlename);
if ~isempty(transitions)
    hline(transitions+.5,'k')
    vline(transitions+.5,'k')
end

set(gcf,'Color',[1 1 1])