function make_block_diagram(recolored_allcolumns_file,thresholds)

colors = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;.95 .95 .75;0 .4 0;.5 .5 .5];

reorder = [1 16 7 15 8 13 14 3 5 10 17 6 11 12 9 4 2 18];

if ischar(thresholds)
    thresholds = load(thresholds);
end

reorderedcolors = colors(reorder,:);

origmatrix = ft_read_cifti_mod(recolored_allcolumns_file); origmatrix = origmatrix.data;
origmatrix(origmatrix<1) = 18; origmatrix(origmatrix>17) = 18;

matrix = ones(size(origmatrix)) .* 18;
for i = 1:17
    matrix(origmatrix==i) = find(reorder==i);
end

matrix = matrix(:,[end:-1:1]);
thresholds = thresholds([end:-1:1]);

[~,modefreq]=mode(matrix,1);
allmerged_inds = find(modefreq==size(matrix,1));
if ~isempty(allmerged_inds)
matrix = matrix(:,allmerged_inds(end):end);
thresholds = thresholds(allmerged_inds(end):end);
end

sortedmat = sortrows(matrix,[1:size(matrix,2)]);

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
imagesc(sortedmat); colormap(reorderedcolors);

set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
xticks = get(gca,'XTick');
for i = 1:length(xticks)
    string = sprintf('%.3f',thresholds(xticks(i)));
    threshlabel(i,:) = string(2:end);
end
set(gca,'YTick',[])
set(gca,'XTickLabel',threshlabel)

dotsloc = strfind(recolored_allcolumns_file,'.');
outname = [recolored_allcolumns_file(1:dotsloc(1)-1) '_block_diagram'];
try
    export_fig(gca,[outname '.pdf'],'-nocrop')
catch
    savefig(gcf,[outname '.fig'])
end

