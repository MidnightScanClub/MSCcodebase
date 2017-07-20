file = '/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/convergence/similarity_metrics.mat';

load(file)
addeddata = load(file);
names = fieldnames(addeddata);
clear addeddata

for j = 1:length(names)
        if ~strcmp(names{j},'PC_all')
            eval([names{j} '(' names{j} '==0) = NaN;']);
        end
end


%%
MSCnums = 1:10;
datalength_totest = [2.5 5 10 : 10 : 100];
outfolder = pwd;
colors = [0 0 0; .9 .9 0; 0 1 0 ; 1 0 0; 0 0 1; .2 1 1; 1 0 1; .7 .7 .7; 0 .6 .6 ; 1 .5 0];

h = figure('Color','white','position',[1982 478 1352 804],'DefaultAxesFontSize',40);
hold on

final_similarities_toplot = zeros(length(MSCnums),length(datalength_totest));

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(corrmat_similarity(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(corrmat_similarity(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(corrmat_similarity,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(MSCnum,:),'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3,'FontName','Helvetica')
title(gca,['Connectivity Matrix'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50,'FontName','Helvetica')
ylabel('Correlation to other half','Fontweight','bold','FontSize',50,'FontName','Helvetica')
ylim([.1 .95])
legend(legendnames,'Location','SouthEast')

try export_fig(gca,[outfolder '/Correlation Similarity.pdf'])
catch
    savefig(gcf,[outfolder '/Correlation Similarity.fig'])
end




h = figure('Color','white','position',[1982 478 1352 804]);
hold on



for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(community_dice(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(community_dice(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(community_dice,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(MSCnum,:),'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Network Assignment'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('Dice to other half','Fontweight','bold','FontSize',50)
ylim([.1 .95])
legend(legendnames,'Location','SouthEast')

try export_fig(gca,[outfolder '/Community Detection Overlap.pdf'])
catch
    savefig(gcf,[outfolder '/Community Detection Overlap.fig'])
end



h = figure('Color','white','position',[1982 478 1352 804]);
hold on

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(PC_similarity(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(PC_similarity(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(PC_similarity,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(MSCnum,:),'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Participation Coefficient'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('Correlation to other half','Fontweight','bold','FontSize',50)
ylim([.1 .95])
legend(legendnames,'Location','SouthEast')

try export_fig(gca,[outfolder '/PC Similarity.pdf'])
catch
    savefig(gcf,[outfolder '/PC Similarity.fig'])
end




h = figure('Color','white','position',[1982 478 1352 804]);
hold on

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(GEff_delta(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(GEff_delta(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(GEff_delta,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,meancorrmat.*100,'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Global Efficiency'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('% Difference from other half','Fontweight','bold','FontSize',50)
ylim([0 20])
legend(legendnames,'Location','NorthEast')

try export_fig(gca,[outfolder '/Global Efficiency Difference.pdf'])
catch
    savefig(gcf,[outfolder '/Global Efficiency Difference.fig'])
end




h = figure('Color','white','position',[1982 478 1352 804]);
hold on

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(Modularity_delta(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(Modularity_delta(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(Modularity_delta,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,meancorrmat.*100,'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Modularity'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('% Difference from other half','Fontweight','bold','FontSize',50)
ylim([0 20])
legend(legendnames,'Location','NorthEast')

try export_fig(gca,[outfolder '/Modularity Difference.pdf'])
catch
    savefig(gcf,[outfolder '/Modularity Difference.fig'])
end




h = figure('Color','white','position',[1982 478 1352 804]);
hold on

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(GEff_all(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(GEff_all(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(GEff_all,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(MSCnum,:),'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Global Efficiency'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('Global Efficiency value','Fontweight','bold','FontSize',50)
%ylim([0 20])
legend(legendnames,'Location','NorthEast')

try export_fig(gca,[outfolder '/Global Efficiency value.pdf'])
catch
    savefig(gcf,[outfolder '/Global Efficiency value.fig'])
end




h = figure('Color','white','position',[1982 478 1352 804]);
hold on

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    legendnames{MSCnum} = MSCname;
    meancorrmat = nanmean(Modularity_all(:,:,MSCnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(Modularity_all(:,:,MSCnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(Modularity_all,0,1);
    final_similarities_toplot(MSCnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(MSCnum,:),'Color',colors(MSCnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3)
title(gca,['Modularity'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50)
ylabel('Modularity Value','Fontweight','bold','FontSize',50)
%ylim([0 20])
legend(legendnames,'Location','NorthEast')

try export_fig(gca,[outfolder '/Modularity value.pdf'])
catch
    savefig(gcf,[outfolder '/Modularity value.fig'])
end