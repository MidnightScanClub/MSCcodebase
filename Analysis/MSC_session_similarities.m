
warning off

run_new = 1;

MSCnums = 1:10;

groupwidth = 1;

outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/corrmats/'];
mkdir(outfolder);


parcels_LR = ['/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii'];
parcels_struct = ft_read_cifti_mod(parcels_LR);
parcels = parcels_struct.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

all_corrmats = zeros(nnz(triu(true(length(parcelIDs)),1)),0);

transitions = [];
colors = [1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.85 .85 .85;.5 .5 .3;.8 .35 .5;.5 .75 .2];
sess_sub = [];
legendnames = cell(0,1);
subcounter = 0;



for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    this_outfolder = [outfolder '/' MSCname];
    mkdir(this_outfolder)
    
    tmaskfile = ['/data/nil-bluearc/GMT/Evan/MSC/subjects/' MSCname '_TMASKLIST.txt'];
    [subjectlist, tmask_list] = textread(tmaskfile,'%s %s');
    subcounter = subcounter + 1;
    for s = 1 : length(subjectlist)
        
        this_corrmat = smartload([this_outfolder '/corrmat_sess' num2str(s) '.mat']);
        
        all_corrmats(:,end+1) = this_corrmat(triu(true(size(this_corrmat)),1));
        
        sess_sub(end+1) = subcounter;
        
        
    end
    legendnames{end+1} = MSCname;
    
    
    if floor(MSCnum) == (length(MSCnums) ./ 2)
        
    end
        
    
    
end

subcounter = subcounter + 1;
avg_index = subcounter;
group_corrmat = smartload('/data/nil-bluearc/GMT/Evan/MSC/corrmats/MSCavg_corrmat.mat');
for s = 1:groupwidth
    all_corrmats(:,end+1) = group_corrmat(triu(true(size(group_corrmat)),1));
    sess_sub(end+1) = subcounter;
end
legendnames{end+1} = ['MSCavg'];


corrmat_similarities = paircorr_mod(all_corrmats);
corrmat_similarities(sess_sub==avg_index,sess_sub==avg_index) = 0;
similarity_figmaker(corrmat_similarities,sess_sub',[.2 .9],legendnames,'Session Similarity')

