MSCnums = 1:10;
edgethresh = .5;
xdist = 30;
thresholds = [.003 .004 .005:.005:.05];


%% Analyses to run
run_vertexwise_infomap = 1;
run_parcellation = 1;
make_parcel_corrmats = 1;
make_parcel_distmats = 1;
run_parcel_infomap = 1;
run_spring_embedding = 1;




for MSCnum = MSCnums
    
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    %%
    infomap_outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/infomap/' MSCname '_infomap_p003_p005_p05/'];
    parcellation_outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/parcels/'];
    springembed_outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/spring_embed/'];
    parcelinfomap_outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/parcels/' MSCname '_parcels_LR_infomap_p003_p05/'];
    surfdir = ['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/'];
    dmatname = [infomap_outfolder '/distmat_surf_geodesic_vol_euc_xhem_large_uint8.mat'];
        
    
    
    tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
    [subjectlist, tmask_list] = textread(tmaskfile,'%s %s');
    ciftifiles = cell(length(subjects),1);
    for s = 1:length(subjectlist)
        ciftifiles{s} = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_smallwall_native_freesurf/' subjectlist{s} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    end
        
    
    ciftidata = ['/data/nil-bluearc/GMT/Evan/MSC/subjects/' MSCname '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    parcellation_file = [parcellation_outfolder '/' MSCname '_parcels_LR.dtseries.nii'];
    
    mkdir(infomap_outfolder);
    mkdir(springembed_outfolder)
    
    %% vertexwise infomap
    if run_vertexwise_infomap
        
        distances = smartload(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/normalwall_distmat_333_native_freesurf/distmat_surf_geodesic_vol_euc.mat']);
        distances = uint8(distances);
        distances(1:29696,29697:59412) = 255;
        distances(29697:59412,1:29696) = 255;
        save(dmatname,'distances','-v7.3');
        clear distances
        
        cd(infomap_outfolder);
        
        for s = 1:length(subjectlist)
            tmask = load(tmask_list{s}); 
            data = ft_read_cifti_mod(ciftifiles{s});
            if s==1
                alldata = data.data(:,logical(tmask));
            else
                alldata = [alldata data.data(:,logical(tmask))];
            end
            data.data = [];
        end
        
        corrmat = paircorr_mod(alldata');
        clear alldata
        corrmat(isnan(corrmat)) = 0;
        corrmat = FisherTransform(corrmat);
        
        structure_indices = data.brainstructure;
        structure_indices = structure_indices(structure_indices>0);
        structure_indices = (structure_indices > 2) +1;
        
        dmatname = [infomap_outfolder '/distmat_surf_geodesic_vol_euc_xhem_large_uint8.mat'];
        
        Run_Infomap_2015(corrmat, dmatname, xdist, thresholds, 0, infomap_outfolder, 8, structure_indices);
        clear corrmat
        
        communities = modify_clrfile('simplify','rawassn.txt',400);
        regularized = regularize(communities);
        data.data = regularized;
        ft_write_cifti_mod([MSCname '_rawassn_minsize400_regularized'],data);
        try movefile('rawassn_minsize400_regularized.dtseries.nii',[MSCname '_rawassn_minsize400_regularized.dtseries.nii']); catch; end
        consensus_maker_knowncolors([MSCname '_rawassn_minsize400_regularized.dtseries.nii'],[],[],1)
        make_block_diagram([MSCname '_rawassn_minsize400_regularized_allcolumns_recolored.dtseries.nii'],thresholds)
        cifti_to_border_v2([MSCname '_rawassn_minsize400_regularized_recolored.dscalar.nii'],1,1,'default')
        
    end
    
    %% parcellation
    if run_parcellation
        
        mkdir(parcellation_outfolder)
        mkdir([parcellation_outfolder '/' MSCname '/'])
        cd([parcellation_outfolder '/' MSCname '/'])
        
        for s = 1:length(subjectlist)
            tmask = load(tmask_list{s}); 
            data = ft_read_cifti_mod(ciftifiles{s});
            if s==1
                alldata = data.data(:,logical(tmask));
            else
                alldata = [alldata data.data(:,logical(tmask))];
            end
            data.data = [];
        end
        
        corrmat = paircorr_mod(alldata');
        clear alldata
        
        corrmat(isnan(corrmat)) = 0;
        corrmat = FisherTransform(corrmat);
        
        data.data = corrmat;
        clear corrmat;
        
        surface_parcellation_singlesub(MSCname,data,surfdir,1,0,[parcellation_outfolder '/' MSCname '/'])
        
        parcel_creator_cifti('corrofcorr_allgrad_LR_subcort_smooth2.55_wateredge_avg.dtseries.nii',[MSCname '_parcels'],edgethresh,ciftidata)
        movefile([parcellation_outfolder '/' MSCname '/' MSCname '_parcels_edgethresh_' num2str(edgethresh) '.dtseries.nii'],parcellation_file);
        try delete([parcellation_outfolder '/corrofcorr_allgrad_LR_subcort_smooth2.55.dtseries.nii']); catch; end
        
    end
    
    %% parcel corrmats
    if make_parcel_corrmats
        
        parcels = ft_read_cifti_mod(parcellation_file);
        parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
        
        for s = 1:length(subjectlist)
            tmask = load(tmask_list{s}); 
            data = ft_read_cifti_mod(ciftifiles{s});
            if s==1
                alldata = data.data(:,logical(tmask));
            else
                alldata = [alldata data.data(:,logical(tmask))];
            end
            data.data = [];
        end
        
        tcs = zeros(size(alldata,2),length(parcelIDs));
        for IDnum = 1:length(parcelIDs)
            tcs(:,IDnum) = mean(alldata(parcels.data==parcelIDs(IDnum),:),1);
        end
        clear alldata
        corrmat = paircorr_mod(tcs);
        corrmat(isnan(corrmat)) = 0;
        corrmat = FisherTransform(corrmat);
        
        save([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat'],'corrmat');
        
    end
    
    
    %% parcel distmats
    if make_parcel_distmats
        
        distances = smartload(dmatname);
        
        parcels = ft_read_cifti_mod([parcellation_outfolder '/' MSCname '_parcels_LR.dtseries.nii']);
        parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
        
        parcel_centroids = zeros(length(parcelIDs),1);
        for IDnum = 1:length(parcelIDs)
            parcelinds = find(parcels.data==parcelIDs(IDnum));
            within_parcel_distances = distances(parcelinds,parcelinds);
            [~,centroidind] = min(sum(within_parcel_distances,2));
            parcel_centroids(IDnum) = parcelinds(centroidind);
        end
        parcel_distances = distances(parcel_centroids,parcel_centroids);
        save([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat'],'parcel_distances');
        clear distances        
        
    end
    
    %% parcel infomap
    if run_parcel_infomap
        
        mkdir(parcelinfomap_outfolder)
        cd(parcelinfomap_outfolder)
        corrmat = smartload([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat']);
        distances = smartload([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat']);
        
        
        Run_Infomap_nopar(corrmat, distances, xdist, thresholds, 0, parcelinfomap_outfolder);
        communities = modify_clrfile('simplify','rawassn.txt',4);
        communities = load('rawassn_minsize4.txt');
        regularized = regularize(communities);
        
        parcels = ft_read_cifti_mod(parcellation_file);
        parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
        
        regularized_out = parcels;
        regularized_out.data = zeros(size(parcels.data,1),size(regularized,2));
        for IDnum = 1:length(parcelIDs)
            regularized_out.data(parcels.data==parcelIDs(IDnum),:) = repmat(regularized(IDnum,:),nnz(parcels.data==parcelIDs(IDnum)),1);
        end
        ft_write_cifti_mod('rawassn_minsize4_regularized.dtseries.nii',regularized_out);
        consensus_maker_knowncolors('rawassn_minsize4_regularized.dtseries.nii',[],[infomap_outfolder '/' MSCname '_rawassn_minsize400_regularized_recolored.dscalar.nii'],1,[],[],parcellation_file)
        
    end
    
    
    
    %% spring embedding
    if run_spring_embedding
        
        cd(springembed_outfolder)
        corrmat = smartload([parcellation_outfolder '/' MSCname '_parcel_corrmat.mat']);
        distances = smartload([parcellation_outfolder '/' MSCname '_parcel_distances_xhemlarge.mat']);
        consensus = load([parcelinfomap_outfolder '/rawassn_minsize4_regularized_recolored.txt']);
        spring_embedding_func_easy_crossthresh(corrmat,consensus,1,25,distances,xdist,thresholds,[MSCname '_spring_embed']);
        close all
        
    end
    
    
    
    
end




    








