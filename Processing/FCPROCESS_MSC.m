function [allends switches] = FCPROCESS_MSC(datafile,tdir,tmasktype,cleanup,roilist,varargin)

% jdp 2/22/2012
%
% IMPORTANT: fcp_process removes the /newtargetdir/vcdir/ folder when it
% begins processing a subject, so DO NOT set the targetdir to the directory
% where your 333 BOLD data exist - use a new directory, specific to FC
% data. Since fcp_initial_process is a wrapper for fcp_process, this goes
% for that script too.
%
% This script is the fcprocessing script for the Petersen lab.
%
% The datalist specifies the data to process, and is a 5-column file:
% srcdir vcnum prmfile TR skipframes
% e.g.:
% /data/src/ vc11111 /data/prms/vc11111.prm 2.5 5
% /data/src/ vc22222 /data/prms/vc22222.prm 2.2 5
%
% We presume the file structure for fMRI data:
% /srcdir
%   /vcnum
%     /boldX
%     /boldY
%     /atlas
%
% Where the bold# for resting runs are specified in the prmfile, which
% contains the single line:
% set boldruns = (X Y whatever the bold run numbers are)
% e.g: set boldruns = (2 4 6)
%
% The srcdir and prmfile need to be hard paths
%
% The TR is obvious, and the skipframes are the number of frames to skip at
% the beginning of each run (5 is generous but common).
%
% targetdir: where to write files to
% tmasktype can be one of two things:
%     'ones' - just skips the skipframes at the start of each run
%     tmasklist - a tmasklist from COHORTSELECOR.m
% cleanup: delete all files but QC files and final fcimage (0=no. 1=yes)
% freesurferlist: a file with 2 columns, the first the subject ID, the second 
% 	where the freesurfer segmentation exists that was used for MAKE_FS_MASKS.csh
% e.g.,
% 	vc11111	/data/cn4/segmentation/freesurfer5_supercomputer
%	vc22222	/data/cn4/segmentation/freesurfer5_supercomputer
%
% The processing order is:
%    demean/detrend (mask)
%    extract nuisance signals
%    multiple regression (mask)
%    *interpolate* (mask)
%    temporal filter (butter1 filtfilt low-pass)
%    demean/detrend (mask)
%    spatial blur (gauss_4dfp)
%
%
%
% IMPORTANT: FCPROCESS removes the /newtargetdir/vcdir/ folder when it
% begins processing a subject, so DO NOT set the newtargetdir to the directory
% where your 333 BOLD data exist - use a new directory, specific to FC
% data.
%
% FCPROCESS(datalist,targetdir,tmasktype,cleanup,roilist,freesurferlist)
% FCPROCESS('datalist.txt','/put/data/here/','ones',1,'bigbrain264.txt','fsinfo.txt');
% FCPROCESS('datalist.txt','/put/data/here/','specificmasks.txt',1,'fab5roiz.txt','fsinfo.txt');






%% IMPORTANT VARIABLES
avidir = '/data/cninds01/data2/nil-tools';
releasedir = '/data/petsun4/data1/solaris';
roidir='/home/usr/fidl/lib/';
startdir=pwd;
FDoutlier=0.25; % flag frames over this
GLMoutlierlim=2; % flag frames over this z-score of global signal
GREYoutlierlim=2;
voxdim='333';
switches.WMero=4;
switches.CSFero=1;
set(0, 'DefaultFigureVisible', 'off'); 

% [rois x y z] = textread(roilist,'%s%f%f%f');
% roimasks = zeros(147456,size(rois,1));
% %roimasks = zeros(147456,size(rois,1));
% for i=1:size(rois,1)
%     i
%     roimasks_orig(:,i)=read_4dfpimg_HCP(rois{i,1});
% end
[rois] = textread(roilist,'%s');
roimasks_orig = read_4dfpimg_HCP(rois{1});

switches.doregression=0;
switches.regressiontype=1;
switches.motionestimates=1;
switches.WM=0;
switches.V=0;
switches.GS=0;
switches.dointerpolate=0;
switches.dobandpass=0;
switches.temporalfiltertype=3;
switches.lopasscutoff=0;
switches.hipasscutoff=0;
switches.order=0;
switches.doblur=0;
switches.blurkernel=0;


%% CHECK THAT BASIC FILES ARE ACTUALLY THERE

% read in the subject data including location, vcnum, and boldruns
[df.filepath df.vcnum df.prmfile df.TR df.TRskip] = textread(datafile,'%s%s%s%f%d');
fprintf('CHECKING THAT DATA EXIST\n');
pause(2);
numsubs=size(df.vcnum,1);

for i=1:numsubs
    QC(i).vcnum=df.vcnum{i,1};
    QC(i).TR=df.TR(i,1);
    QC(i).TRskip=df.TRskip(i,1);
end

%% CHECK TEMPORAL MASKS
fprintf('CHECKING TEMPORAL MASKS\n');
switch tmasktype
    case 'ones'
    otherwise
        [tm.vcnum tm.tmaskfiles]=textread(tmasktype,'%s%s');
        if ~isequal(tm.vcnum,df.vcnum)
            error('masklist vcnumbers do not match datalist vcnumbers');
        end
end


%% READ IN FREESURFER
if ~isempty(varargin)
    [FSfile.vc FSfile.dir]=textread(varargin{1,1},'%s%s');
end
errorout=0;
for i=1:numsubs
    found=1;
    %for j=1:length(FSfile.dir)
        %if isequal(FSfile.vc{j,1},QC(i).vcnum) || isequal([FSfile.vc{j,1} '_2'],[QC(i).vcnum])
            QC(i).fsdir=[FSfile.dir{1} '/' FSfile.vc{1}];
          %  found=found+1;
       % end
   % end
%     if found==1
%     elseif found>1
%         errorout=1;
%         fprintf('%s has multiple matches to FSdirs\n',QC(i).vcnum);
%     elseif found<1
%         errorout=1;
%         fprintf('%s has zero matches to FSdirs\n',QC(i).vcnum);
%     end
end
% if errorout
%     error();
% end


%% SET SWITCHES

needcorrectinput=1;
while needcorrectinput
    switches.doregression=input('Do you want to regress nuisance signals? (1=y; 0=no): ');
    switch switches.doregression
        case 1
            needcorrectinput=0;
            
            needcorrectinput2=1;
            while needcorrectinput2
                switches.regressiontype=input('Regression: 1 - Freesurfer seeds: ');
                switch switches.regressiontype
                    
                    % classic or freesurfer seeds
                    case {1}
                        needcorrectinput2=0;
                        
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.motionestimates=input('Do you want to regress motion estimates and derivatives? (0=no; 1=R,R`; 2=FRISTON; 20=R,R`,12rand; 3:volt3): ');
                            switch switches.motionestimates
                                case {0,1,2,20,3}
                                    needcorrectinput1=0;
                            end
                        end
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.WM=input('Do you want to regress white matter signals and derivatives? (1=y; 0=no): ');
                            switch switches.WM
                                case {0,1}
                                    needcorrectinput1=0;
                            end
                        end
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.V=input('Do you want to regress ventricular signals and derivatives? (1=y; 0=no): ');
                            switch switches.V
                                case {0,1}
                                    needcorrectinput1=0;
                            end
                        end
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.GS=input('Do you want to regress global signal and derivative? (1=y; 0=no): ');
                            switch switches.GS
                                case {0,1}
                                    needcorrectinput1=0;
                            end
                        end
                        
                        % user seeds
                    case 2
                        needcorrectinput2=0;
                        switches.nus4dfplistfile=input('Enter the hard path to the list of nuisance regressor 4dfps ','s');
                        [switches.nus4dfpvcnum switches.nus4dfplist] = textread(switches.nus4dfplistfile,'%s%s');
                        if ~isequal(switches.nus4dfpvcnum,prepstem)
                            error('Nusiance 4dfp vcnums do not match the vc ordering of the datalist');
                        end
                        
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.motionestimates=input('Do you want to also regress motion estimates and derivatives? (1=y; 0=no): ');
                            switch switches.motionestimates
                                case {0,1}
                                    needcorrectinput1=0;
                            end
                        end
                        
                        % user txt file
                    case 3
                        needcorrectinput2=0;
                        switches.nustxtlistfile=input('Enter the hard path to the list of nuisance regressor files ','s');
                        [switches.nustxtvcnum switches.nustxtlist] = textread(switches.nustxtlistfile,'%s%s');
                        if ~isequal(switches.nustxtvcnum,prepstem)
                            error('Nusiance txt vcnums do not match the vc ordering of the datalist');
                        end
                        
                        needcorrectinput1=1;
                        while needcorrectinput1
                            switches.motionestimates=input('Do you want to also regress motion estimates and derivatives? (1=y; 0=no): ');
                            switch switches.motionestimates
                                case {0,1}
                                    needcorrectinput1=0;
                            end
                        end
                        
                    case 9
                        needcorrectinput=0;
                        needcorrectinput2=0;
                        switches.motionestimates=1;
                        switches.WMV=1.
                        switches.GS=1;
                        switches.dicesize=input('What size dice to use? (# voxels): ');
                        switches.mindicevox=input('What is minimum # voxels needed in cubes?: ');
                        switches.tcchop=input('What size timeseries to use? (TRs): ');
                        switches.varexpl=input('How much variance should SVD explain?: ');
                        switches.WMero=input('How many WM erosions (4 recommended)?: ');
                        switches.CSFero=input('How many CSF erosions (1 recommended)?: ');
                        switches.sdval=input('What s.d. threshold to form nuisance mask? ');
                        
                end
            end
        case 0
            needcorrectinput=0;
    end
end

% interpolate
needcorrectinput=1;
while needcorrectinput
    switches.dointerpolate=input('Do you want to interpolate over motion epochs? (1=y; 0=no): ');
    switch switches.dointerpolate
        case {0,1}
            needcorrectinput=0;
    end
end

% temporal filter lowpass
needcorrectinput=1;
while needcorrectinput
    switches.dobandpass=input('Do you want to temporally filter the data (1=y; 0=no): ');
    switch switches.dobandpass
        case 1
            needcorrectinput=0;
        case 0
            needcorrectinput=0;
    end
end

if switches.dobandpass
    needcorrectinput=1;
    while needcorrectinput
        switches.temporalfiltertype=input('What type of filter: (1) lowpass, (2) highpass (3) bandpass: ');
        switch switches.temporalfiltertype
            case {1,2,3}
                needcorrectinput=0;
        end
    end
    
    if switches.temporalfiltertype
        switch switches.temporalfiltertype
            case 1
                needcorrectinput1=1;
                while needcorrectinput1
                    switches.lopasscutoff=input('What low-pass cutoff is desired (in Hz; .08 is standard): ');
                    if isnumeric(switches.lopasscutoff)
                        needcorrectinput1=0;
                    end
                end
            case 2
                needcorrectinput1=1;
                while needcorrectinput1
                    switches.hipasscutoff=input('What high-pass cutoff is desired (in Hz; .009 is standard): ');
                    if isnumeric(switches.hipasscutoff)
                        needcorrectinput1=0;
                    end
                end
            case 3
                needcorrectinput1=1;
                while needcorrectinput1
                    switches.lopasscutoff=input('What low-pass cutoff is desired (in Hz; .08 is standard): ');
                    if isnumeric(switches.lopasscutoff)
                        needcorrectinput1=0;
                    end
                end
                needcorrectinput1=1;
                while needcorrectinput1
                    switches.hipasscutoff=input('What high-pass cutoff is desired (in Hz; .009 is standard): ');
                    if isnumeric(switches.hipasscutoff)
                        needcorrectinput1=0;
                    end
                end
                if switches.lopasscutoff <= switches.hipasscutoff
                    fprintf('Low-pass cutoff must be higher than the high-pass cutoff\n');
                end
        end
        needcorrectinput1=1;
        while needcorrectinput1
            switches.order=input('What filter order is desired (1 is standard): ');
            if isnumeric(switches.order)
                needcorrectinput1=0;
            end
        end
        
    end
    
end

% blurring
needcorrectinput=1;
while needcorrectinput
    switches.doblur=input('Do you want to spatially blur the data (1=y; 0=no): ');
    switch switches.doblur
        case 1
            needcorrectinput=0;
            needcorrectinput1=1;
            while needcorrectinput1
                switches.blurkernel=input('What blurring kernel do you want (in mm; 6 is standard): ');
                if isnumeric(switches.blurkernel)
                    needcorrectinput1=0;
                    
                    fprintf('Blur kernel is %d mm\n',switches.blurkernel);
                    
                    
                end
            end
        case 0
            needcorrectinput=0;
    end
end

if switches.doblur
    blursize=2*log(2)/pi*10/switches.blurkernel;
    fprintf('gauss_4dfp set to blur at %d; default is .735452\n',blursize);
end

fprintf('\n\n\n\n*** HERE ARE YOUR SETTINGS ***:\n')
fprintf('switches.doregression (1=yes;0=no): %d\n',switches.doregression);
fprintf('switches.regressiontype (1=freesurfer): %d\n',switches.regressiontype);
fprintf('switches.motionestimates (0=no; 1=R,R`; 2=FRISTON; 20=R,R`,12rand): %d\n',switches.motionestimates);
fprintf('switches.WM (1=regress;0=no): %d\n',switches.WM);
fprintf('switches.V (1=regress;0=no): %d\n',switches.V);
fprintf('switches.GS (1=regress;0=no): %d\n',switches.GS);
fprintf('switches.dointerpolate (1=yes;0=no): %d\n',switches.dointerpolate);
fprintf('switches.dobandpass (1=yes;0=no): %d\n',switches.dobandpass);
fprintf('switches.temporalfiltertype (1=lowpass;2=hipass;3=bandpass): %d\n',switches.temporalfiltertype);
fprintf('switches.lopasscutoff (in Hz; 0.08 is typical): %g\n',switches.lopasscutoff);
fprintf('switches.hipasscutoff (in Hz; 0.009 is typical): %g\n',switches.hipasscutoff);
fprintf('switches.order (1 is typical): %g\n',switches.order);
fprintf('switches.doblur (1=yes;0=no): %d\n',switches.doblur);
fprintf('switches.blurkernel (in mm; 6 is typical): %d\n',switches.blurkernel);

cont=input('\nDo you wish to continue with these settings? (1=yes) ' );
if cont~=1
    error('Quitting');
end

tic
%% CHECK FOR PRESENCE OF STRUCTURAL AND BOLD FILES

% CHECK FOR BOLD DATA
for i=1:numsubs
    fprintf('CHECKING BOLD RUNS\t%d\t%s\n',i,QC(i).vcnum);
    % determine which BOLD runs are rest runs
    clear tempruns; [trash tempruns] = system([ 'awk -F "(" ''{print$2}'' ' df.prmfile{i,1} ' | awk -F ")" ''{print $1}''' ]);
    QC(i).restruns = strread(tempruns,'%s','delimiter',' ');
   
    
    % cycle through each BOLD run
    for j=1:size(QC(i).restruns,1)
        
        % set original data
        boldmat{i,j} = cell2mat(strcat(df.filepath{i,1},'/',QC(i).vcnum,'/bold',QC(i).restruns(j,:),'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d.mat'));
        boldimg{i,j} = cell2mat(strcat(df.filepath{i,1},'/',QC(i).vcnum,'/bold',QC(i).restruns(j,:),'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d_uwrp_atl.4dfp.img'));
        boldifh{i,j} = cell2mat(strcat(df.filepath{i,1},'/',QC(i).vcnum,'/bold',QC(i).restruns(j,:),'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d_uwrp_atl.4dfp.ifh'));
        
        % check that original data is present
        if ~exist(boldmat{i,j})
            error('\t%s is not found\n',boldmat{i,j});
        end
        if ~exist(boldimg{i,j})
            error('\t%s is not found\n',boldimg{i,j});
        end
        if ~exist(boldifh{i,j})
            error('\t%s is not found\n',boldifh{i,j});
        end
    end
    
    anataveimg{i,1} = [ df.filepath{i,1} '/' QC(i).vcnum '/unwarp_mean/' QC(i).vcnum '_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.img'];
    anataveifh{i,1} = [ df.filepath{i,1} '/' QC(i).vcnum '/unwarp_mean/' QC(i).vcnum '_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.ifh'];
    
    mprimg{i,1}=[ df.filepath{i,1} '/' QC(i).vcnum '/atlas/' QC(i).vcnum '_mpr_n1_' voxdim '_t88.4dfp.img'];
    mprifh{i,1}=[ df.filepath{i,1} '/' QC(i).vcnum '/atlas/' QC(i).vcnum '_mpr_n1_' voxdim '_t88.4dfp.ifh'];
  %  mprimg{i,1}(end)=[];
  %  mprifh{i,1}(end)=[];
    
    if ~exist(anataveimg{i,1})
        error('\t%s is not found\n',anataveimg{i,1});
    end
    if ~exist(mprimg{i,1})
        error('\t%s is not found\n',mprimg{i,1});
    end
    
end


%% LINKING TO BOLD DATA

% prepare output directory
if ~exist(tdir)
    mkdir(tdir)
end
fprintf('PREPARING OUTPUT DIRECTORIES\n');
fprintf('LINKING BOLD DATA\n');
pause(2);
for i=1:numsubs
    
    % prepare target subject directory
    QC(i).subdir=[tdir '/' QC(i).vcnum ];
    if exist(QC(i).subdir)
        rmdir(QC(i).subdir,'s');
    end
    mkdir(QC(i).subdir);
    
    % make links to atlas and seed data
    QC(i).subatlasdir=[QC(i).subdir '/atlas/'];
    if ~exist(QC(i).subatlasdir)
        mkdir(QC(i).subatlasdir);
    end
    
    % set symbolic link to structural data
    tanataveimg{i,1} = [ QC(i).subatlasdir '/' QC(i).vcnum '_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.img'];
    tanataveifh{i,1} = [ QC(i).subatlasdir '/' QC(i).vcnum '_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.ifh'];
    system([ 'ln -s ' anataveimg{i,1} ' ' tanataveimg{i,1}]);
    system([ 'ln -s ' anataveifh{i,1} ' ' tanataveifh{i,1}]);
    
    tmprimg{i,1} = [ QC(i).subatlasdir '/' QC(i).vcnum '_mpr_n1_' voxdim '_t88.4dfp.img'];
    tmprifh{i,1} = [ QC(i).subatlasdir '/' QC(i).vcnum '_mpr_n1_' voxdim '_t88.4dfp.ifh'];
    system([ 'ln -s ' mprimg{i,1} ' ' tmprimg{i,1}]);
    system([ 'ln -s ' mprifh{i,1} ' ' tmprifh{i,1}]);
    
    system([ 'cp ' df.prmfile{i,1} ' ' QC(i).subdir '/fc.params' ]);
    
    QC(i).MPRAGE=read_4dfpimg_HCP(mprimg{i,1});
    QC(i).ANATAV=read_4dfpimg_HCP(anataveimg{i,1});
    
    % cycle through each BOLD run
    for j=1:size(QC(i).restruns,1)
        
        % prepare and enter targetsubbolddir
        QC(i).subbolddir{j}=cell2mat(strcat(QC(i).subdir,'/bold',QC(i).restruns(j,:)));
        if ~exist(QC(i).subbolddir{j})
            mkdir(QC(i).subbolddir{j});
        end
        cd(QC(i).subbolddir{j});
        
        % set symbolic link to original data
        tboldmat{i,j} = cell2mat(strcat(QC(i).subbolddir{j},'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d.mat'));
        tboldimg{i,j} = cell2mat(strcat(QC(i).subbolddir{j},'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d_uwrp_atl.4dfp.img'));
        tboldifh{i,j} = cell2mat(strcat(QC(i).subbolddir{j},'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d_uwrp_atl.4dfp.ifh'));
        system([ 'ln -s ' boldmat{i,j} ' ' tboldmat{i,j}]);
        system([ 'ln -s ' boldimg{i,j} ' ' tboldimg{i,j}]);
        system([ 'ln -s ' boldifh{i,j} ' ' tboldifh{i,j}]);
    end
end


%% SET UP DUMMY HOLDING VARIABLE FOR LAST PROCESSED DATA

%%%% To make this modular, we are using LASTIMG to index the previous
%%%% processed image. This allows us to pick and move processing steps

stage=1;
LASTCONC{stage}= '333';
allends='333';

% initialize LASTIMG to the 333 image
for i=1:numsubs
    for j=1:size(QC(i).restruns,1)
        LASTIMG{i,j,stage}=cell2mat(strcat(QC(i).subbolddir{j},'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d_uwrp_atl'));
    end
    cd(QC(i).subdir);
    
    % write a concfile
    fid=fopen([LASTCONC{stage} '.conc'],'w');
    fprintf(fid,'number_of_files: %d\n',size(QC(i).restruns,1));
    for j=1:size(QC(i).restruns,1)
        fprintf(fid,'\tfile:%s\n',[LASTIMG{i,j,stage} '.4dfp.img' ]);
    end
    fclose(fid);
end


%% COMPUTE DEFINED VOXELS FOR THE BOLD RUNS

for i=1:numsubs
    cd(QC(i).subdir);
    fprintf('COMPUTE DEFINED VOXELS\t%d\t%s\n',i,QC(i).vcnum);
    commands=['compute_defined_4dfp ' QC(i).subdir '/' LASTCONC{stage} '.conc >/dev/null' ];
    system(commands);
    pause(5);
    QC(i).DFNDVOXELS=read_4dfpimg_HCP([LASTCONC{stage} '_dfnd.4dfp.img']);
end


%% CHECK FOR NUISANCE SEEDS

%if switches.doregression
    
    needtostop=0;
    switch switches.regressiontype
%         case 0 % classic WashU seeds in WM and V
%             
%             for i=1:numsubs
%                 QC(i).GLMmaskfile=['//data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/freesurfer_washu/FREESURFER_fs_LR/sub013/7112b_fs_LR/brainmask_fs_222.4dfp.img'];% ['/home/usr/fidl/lib/glm_atlas_mask_' voxdim '.4dfp.img'];
%                 QC(i).WMmaskfile=['/home/usr/fidl/lib/small_WM.4dfp.img'];
%                 QC(i).CSFmaskfile=['/home/usr/fidl/lib/711-2B_lat_vent_' voxdim '.4dfp.img'];
%             end
%             
%             for i=1:numsubs
%                 fprintf('CHECKING NUISANCE SEEDS\t%d\t%s\n',i,QC(i).vcnum);
%                 needtostop=0;
%                 
%                 tmpmask=read_4dfpimg_HCP(QC(i).GLMmaskfile);
%                 tmpmask=tmpmask & QC(i).DFNDVOXELS;
%                 if nnz(tmpmask)==0
%                     fprintf('%s is empty!\n',QC(i).GLMmaskfile);
%                     needtostop=1;
%                 else
%                     QC(i).GLMMASK=~~tmpmask;
%                 end
%                 
%                 tmpmask=read_4dfpimg_HCP(QC(i).WMmaskfile);
%                 tmpmask=tmpmask & QC(i).DFNDVOXELS;
%                 if nnz(tmpmask)==0
%                     fprintf('%s is empty!\n',QC(i).WMmaskfile);
%                     needtostop=1;
%                 else
%                     QC(i).WMMASK=~~tmpmask;
%                 end
%                 
%                 tmpmask=read_4dfpimg_HCP(QC(i).CSFmaskfile);
%                 tmpmask=tmpmask & QC(i).DFNDVOXELS;
%                 if nnz(tmpmask)==0
%                     fprintf('%s is empty!\n',QC(i).CSFmaskfile);
%                     needtostop=1;
%                 else
%                     QC(i).CSFMASK=~~tmpmask;
%                 end
%                 
%                 if needtostop
%                     error('Fix empty masks.\n');
%                 end
%             end
%             
%             % create a holding variable for seeds and labels
%             for i=1:numsubs
%                 seedcount=0;
%                 QC(i).seedhold=[];
%                 if switches.GS
%                     seedcount=seedcount+1;
%                     QC(i).seedhold=[QC(i).seedhold QC(i).GLMMASK];
%                     QC(i).seedholdlabel{seedcount}='WB';
%                 end
%                 if switches.WM
%                     seedcount=seedcount+1;
%                     QC(i).seedhold=[QC(i).seedhold QC(i).WMMASK];
%                     QC(i).seedholdlabel{seedcount}='WM';
%                 end
%                 if switches.V
%                     seedcount=seedcount+1;
%                     QC(i).seedhold=[QC(i).seedhold QC(i).CSFMASK];
%                     QC(i).seedholdlabel{seedcount}='CSF';
%                 end
%             end
            
        case {0,1,9} % freesurfer masks of WM and V
            disp('hello')
            % set basic names
            for i=1:numsubs
                QC(i).GLMmaskfile=['/data/cn4/laumannt/Standard/glm_atlas_mask_' voxdim '.4dfp.img'];
                QC(i).WMmaskfile=[ QC(i).fsdir '/nusmask/aparc+aseg_cerebralwm_ero' num2str(switches.WMero) '_mask_' voxdim '.4dfp.img'];
                QC(i).CSFmaskfile=[ QC(i).fsdir '/nusmask/aparc+aseg_CSF_ero' num2str(switches.CSFero) '_mask_' voxdim '.4dfp.img'];
                QC(i).WBmaskfile=[ QC(i).fsdir '/nusmask/aparc+aseg_brainmask_mask_' voxdim '.4dfp.img'];
                QC(i).GREYmaskfile=[ QC(i).fsdir '/nusmask/aseg_GM_mask_' voxdim '.4dfp.img'];
            end
    
            % check for existence of mask files
            for i=1:numsubs
                fprintf('CHECKING NUISANCE SEEDS\t%d\t%s\n',i,QC(i).vcnum);
                needtostop=0;
                
                if ~exist([QC(i).GLMmaskfile])
                    fprintf('No GLMmask found: %s\n',QC(i).GLMmaskfile);
                    needtostop=1;
                end
                
                if ~exist([ QC(i).fsdir '/' ])
                    fprintf('No FREESURFER folder! Subject %s needs to have mpr2caret run\n',QC(i).vcnum);
                    needtostop=1;
                end
                if ~exist(QC(i).WBmaskfile)
                    fprintf('%s missing!\n',QC(i).WBmaskfile);
                    needtostop=1;
                end
                if ~exist(QC(i).GREYmaskfile)
                    fprintf('%s missing!\n',QC(i).GREYmaskfile);
                    needtostop=1;
                end
                if ~exist(QC(i).WMmaskfile)
                    fprintf('%s missing!\n',QC(i).WMmaskfile);
                    needtostop=1;
                end
                if ~exist(QC(i).CSFmaskfile)
                    fprintf('%s missing!\n',QC(i).CSFmaskfile);
                    needtostop=1;
                end
                if needtostop
                    error('Fix the FREESURFER masks.\n');
                end
            end
            
            % ensure masks contain something, relax erosions if not
            for i=1:numsubs
                needtostop=0;
                
                tmpmask=read_4dfpimg(QC(i).GLMmaskfile);
                tmpmask=tmpmask & QC(i).DFNDVOXELS;
                if nnz(tmpmask)==0
                    fprintf('%s is empty!\n',QC(i).GLMmaskfile);
                    needtostop=1;
                else
                    QC(i).GLMMASK=~~tmpmask;
                end
                
                tmpmask=read_4dfpimg_HCP(QC(i).WBmaskfile);
                tmpmask=tmpmask & QC(i).DFNDVOXELS;
                if nnz(tmpmask)==0
                    fprintf('%s is empty!\n',QC(i).WBmaskfile);
                    needtostop=1;
                else
                    QC(i).WBMASK=~~tmpmask;
                end
                
                tmpmask=read_4dfpimg_HCP(QC(i).GREYmaskfile);
                tmpmask=tmpmask & QC(i).DFNDVOXELS;
                if nnz(tmpmask)==0
                    fprintf('%s is empty!\n',QC(i).GREYmaskfile);
                    needtostop=1;
                else
                    QC(i).GMMASK=~~tmpmask;
                end
                
                tmpmask=read_4dfpimg_HCP(QC(i).WMmaskfile);
                tmpmask=tmpmask & QC(i).DFNDVOXELS;
                tmpWMero=switches.WMero;
                while nnz(tmpmask)==0
                    tmpWMero=tmpWMero-1;
                    if tmpWMero==-1
                        fprintf('Subject %s has no white matter voxels!\n',QC(i).vcnum)
                        needtostop=1;
                        break;
                    end
                    fprintf('Reducing WM erosion to %d; zero voxels in %s\n',tmpWMero,QC(i).WMmaskfile);
                    QC(i).WMmaskfile = [ QC(i).fsdir '/nusmask/aparc+aseg_cerebralwm_ero' num2str(tmpWMero) '_mask_' voxdim '.4dfp.img'];
                    tmpmask=read_4dfpimg_HCP(QC(i).WMmaskfile);
                    tmpmask=tmpmask & QC(i).DFNDVOXELS;
                end
                QC(i).WMMASK=~~tmpmask;
                QC(i).WMero=tmpWMero;
                
                tmpmask=read_4dfpimg_HCP(QC(i).CSFmaskfile);
                tmpmask=tmpmask & QC(i).DFNDVOXELS;
                tmpCSFero=switches.CSFero;
                while nnz(tmpmask)==0
                    tmpCSFero=tmpCSFero-1;
                    if tmpCSFero==-1
                        fprintf('Subject %s has no CSF voxels!\n',QC(i).vcnum)
                        needtostop=1;
                        break;
                    end
                    fprintf('Reducing CSF erosion to %d; zero voxels in %s\n',tmpCSFero,QC(i).CSFmaskfile);
                    QC(i).CSFmaskfile = [ QC(i).fsdir '/nusmask/aparc+aseg_cerebralwm_ero' num2str(tmpCSFero) '_mask_' voxdim '.4dfp.img'];
                    tmpmask=read_4dfpimg_HCP(QC(i).CSFmaskfile);
                end
                QC(i).CSFMASK=~~tmpmask;
                QC(i).CSFero=tmpCSFero;
                
                if needtostop
                    error('Fix empty masks.\n');
                end
            end
            
            % create a holding variable for seeds and labels
            for i=1:numsubs
                seedcount=0;
                QC(i).seedhold=[];
                if switches.GS
                    seedcount=seedcount+1;
                    QC(i).seedhold=[QC(i).seedhold QC(i).WBMASK];
                    QC(i).seedholdlabel{seedcount}='WB';
                end
                if switches.WM
                    seedcount=seedcount+1;
                    QC(i).seedhold=[QC(i).seedhold QC(i).WMMASK];
                    QC(i).seedholdlabel{seedcount}='WM';
                end
                if switches.V
                    seedcount=seedcount+1;
                    QC(i).seedhold=[QC(i).seedhold QC(i).CSFMASK];
                    QC(i).seedholdlabel{seedcount}='CSF';
                end
            end
            
        case 2 % external 4dfp of regressor ROIs
            
            
        case 3 % external txt file
            
        otherwise
    end
%end


%% LINK NUISANCE SEEDS

for i=1:numsubs
    
    QC(i).subseeddir=[QC(i).subdir '/nusseeds/'];
    if ~exist(QC(i).subseeddir)
        mkdir(QC(i).subseeddir);
    end
    
    cd(QC(i).subdir);
    system(['mv ' LASTCONC{stage} '_dfnd.* ' QC(i).subseeddir ]);
    
 %   if switches.doregression
        switch switches.regressiontype
%             case 0
%                 [pth tmpseedimg]=fileparts(QC(i).GLMmaskfile);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
%                 [pth tmpseedimg]=fileparts(QC(i).WMmaskfile);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
%                 [pth tmpseedimg]=fileparts(QC(i).CSFmaskfile);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
%                 system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
%                 write_4dfpimg(QC(i).seedhold,[ QC(i).subseeddir '/usedseeds.4dfp.img'],'bigendian');
%                 write_4dfpifh_222([ QC(i).subseeddir '/usedseeds.4dfp.img'],size(QC(i).seedhold,2),'bigendian',128,128,75,2,2,2);
%                 
            case {0,1,9}
                [pth tmpseedimg]=fileparts(QC(i).GLMmaskfile);
                system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
                system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
                [pth tmpseedimg]=fileparts(QC(i).WMmaskfile);
                system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
                system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
                [pth tmpseedimg]=fileparts(QC(i).CSFmaskfile);
                system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
                system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
                [pth tmpseedimg]=fileparts(QC(i).GREYmaskfile);
                system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
                system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
                [pth tmpseedimg]=fileparts(QC(i).WBmaskfile);
                system([ 'ln -s ' pth '/' tmpseedimg '.img ' QC(i).subseeddir '/' tmpseedimg '.img' ]);
                system([ 'ln -s ' pth '/' tmpseedimg '.ifh ' QC(i).subseeddir '/' tmpseedimg '.ifh' ]);
                write_4dfpimg(QC(i).seedhold,[ QC(i).subseeddir '/usedseeds.4dfp.img'],'bigendian');
                write_4dfpifh([ QC(i).subseeddir '/usedseeds.4dfp.img'],size(QC(i).seedhold,2),'bigendian',128,128,75,2,2,2);
                
            case 2
                
        end
    end
%end


%% CALCULATE SUBJECT MOVEMENT
for i=1:numsubs
    for j=1:size(QC(i).restruns,1)
        
        cd(QC(i).subbolddir{j});
        fprintf('CALCULATING MOTION\t%d\t%s\trun%s\n',i,QC(i).vcnum,cell2mat(QC(i).restruns(j,:)));
        
        % obtain realignment parameters
        system([ 'mat2dat -R -D -L ' tboldmat{i,j} ' >/dev/null' ]);
        pause(5);
        
        rdatfile{i,j}= cell2mat(strcat(df.filepath{i,1},'/',QC(i).vcnum,'/movement/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d.rdat'));
        ddatfile{i,j}= cell2mat(strcat(df.filepath{i,1},'/',QC(i).vcnum,'/movement/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_faln_dbnd_xr3d.ddat'));
        
        % convert realignment parameters to FD, write them
        [rmstotal{i,j} rmstrans rmsrot rmscol mvm{i,j}] = rdat_calculations(rdatfile{i,j},QC(i).TRskip,50);
        [drmstotal{i,j} drmstrans drmsrot drmscol ddt_mvm{i,j}] = rdat_calculations(ddatfile{i,j},QC(i).TRskip,50);
        FD{i,j}=(sum(abs(ddt_mvm{i,j}),2));
        dlmwrite('FDvals.txt',FD{i,j},'\t');
        
        % write the movement parameters
        dlmwrite('movement.txt',mvm{i,j},'\t');
        dlmwrite('ddt_movement.txt',ddt_mvm{i,j},'\t');
        
        % detrend the movement
        system(['cat movement.txt | gawk -f ' releasedir '/trendout.awk >! movement_detrended.txt']);
        system(['cat ddt_movement.txt | gawk -f ' releasedir '/trendout.awk >! ddt_movement_detrended.txt']);
        
        % load in the detrended regressors
        mvm_detrend{i,j}=load('movement_detrended.txt');
        ddt_mvm_detrend{i,j}=load('ddt_movement_detrended.txt');
    end
    
    % WRITE TOTAL DATA FOR EACH SUBJECT
    
    cd(QC(i).subdir);
    
    % write the total movement data
    QC(i).MVM=[];
    QC(i).ddtMVM=[];
    QC(i).DTMVM=[];
    QC(i).ddtDTMVM=[];
    QC(i).FD=[];
    for j=1:size(QC(i).restruns,1)
        QC(i).MVM=[QC(i).MVM; mvm{i,j}];
        QC(i).ddtMVM=[QC(i).ddtMVM; ddt_mvm{i,j}];
        QC(i).DTMVM=[QC(i).DTMVM; mvm_detrend{i,j}];
        QC(i).ddtDTMVM=[QC(i).ddtDTMVM; ddt_mvm_detrend{i,j}];
        QC(i).FD=[QC(i).FD; FD{i,j}];
    end
    dlmwrite('total_movement.txt',QC(i).MVM,'\t');
    dlmwrite('total_ddt_movement.txt',QC(i).ddtMVM,'\t');
    dlmwrite('total_movement_detrended.txt',QC(i).DTMVM,'\t');
    dlmwrite('total_ddt_movement_detrended.txt',QC(i).ddtDTMVM,'\t');
    dlmwrite('total_FD.txt',QC(i).FD,'\t');
    
    QC(i).MVMrms=array_calc_rms(QC(i).MVM);
    QC(i).DTMVMrms=array_calc_rms(QC(i).DTMVM);
    QC(i).FDbar=mean(QC(i).FD);
    QC(i).switches=switches;
    
end


%% DEFINE RUN BORDERS
for i=1:numsubs
    trpos=0;
    tr(i).tot=numel(QC(i).FD);
    for j=1:size(QC(i).restruns,1)
        tr(i).runtrs(j)=numel(FD{i,j});
        tr(i).start(j,1:2)=[trpos+1 trpos+tr(i).runtrs(j)];
        trpos=trpos+tr(i).runtrs(j);
        QC(i).runborders(j,:) = [j tr(i).start(j,1:2)];
    end
   % QC(i).runborders=[QC(i).restruns' tr(i).start];
    cd(QC(i).subdir);
    dlmwrite('runborders.txt',QC(i).runborders,'\t');
end

[sortTR sortsubs]=sort(cell2mat({tr.tot}));





%% ASSEMBLE TEMPORAL MASKS
switch tmasktype
    case 'ones'
        for i=1:numsubs
            for j=1:size(QC(i).restruns,1)
                QC(i).runtmask{j}=ones(size(FD{i,j},1),1);
                QC(i).runtmask{j}(1:QC(i).TRskip)=0;
            end
            QC(i).tmask=[];
            for j=1:size(QC(i).restruns,1)
                QC(i).tmask=[QC(i).tmask; QC(i).runtmask{j}];
            end
        end
    otherwise
        for i=1:numel(tm.tmaskfiles)
            fprintf('GETTING TMASK FILE\t%d\t%s\t%s\n',i,QC(i).vcnum,tm.tmaskfiles{i,1});
            QC(i).tmask=load(tm.tmaskfiles{i,1});
            if ~isequal(numel(QC(i).tmask),numel(QC(i).FD))
                error('tmask length does not match the subject data');
            end
            
            for j=1:size(QC(i).restruns,1)
                QC(i).runtmask{j}=QC(i).tmask(QC(i).runborders(j,2):QC(i).runborders(j,3));
            end
        end
end

for i=1:numsubs
    cd(QC(i).subdir);
    dlmwrite('total_tmask.txt',QC(i).tmask,'\t');
    tmask2format(QC(i).tmask,'total_tmask_avi.txt');
    QC(i).tmask=~~QC(i).tmask;
    for j=1:size(QC(i).restruns,1)
        cd(QC(i).subbolddir{j});
        dlmwrite('tmask.txt',QC(i).runtmask{j},'\t');
        tmask2format(QC(i).runtmask{j},'tmask_avi.txt');
        QC(i).runtmask{j}=~~QC(i).runtmask{j};
    end
end


%% FUNCTIONAL CONNECTIVITY PROCESSING
bigstuff=1; % this saves voxelwise timecourses over processing.
skipvox=15; % downsample grey matter voxels for visuals.
set(0, 'DefaultFigureVisible', 'off'); 
for f=1:numsubs
    
    % orders subjects in decreasing order to minimize memory fragmentation
    i=sortsubs(f);
    
    fprintf('FCPROCESSING SUBJECT %d %s\n',f,QC(i).vcnum);
    pause(2);
    
    %Select voxels in glmmask
    QC(i).CSFMASK_glmmask = QC(i).CSFMASK(logical(QC(i).GLMMASK));
    QC(i).WMMASK_glmmask = QC(i).WMMASK(logical(QC(i).GLMMASK));
    QC(i).GMMASK_glmmask = QC(i).GMMASK(logical(QC(i).GLMMASK));
    QC(i).WBMASK_glmmask = QC(i).WBMASK(logical(QC(i).GLMMASK));
    QC(i).GLMMASK_glmmask = QC(i).GLMMASK(logical(QC(i).GLMMASK));
    roimasks =roimasks_orig(logical(QC(i).GLMMASK),:);
   
    %%%
    % THE PROCESSING BEGINS
    %%%
    
    stage=1;
    LASTCONC{stage}= '333';
    allends='333';
    for j=1:size(QC(i).restruns,1)
        LASTIMG{i,j,stage}=cell2mat(strcat(QC(i).subbolddir{j},'/',QC(i).vcnum,'_b',QC(i).restruns(j,:),'_xr3d_uwrp_atl'));
    end
    
    cd(QC(i).subdir);
    
    % obtain the raw images
    ending='333';
    [bolds]=conc2bolds([LASTCONC{stage} '.conc']);
    [tempimg]=bolds2img(bolds,tr(i).tot,tr(i).start,QC(i).GLMMASK);
    
    [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
    QC(i).process{stage}=ending;
    QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
    dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
    if bigstuff
        tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
        QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
        QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
        QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
    end
    makepictures(QC(i),stage,switches,[700:200:1300],[0:50:100],50);
    saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% 0-mean, detrend %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    stage=stage+1;
    ending='zmdt';
    allends=[allends '_' ending];
    for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
    end
    LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
    
    % for each BOLD run
    for j=1:size(QC(i).restruns,1)
        fprintf('\tDEMEAN DETREND\trun%s\n',cell2mat(QC(i).restruns(j,:)));
        tic;
        temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
        temprun=demean_detrend(temprun,QC(i).runtmask{j});
        tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3))=temprun;
        toc;
    end
    
    [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
    QC(i).process{stage}=ending;
    QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
    dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
    if bigstuff
        tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
        QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
        QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
        QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
    end
    makepictures(QC(i),stage,switches,[-20:20:20],[0:50:100],50);
    saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% MULTIPLE REGRESSION %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % load the image in question, including all BOLD runs
    fprintf('\tNUISANCE REGRESSION\n');
    stage=stage+1;
    ending='resid';
    allends=[allends '_' ending];
    for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
        end
    LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
    
    
    if switches.doregression
        
        % get the movement-based regressors
        switch switches.motionestimates
            case 0 %
                QC(i).mvmregs=[];
                QC(i).mvmlabels={''};
            case 1 % R,R`                   LAB CLASSIC
                QC(i).mvmregs=[QC(i).DTMVM QC(i).ddtDTMVM];
                QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','X`','Y`','Z`','pitch`','yaw`','roll`'};
            case 2 % R,R^2,R-1,R-1^2       FRISTON
                frist1=circshift(QC(i).DTMVM,[1 0]);
                frist1(1,:)=0;
                QC(i).mvmregs=[QC(i).DTMVM (QC(i).DTMVM.^2) frist1 frist1.^2 ];
                QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','sqrX','sqrY','sqrZ','sqrpitch','sqryaw','sqrroll`','Xt-1','Yt-1','Zt-1','pitcht-1','yawt-1','rollt-1`','sqrXt-1','sqrYt-1','sqrZt-1','sqrpitcht-1','sqryawt-1','sqrrollt-1`'};
            case 20 % R,R`,12               CONTROL for 24 parameters
                QC(i).mvmregs=[QC(i).DTMVM QC(i).ddtDTMVM rand(size(QC(i).DTMVM,1),12) ];
                QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','X`','Y`','Z`','pitch`','yaw`','roll`','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12'};
            case 3
                frist=[QC(i).DTMVM];
                frist1=circshift(QC(i).DTMVM,[1 0]);
                frist1(1,:)=0;
                frist2=circshift(QC(i).DTMVM,[2 0]);
                frist2(1:2,:)=0;
                QC(i).mvmregs=[QC(i).DTMVM (QC(i).DTMVM.^2) frist1 frist1.^2 frist2 frist2.^2 ];
                QC(i).mvmlabels={'X','Y','Z','pitch','yaw','roll`','sqrX','sqrY','sqrZ','sqrpitch','sqryaw','sqrroll`','Xt-1','Yt-1','Zt-1','pitcht-1','yawt-1','rollt-1`','sqrXt-1','sqrYt-1','sqrZt-1','sqrpitcht-1','sqryawt-1','sqrrollt-1`','Xt-2','Yt-2','Zt-2','pitcht-2','yawt-2','rollt-2`','sqrXt-2','sqrYt-2','sqrZt-2','sqrpitcht-2','sqryawt-2','sqrrollt-2`'};
        end
        
        
        % get the signal regressors
        QC(i).sigregs=[];
        QC(i).siglabels=[];
        switch switches.regressiontype
            case {0,1}
                if switches.GS
                    if nnz(QC(i).GLMMASK_glmmask)
                        sig=mean(tempimg(QC(i).GLMMASK_glmmask,:))';
                        QC(i).sigregs=[QC(i).sigregs sig];
                        QC(i).siglabels=[QC(i).siglabels {'WB'}];
                    else
                        fprintf('Whole brain mask empty: %s\n',QC(i).GLMmaskfile);
                    end
                end
                if switches.WM
                    if nnz(QC(i).WMMASK_glmmask)
                        sig=mean(tempimg(QC(i).WMMASK_glmmask,:))';
                        QC(i).sigregs=[QC(i).sigregs sig];
                        QC(i).siglabels=[QC(i).siglabels {'WM'}];
                    else
                        fprintf('Whole brain mask empty: %s\n',QC(i).WMMmaskfile);
                    end
                end
                if switches.V
                    if nnz(QC(i).CSFMASK_glmmask)
                        sig=mean(tempimg(QC(i).CSFMASK_glmmask,:))';
                        QC(i).sigregs=[QC(i).sigregs sig];
                        QC(i).siglabels=[QC(i).siglabels {'V'}];
                    else
                        fprintf('Whole brain mask empty: %s\n',QC(i).CSFmaskfile);
                    end
                end
                if ~isempty(QC(i).sigregs)
                    QC(i).sigregs=[QC(i).sigregs [repmat(0,[1 size(QC(i).sigregs,2)]); diff(QC(i).sigregs)]];
                    kk=numel(QC(i).siglabels);
                    for k=1:kk
                        QC(i).siglabels{k+kk}=[ QC(i).siglabels{k} '`'];
                    end
                end
        end
        
        QC(i).nuisanceregressors=[QC(i).mvmregs QC(i).sigregs];
        QC(i).nuisanceregressorlabels=[QC(i).mvmlabels QC(i).siglabels];
        dlmwrite('total_nuisance_regressors.txt',QC(i).nuisanceregressors,'\t');
        
        subplot(8,1,8);
        imagesc(zscore(QC(i).nuisanceregressors)',[-2 2]); ylabel('REGS');
        saveas(gcf,'total_nuisance_regressors.tiff','tiff');
        
        % write the correlations of the nuisance regressors
        clf;
        h=imagesc(triu(corrcoef(QC(i).nuisanceregressors),1),[-.5 1]);
        colorbar;
        saveas(gcf,'total_nuisance_regressors_correlations.tiff','tiff');
        close;
        dlmwrite('total_nuisance_regressors_correlations.txt',corrcoef(QC(i).nuisanceregressors),'\t');
        close;
        
        tic
        [tempimg zb regsz]=regress_nuisance(tempimg,QC(i).nuisanceregressors,QC(i).tmask);
        toc
        
        QC(i).nuisanceregressors_ZSCORE=regsz;
        
        % write the beta images
        zb_out = zeros(size(QC(i).GLMMASK,1),size(zb,2)); %Put back in volume space
        zb_out(logical(QC(i).GLMMASK),:) = zb; %Put back in volume space
        clear zb
        write_4dfpimg(zb_out,'beta.4dfp.img','bigendian');
        write_4dfpifh('beta.4dfp.img',size(zb_out,2),'bigendian');
        
        
        [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
        QC(i).process{stage}=ending;
        QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
        dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
        if bigstuff
            tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
            QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
            QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
            QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
        end
        makepictures(QC(i),stage,switches,[-20:20:20],[0:50:100],50);
        saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % INTERPOLATION
    %%%%%%%%%%%%%%%%%%%%%%%%
    
%     if switches.dointerpolate
%         
%         stage=stage+1;
%         ending='ntrpl';
%         allends=[allends '_' ending];
%         for j=1:size(QC(i).restruns,1)
%             LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
%         end
%         LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
%         
%         % for each BOLD run
%         for j=1:size(QC(i).restruns,1)
%             fprintf('\tINTERPOLATE\trun%d\n',QC(i).restruns(j,:));
%             tic;
%             temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
%             [tchunksize tstartnums tendnums] = maskgaps(QC(i).runtmask{j});
%             for m=1:size(temprun,1) % for each voxel replace with interpolation
%                 [temprun(m,:)]=ts_interpolate(temprun(m,:),tstartnums,tendnums);
%             end
%             tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3))=temprun;
%             toc;
%         end
%         
%         [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
%         QC(i).process{stage}=ending;
%         QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
%         dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
%         if bigstuff
%             tmptcs=single(tempimg(QC(i).GMMASK,:));
%             QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
%             QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK,:));
%             QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK,:));
%         end
%         makepictures(QC(i),stage,switches,[-20:20:20],[0:50:100],50);
%         saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
%         
%     end

    if switches.dointerpolate
        
        stage=stage+1;
        ending='ntrpl';
        allends=[allends '_' ending];
        for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
        end
        LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
        
        
        % for each BOLD run
        for j=1:size(QC(i).restruns,1)
            fprintf('\tINTERPOLATE\trun%s\n',cell2mat(QC(i).restruns(j,:)));
            tic;
            temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
            ofac=8;
            hifac=1;
            TRtimes=([1:size(temprun,2)]')*QC(i).TR;
            
            if numel(TRtimes)<150
                voxbinsize=5000;
            elseif (numel(TRtimes)>=150 && numel(TRtimes)<500)
                voxbinsize=500;
            elseif numel(TRtimes)>=500
                voxbinsize=50;
            end
            fprintf('INTERPOLATION VOXBINSIZE: %d\n',voxbinsize);
            voxbin=1:voxbinsize:size(temprun,1);
            voxbin=[voxbin size(temprun,1)];
            
            temprun=temprun';
            tempanish=zeros(size(temprun,1),size(temprun,2));
            
            % gotta bin by voxels: 5K is ~15GB, 15K is ~40GB at standard
            % run lengths. 5K is ~15% slower but saves 2/3 RAM, so that's
            % the call.
            disp(['size ' num2str(size(voxbin))])
            matlabpool open 4
            parfor v=1:numel(voxbin)-1 % this takes huge RAM if all voxels
                %disp(['binnum ' num2str(v)])
                temp = temprun(~~QC(i).runtmask{j},voxbin(v):voxbin(v+1));
                tempanish_piece{v}=getTransform(TRtimes(~~QC(i).runtmask{j}),temp,TRtimes,QC(i).TR,ofac,hifac);

               % tempanish(:,voxbin(v):voxbin(v+1))=getTransform(TRtimes(~~QC(i).runtmask{j}),temp),TRtimes,QC(i).TR,ofac,hifac);
            end
            matlabpool close
            
            for v = 1:numel(voxbin)-1
                tempanish(:,voxbin(v):voxbin(v+1)) = tempanish_piece{v};
            end
            clear tempanish_piece;
    %        matlabpool open 4
    %        parfor v=1:size(temprun,2) % this takes huge RAM if all voxels
    %            tempanish(:,v)=getTransform(TRtimes(~~QC(i).runtmask{j}),temprun(~~QC(i).runtmask{j},v),TRtimes,QC(i).TR,ofac,hifac);
    %        end
    %        matlabpool close
            
            tempanish=tempanish';
            temprun=temprun';
            
            temprun(:,~QC(i).runtmask{j})=tempanish(:,~QC(i).runtmask{j});
            tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3))=temprun;
            toc;
        end
        
        [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
        QC(i).process{stage}=ending;
        QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
        dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
        if bigstuff
            tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
            QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
            QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
            QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
        end
        makepictures(QC(i),stage,switches,[-20:20:20],[0:50:100],50);
        saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEMPORAL FILTER %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if switches.dobandpass
        
        stage=stage+1;
        ending='bpss';
        allends=[allends '_' ending];
        for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
        end
        LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
        
        filtorder=switches.order;
        switch switches.temporalfiltertype
            case 1
                lopasscutoff=switches.lopasscutoff/(0.5/QC(i).TR); % since TRs vary have to recalc each time
                [butta buttb]=butter(filtorder,lopasscutoff,'low');
            case 2
                hipasscutoff=switches.hipasscutoff/(0.5/QC(i).TR); % since TRs vary have to recalc each time
                [butta buttb]=butter(filtorder,hipasscutoff,'high');
            case 3
                lopasscutoff=switches.lopasscutoff/(0.5/QC(i).TR); % since TRs vary have to recalc each time
                hipasscutoff=switches.hipasscutoff/(0.5/QC(i).TR); % since TRs vary have to recalc each time
                [butta buttb]=butter(filtorder,[hipasscutoff lopasscutoff]);
        end
        
        for j=1:size(QC(i).restruns,1)
            fprintf('\tTEMPORAL FILTER\trun%s\n',cell2mat(QC(i).restruns(j,:)));
            tic;
            temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
            temprun=temprun';
            size_temprun = size(temprun); 
            pad = 1000; 
            temprun = cat(1, zeros(pad, size_temprun(2)), temprun, zeros(pad, size_temprun(2))); 
            [temprun]=filtfilt(butta,buttb,double(temprun)); 
            temprun = temprun(pad+1:end-pad, 1:size_temprun(2)); 
            temprun=temprun';
            tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3))=temprun;
            toc;
        end
        
        [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
        QC(i).process{stage}=ending;
        QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
        dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
        if bigstuff
            tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
            QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
            QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
            QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
        end
        makepictures(QC(i),stage,switches,[-20:20:20],[0:10:20],10);
        saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
        
        % create temporal mask based on filter properties
        filtertrim=15; % TRs at beginning and end of run to ignore due to IIR zero-phase filter
        QC(i).filtertmask=[];
        
        for j=1:size(QC(i).restruns,1)
            QC(i).runfiltertmask{j}=QC(i).runtmask{j};
            QC(i).runfiltertmask{j}(1:filtertrim)=0;
            QC(i).runfiltertmask{j}(end-filtertrim+1:end)=0;
            QC(i).runfiltertmask{j}=QC(i).runfiltertmask{j} & QC(i).runtmask{j};
            QC(i).filtertmask=[QC(i).filtertmask; QC(i).runfiltertmask{j}];
        end
                     
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% 0-mean, detrend %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    stage=stage+1;
    ending='zmdt';
    allends=[allends '_' ending];
    for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
        end
    LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
    
    % for each BOLD run
    for j=1:size(QC(i).restruns,1)
        fprintf('\tDEMEAN DETREND\trun%s\n',cell2mat(QC(i).restruns(j,:)));
        tic;
        temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
        if switches.dobandpass
            temprun=demean_detrend(temprun,QC(i).runfiltertmask{j});
        else
            temprun=demean_detrend(temprun,QC(i).runtmask{j});
        end
        tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3))=temprun;
        toc;
    end
    
    [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
    QC(i).process{stage}=ending;
    QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
    dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
    if bigstuff
        tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
        QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
        QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
        QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
    end
    makepictures(QC(i),stage,switches,[-20:20:20],[0:10:20],10);
    saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% SPATIAL BLUR %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    if switches.doblur
        stage=stage+1;
        ending=[ 'g' num2str(round(blursize*10)) ];
        allends=[allends '_' ending];
        LASTCONC{stage}=[LASTCONC{stage-1} '_' ending];
        for j=1:size(QC(i).restruns,1)
            LASTIMG{i,j,stage}=[ LASTIMG{i,j,stage-1} '_' ending ];
        end
        
        tic
        fprintf('\tSPATIAL BLUR\n');
        
        % write BOLD files from previous stage b/c use 4dfp tool
        for j=1:size(QC(i).restruns,1)
            temprun=tempimg(:,QC(i).runborders(j,2):QC(i).runborders(j,3));
            temprun_out = zeros(size(QC(i).GLMMASK,1),size(temprun,2)); %Put back in volume space
            temprun_out(logical(QC(i).GLMMASK),:) = temprun; %Put back in volume space
            clear temprun
            write_4dfpimg(temprun_out,[ LASTIMG{i,j,stage-1} '.4dfp.img'],'bigendian');
            write_4dfpifh([ LASTIMG{i,j,stage-1} '.4dfp.img'],size(temprun_out,2),'bigendian');
        end
        % write a new concfile
        fid=fopen([LASTCONC{stage-1} '.conc'],'w');
        fprintf(fid,'number_of_files: %s\n',size(QC(i).restruns,1));
        for j=1:size(QC(i).restruns,1)
            fprintf(fid,'\tfile:%s\n',[LASTIMG{i,j,stage-1} '.4dfp.img' ]);
        end
        fclose(fid);
        
        % run the spatial blur
        system(['gauss_4dfp ' LASTCONC{stage-1} '.conc ' num2str(blursize) ' >/dev/null']);
        pause(5);
        [bolds]=conc2bolds([LASTCONC{stage} '.conc']);
        [tempimg]=bolds2img(bolds,tr(i).tot,tr(i).start);
        toc
        
        [QC]=tempimgsignals(QC,i,tempimg,switches,stage);
        QC(i).process{stage}=ending;
        QC(i).tc(:,:,stage)=gettcs(roimasks,tempimg);
        dlmwrite(['total_DV_' LASTCONC{stage} '.txt'],QC(i).DV_GLM(:,stage));
        if bigstuff
            tmptcs=single(tempimg(QC(i).GMMASK_glmmask,:));
            QC(i).GMtcs(:,:,stage)=tmptcs(1:skipvox:end,:);
            QC(i).WMtcs(:,:,stage)=single(tempimg(QC(i).WMMASK_glmmask,:));
            QC(i).CSFtcs(:,:,stage)=single(tempimg(QC(i).CSFMASK_glmmask,:));
        end
        makepictures(QC(i),stage,switches,[-20:20:20],[0:10:20],10);
        saveas(gcf,[ num2str(i) '_' QC(i).vcnum '_' num2str(stage) '_' allends '.tiff'],'tiff');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONCATENATE INTO SINGLE IMAGE %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    fprintf('\tCONCATENATE AND CLEANUP\n');
    cd(QC(i).subdir);
    tempimg_out = zeros(size(QC(i).GLMMASK,1),size(tempimg,2)); %Put back in volume space
    tempimg_out(logical(QC(i).GLMMASK),:) = tempimg; %Put back in volume space
    clear tempimg
    write_4dfpimg(tempimg_out,[ QC(i).vcnum '_' allends '.4dfp.img'],'bigendian');
    write_4dfpifh([ QC(i).vcnum '_' allends '.4dfp.img'],size(tempimg_out,2),'bigendian');
    
    % create a final DV file
    system(['cp total_DV_' LASTCONC{stage} '.txt total_DV_final.txt']);
    
    %%% CLEAN UP %%%
    if cleanup
        system(['rm *.conc*']);
        for j=1:size(QC(i).restruns,1)
            cd(QC(i).subbolddir{j});
            for k=2:stage % all intermediate files
                system(['rm ' LASTIMG{i,j,k} '.*' ]);
            end
        end
    end
    
    % write summary log file
    cd(QC(i).subdir);
    fid=fopen('processing_log.txt','w');
    fprintf(fid,'Subject: %s\n',QC(i).vcnum);
    fprintf(fid,'Sourcedir: %s\n',df.filepath{i});
    fprintf(fid,'Sourceprm: %s\n',df.prmfile{i});
    fprintf(fid,'set boldruns = (');
    for j=1:size(QC(i).restruns,1)
        fprintf(fid,'%d',cell2mat(QC(i).restruns(j,:)));
        if j~=size(QC(i).restruns,1)
            fprintf(fid,' ');
        end
    end
    fprintf(fid,')\n');
    fprintf(fid,'tmasktype: %s\n',tmasktype);
    if ~isequal(tmasktype,'ones')
        fprintf(fid,'%s',tm.tmaskfiles{i,1});
    end
    fprintf(fid,'switches.doregression (1=yes;0=no): %d\n',switches.doregression);
    fprintf(fid,'switches.regressiontype (0=classic;1=freesurfer;2=user4dfp;3:usertxt): %d\n',switches.regressiontype);
    fprintf(fid,'switches.motionestimates (0=no; 1=R,R`; 2=FRISTON; 20=R,R`,12rand): %d\n',switches.motionestimates);
    fprintf(fid,'switches.WM (1=regress;0=no): %d\n',switches.WM);
    fprintf(fid,'switches.V (1=regress;0=no): %d\n',switches.V);
    fprintf(fid,'switches.GS (1=regress;0=no): %d\n',switches.GS);
    fprintf(fid,'switches.dointerpolate (1=yes;0=no): %d\n',switches.dointerpolate);
    fprintf(fid,'switches.dobandpass (1=yes;0=no): %d\n',switches.dobandpass);
    fprintf(fid,'switches.temporalfiltertype (1=lowpass;2=hipass;3=bandpass): %d\n',switches.temporalfiltertype);
    fprintf(fid,'switches.lopasscutoff (in Hz; 0.08 is typical): %g\n',switches.lopasscutoff);
    fprintf(fid,'switches.hipasscutoff (in Hz; 0.009 is typical): %g\n',switches.hipasscutoff);
    fprintf(fid,'switches.order (1 is typical): %g\n',switches.order);
    fprintf(fid,'switches.doblur (1=yes;0=no): %d\n',switches.doblur);
    fprintf(fid,'switches.blurkernel (in mm; 6 is typical): %d\n',switches.blurkernel);
    fprintf(fid,'PROCESSING STEPS\n');
    for j=1:numel(LASTCONC)
        fprintf(fid,'%s\n',LASTCONC{j});
    end
    fclose(fid);
    
    % write a tiff of the QC figures for this subject
    close;
    QC(i).numTRs=numel(QC(i).FD);
    subplot(3,1,1);
    plot(QC(i).FD,'r');
    xlim([0 QC(i).numTRs]);
    ylim([0 1]); ylabel('FD');
    subplot(3,1,2); [figaxis]=plotyy(1:QC(i).numTRs,QC(i).DV_GLM(:,1),1:QC(i).numTRs,QC(i).DV_GLM(:,end));
    set(figaxis(1),'ylim',[0 50],'ytick',[0 50],'xlim',[0 QC(i).numTRs]);
    set(figaxis(2),'ylim',[0 5],'ytick',[0 5],'xlim',[0 QC(i).numTRs]);
    ylabel('DV');
    r1=corrcoef([QC(i).FD(QC(i).tmask) QC(i).DV_GLM(QC(i).tmask,1)]);
    r1=r1(2,1);
    r2=corrcoef([QC(i).FD(QC(i).tmask) QC(i).DV_GLM(QC(i).tmask,end)]);
    r2=r2(2,1);
    legend({['333 r=' num2str(r1,2) ],['final r=' num2str(r2,2)]});
    hold on;
    subplot(3,1,3); imagesc(QC(i).tmask');
    saveas(gcf,'total_QC.tiff','tiff');
    close;
    
end


%% RETURN TO THE STARTING DIRECTORY AND SAVE LABELS
cd(tdir);
% write a corrfile
fid=fopen('corrfile.txt','w');
for i=1:numsubs
    fprintf(fid,'%s\t%s\t%s\n',[QC(i).subdir '/' QC(i).vcnum '_' allends '.4dfp.img'],[QC(i).subdir '/total_tmask.txt'],QC(i).vcnum);
end
fclose(fid);
save('QC.mat','QC','-v7.3');
cd(startdir);

toc
%exit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bolds]=conc2bolds(concfile)

system(['cat ' concfile ' | grep file: | awk -F: ''{print $2}'' >! tmpAB']);
[bolds]=textread('tmpAB','%s');
system('rm tmpAB');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [fcimg]=bolds2img(bolds,trtot,trborders)
% vox=147456;
% fcimg=zeros(vox,trtot);
% for j=1:size(bolds,1)
%     fcimg(:,trborders(j,1):trborders(j,2))=read_4dfpimg_HCP(bolds{j,1});
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fcimg]=bolds2img(bolds,trtot,trborders,GLMmask)

vox = nnz(GLMmask);
%vox=902629;%147456;
fcimg=zeros(vox,trtot);
for j=1:size(bolds,1)
    temp = read_4dfpimg_HCP(bolds{j,1});
%    fcimg(:,trborders(j,1):trborders(j,2))=read_4dfpimg_HCP(bolds{j,1});
    fcimg(:,trborders(j,1):trborders(j,2))=temp(logical(GLMmask),:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rmstotal rmstrans rmsrot rmscol mvm] = rdat_calculations(datfile,varargin)

% jdp 9/15/10
% Here, you pass in an rdat, and get back out the RMS calculations for that file
% RMS is calculated for 3 rotational and 3 translational motion parameters, combined for rotation and translation, and a total RMS is also returned.
% The rdat is also returned in mm as the mvm variable, giving motion values in mm movement.
% The number of frames of a run to skip is set to a default of 5, and the radius for movement calculations is set to 50 (mm), which is standard. Those values can be altered by apssing in additional arguments
%
% USAGE: [rmstotal rmstrans rmsrot rmscol mvm] = rdat_calculations(datfile,*skipframes,radius*)
%

% set default skipframes and radius, but replace with user-defined values values if provided
radius=50;
skipframes=5;
if ~isempty(varargin)
    skipframes=varargin{1,1};
    radius=varargin{1,2};
end

% load the rdat
% fprintf('\n%s\n',datfile);
[pth,fname,ext]=filenamefinder(datfile,'dotsout');
outputfile=[ fname '_rmscalc.txt' ];

% have to strip out a header (has #s at the beginning of each line)
command=[ 'grep -v # ' datfile ];
[trash tempmat]=system(command);
datfile=str2num(tempmat);
if datfile(1,1)~=1
    error('Reading rdatfile isn''t getting a first frame of 1');
end
d=size(datfile);

% convert the rotational mvmt to mm movement
mvm=zeros(d(1),6);
mvm(:,1)=datfile(:,2);
mvm(:,2)=datfile(:,3);
mvm(:,3)=datfile(:,4);
mvm(:,4)=convert_deg_to_mm(datfile(:,5),radius);
mvm(:,5)=convert_deg_to_mm(datfile(:,6),radius);
mvm(:,6)=convert_deg_to_mm(datfile(:,7),radius);

startframe=skipframes+1;
meancol=mean(mvm(startframe:end,:),1);
stdcol=std(mvm(startframe:end,:),1);
[rmstotal rmscol]=calc_rms(mvm(startframe:end,1),mvm(startframe:end,2),mvm(startframe:end,3),mvm(startframe:end,4),mvm(startframe:end,5),mvm(startframe:end,6));
rmstrans=sqrt(sum(rmscol(1,1:3).^2));
rmsrot=sqrt(sum(rmscol(1,4:6).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deg] = convert_deg_to_mm(deg,radius)

deg=deg*(2*radius*pi/360);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totrms rms] = calc_rms(varargin)

% feed in arbitrary number of columns

d=size(varargin);
for i=1:d(2)
    vals(:,i)=varargin{1,i};
    rms(1,i)=sqrt(mean(vals(:,i).^2));
end

totrms=sqrt(sum(rms.^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempbold tempbetas] = demean_detrend(img,varargin)

if ~isnumeric(img)
    [tempbold]=read_4dfpimg_HCP(img); % read image
else
    [tempbold]=img;
    clear img;
end
[vox ts]=size(tempbold);

if ~isempty(varargin)
    tmask=varargin{1,1};
else
    tmask=ones(ts,1);
end

linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
tempboldcell=num2cell(tempbold(:,logical(tmask))',1);
linregcell=repmat({linreg(logical(tmask),:)},[1 vox]);
tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
tempbetas=cell2mat(tempbetas);
tempbetas=tempbetas';
tempintvals=tempbetas*linreg';
tempbold=tempbold-tempintvals;

if nargin==3
    outname=varargin{1,2};
    write_4dfpimg(tempbold,outname,'bigendian');
    write_4dfpifh(outname,size(tempbold,2),'bigendian');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tempimg zb newregs] = regress_nuisance(tempimg,totregs,tot_tmask)

[vox ts]=size(tempimg);
zlinreg=totregs(logical(tot_tmask),:); % only desired data
[zlinreg DMDTB]=demean_detrend(zlinreg'); % obtain fits for desired data
zlinreg=zlinreg';
zstd=std(zlinreg); % calculate std
zmean=mean(zlinreg);
zlinreg=(zlinreg-repmat(zmean,[size(zlinreg,1) 1]))./(repmat(zstd,[size(zlinreg,1) 1])); % zscore

linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
newregs=DMDTB*linreg'; % predicted all regressors demean/detrend
newregs=totregs-newregs'; % these are the demeaned detrended regressors
newregs=(newregs-repmat(zmean,[size(newregs,1) 1]))./(repmat(zstd,[size(newregs,1) 1])); % zscore

% now we have z-scored, detrended good and all regressors.

% demean and detrend the desired data
zmdtimg=tempimg(:,logical(tot_tmask));
[zmdtimg zmdtbetas]=demean_detrend(zmdtimg);

% calculate betas on the good data
tempboldcell=num2cell(zmdtimg',1);
zlinregcell=repmat({zlinreg},[1 vox]);
zb = cellfun(@mldivide,zlinregcell,tempboldcell,'uniformoutput',0);
zb=cell2mat(zb);

% demean and detrend all data using good fits
[zmdttotimg]=zmdtbetas*linreg';
zmdttotimg=tempimg-zmdttotimg;

% calculate residuals on all the data
zb=zb';
tempintvals=zb*newregs';
tempimg=zmdttotimg-tempintvals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tc] = gettcs(roimasks,tempimg)

d=size(roimasks);
for i=1:d(2) % cycle all masks
    tc(:,i)=(mean(tempimg(logical(roimasks(:,i)),:)))';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QC] = tempimgsignals(QC,i,tempimg,switches,stage)

difftempimg=diff(tempimg')';
difftempimg=[zeros(size(difftempimg,1),1) difftempimg];

QC(i).MEAN_CSF(:,stage)=mean(tempimg(QC(i).CSFMASK_glmmask,:))';
QC(i).MEAN_WM(:,stage)=mean(tempimg(QC(i).WMMASK_glmmask,:))';
QC(i).MEAN_GLM(:,stage)=mean(tempimg(QC(i).GLMMASK_glmmask,:))';

QC(i).SD_CSF(:,stage)=std(tempimg(QC(i).CSFMASK_glmmask,:))';
QC(i).SD_WM(:,stage)=std(tempimg(QC(i).WMMASK_glmmask,:))';
QC(i).SD_GLM(:,stage)=std(tempimg(QC(i).GLMMASK_glmmask,:))';

QC(i).SDbar_CSF(:,stage)=mean(QC(i).SD_CSF(:,stage));
QC(i).SDbar_WM(:,stage)=mean(QC(i).SD_WM(:,stage));
QC(i).SDbar_GLM(:,stage)=mean(QC(i).SD_GLM(:,stage));

[jnk tempDV]=array_calc_rms(difftempimg(QC(i).CSFMASK_glmmask,:));
tempDV(QC(i).runborders(:,2))=NaN;
QC(i).DV_CSF(:,stage)=tempDV';
[jnk tempDV]=array_calc_rms(difftempimg(QC(i).WMMASK_glmmask,:));
tempDV(QC(i).runborders(:,2))=NaN;
QC(i).DV_WM(:,stage)=tempDV';
[jnk tempDV]=array_calc_rms(difftempimg(QC(i).GLMMASK_glmmask,:));
tempDV(QC(i).runborders(:,2))=NaN;
QC(i).DV_GLM(:,stage)=tempDV';

QC(i).DVbar_CSF(stage)=nanmean(QC(i).DV_CSF(:,stage));
QC(i).DVbar_WM(stage)=nanmean(QC(i).DV_WM(:,stage));
QC(i).DVbar_GLM(stage)=nanmean(QC(i).DV_GLM(:,stage));

if switches.regressiontype==1 | switches.regressiontype==9
    QC(i).MEAN_GM(:,stage)=mean(tempimg(QC(i).GMMASK_glmmask,:))';
    QC(i).MEAN_WB(:,stage)=mean(tempimg(QC(i).WBMASK_glmmask,:))';
    
    QC(i).SD_GM(:,stage)=std(tempimg(QC(i).GMMASK_glmmask,:))';
    QC(i).SD_WB(:,stage)=std(tempimg(QC(i).WBMASK_glmmask,:))';
    
    QC(i).SDbar_GM(:,stage)=mean(QC(i).SD_GM(:,stage));
    QC(i).SDbar_WB(:,stage)=mean(QC(i).SD_WB(:,stage));
    
    [jnk tempDV]=array_calc_rms(difftempimg(QC(i).GMMASK_glmmask,:));
    tempDV(QC(i).runborders(:,2))=NaN;
    QC(i).DV_GM(:,stage)=tempDV';
    [jnk tempDV]=array_calc_rms(difftempimg(QC(i).WBMASK_glmmask,:));
    tempDV(QC(i).runborders(:,2))=NaN;
    QC(i).DV_WB(:,stage)=tempDV';
    
    QC(i).DVbar_GM(stage)=nanmean(QC(i).DV_GM(:,stage));
    QC(i).DVbar_WB(stage)=nanmean(QC(i).DV_WB(:,stage));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makepictures(QC,stage,switches,rightsignallim,leftsignallim,FDmult)

rylimz=[min(rightsignallim) max(rightsignallim)];
lylimz=[min(leftsignallim) max(leftsignallim)];

numpts=numel(QC.FD);
clf;
subplot(8,1,1:2);
pointindex=1:numpts;
[h a(1) a(2)]=plotyy(pointindex,QC.DV_GLM(:,stage),pointindex,QC.MEAN_GLM(:,stage));
set(h(1),'xlim',[0 numpts],'ylim',lylimz,'ycolor',[0 0 0],'ytick',leftsignallim);
ylabel('B:DV G:GLM');
set(a(1),'color',[0 0 1]);
set(h(2),'xlim',[0 numpts],'ycolor',[0 0 0],'xlim',[0 numel(QC.DV_GLM(:,stage))],'ylim',rylimz,'ytick',rightsignallim);
set(a(2),'color',[0 .5 0]);
axis(h(1)); hold on;
h3=plot([1:numpts],QC.FD*FDmult,'r');
hold off;
axes(h(2)); hold(h(2),'on');
hline(0,'k');
hold(h(2),'off');
set(h(1),'children',[a(1) h3]);
set(h(1),'YaxisLocation','left','box','off');
set(h(2),'xaxislocation','top','xticklabel','');
if switches.regressiontype==1 | switches.regressiontype==9
    subplot(8,1,3);
    imagesc(QC.MVM',[-.5 .5]);
    ylabel('XYZPYR');
    subplot(8,1,4:7);
    imagesc(QC.GMtcs(:,:,stage),rylimz); colormap(gray); ylabel('GRAY');
    subplot(8,1,8);
    imagesc([QC.WMtcs(:,:,stage);QC.CSFtcs(:,:,stage)],rylimz); ylabel('WM CSF');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,f,s,c,tau,w] = getTransform(t,h,TH,Tr,ofac,hifac)

%Input t is a column vector listing the time points for which observations
%are present.  Input h is a matrix with observations in columns and the
%number of rows equals the number the time points.  For our purposes number
%of voxels = number of columns.  Ofac = oversampling frequency (generally
%>=4), hifac = highest frequency allowed.  hifac = 1 means 1*nyquist limit
%is highest frequency sampled.  
%Lasted edited:  Anish Mitra, October 25 2012

N = size(h,1); %Number of time points
T = max(t) - min(t); %Total time span

%calculate sampling frequencies
f = (1/(T*ofac):1/(T*ofac):hifac*N/(2*T)).';

%angular frequencies and constant offsets
w = 2*pi*f;
tau = atan2(sum(sin(2*w*t.'),2),sum(cos(2*w*t.'),2))./(2*w);

%spectral power sin and cosine terms
cterm = cos(w*t.' - repmat(w.*tau,1,length(t)));
sterm = sin(w*t.' - repmat(w.*tau,1,length(t)));

num_tseries = size(h,2); %Number of time series

% D = bsxfun(@minus,h,mean(h));  %This line involved in normalization
D = h;
D = reshape(D,1,N,num_tseries);


%% C_final = (sum(Cmult,2).^2)./sum(Cterm.^2,2);
% This calculation is done by separately for the numerator, denominator,
% and the division

Cmult = bsxfun(@times, cterm,D);
%rewrite the above line with bsxfun to optimize further?
%bsxfun(@power,sum(Cmult,2),2) = sum(Cmult.^2,2) = numerator
% numerator = bsxfun(@power,sum(Cmult,2),2);
numerator = sum(Cmult,2); %Modify the numerator to get the cw term

%sum(bsxfun(@power,Cterm,2),2) = sum(Cterm.^2,2) = denominator
denominator = sum(bsxfun(@power,cterm,2),2); %use Cterm in place of cterm to make it exactly the denominator in the original expression
C_final_new = bsxfun(@rdivide,numerator,denominator);
c = C_final_new;
clear numerator denominator cterm Cmult C_final_new

%% Repeat the above for Sine term
Smult = bsxfun(@times, sterm,D);
% S_final = (sum(Smult,2).^2)./sum(Sterm.^2,2);
% numerator = bsxfun(@power,sum(Smult,2),2);
numerator = sum(Smult,2); %Modify the numerator to get the sw term
denominator = sum(bsxfun(@power,sterm,2),2);
S_final_new = bsxfun(@rdivide,numerator,denominator);
s = S_final_new;
clear numerator denominator sterm Smult S_final_new

% %% Power = C_final + S_final; 
% Power = C_final_new + S_final_new;
% Power_reshaped = reshape(Power,size(Power,1),num_tseries);
% Power_final = bsxfun(@rdivide,Power_reshaped,2*var(h)); %Normalize the power
% clearvars -except Power_final f


% The inverse function to re-construct the original time series
Time = TH';
T_rep = repmat(Time,[size(f,1),1,size(h,2)]);
% T_rep = bsxfun(@minus,T_rep,tau);
% s(200) = s(199);
w = 2*pi*f;
prod = bsxfun(@times,T_rep,w);
sin_t = sin(prod);
cos_t = cos(prod);
sw_p = bsxfun(@times,sin_t,s);
cw_p = bsxfun(@times,cos_t,c);
S = sum(sw_p);
C = sum(cw_p);
H = C + S;
H = reshape(H,size(Time,2),size(h,2));

%Normalize the reconstructed spectrum, needed when ofac > 1
Std_H = std(H);
Std_h = std(h);
norm_fac = Std_H./Std_h;
H = bsxfun(@rdivide,H,norm_fac);




