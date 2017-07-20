function [clrs]=modify_clrfile(type,clrfile,varargin)
%
% Name:modify_clrfile.m
% $Revision: 1.2 $
% $Date: 2011/03/08 20:30:23 $
%
% jdp 10/10/10
% 
% This does a lot of manipulations of clrfiles, and sometimes associated
% files. This script is interactive and will prompt you for choices and
% information as needed.
% 
% Output is saved and named sensibly, depending on what you did.
% 
% USAGE: [clrs]=modify_clrfile(operationtype,clrfile,*variableinput*)
% 
% At the moment, operations it performs include:
%   sort: sorts a clrfile and rearranges an associated .mat and .roi file
%   resort: sorts a naive dataset using the index from another dataset
%   singlecolumn: pulls a single column from a prmfile
%   trim: pulls a set of columns from a clrfile
%   simplify: at every threshold sets all modules with <x members to -1 (a null module)
%   color: puts up the color chart and you pick the rgbs to make custom rgb
%   fixcolor: pulls up an RGB list, shows it, and lets you fix it 
% 
% USAGE: [clrs]=modify_clrfile(operationtype,clrfile,*variableinput*)
% USAGE: [clrs]=modify_clrfile('singlecolumn',clrfile,column)
%   -> columnclrfile
% USAGE: [clrs]=modify_clrfile('trim',clrfile,columnA,columnZ)
%   -> trimmedclrfile
% USAGE: [clrs]=modify_clrfile('simplify',clrfile,mininummodulesize)
%   -> clrfile with modules under minimum assigned to -1 'junk' module
% USAGE: [clrs]=modify_clrfile('sort',clrfile,prmfile,sortbythesecolumns)
%   -> sorted clrfile, roifile, and matfile and new prmfile
%   -> index (feed with resort switch to reorder other datasets)
% USAGE: [clrs]=modify_clrfile('resort',clrfile,prmfile,index)
%   -> sorted clrfile, roifile, and matfile and new prmfile
% USAGE: [clrs]=modify_clrfile('color',clrfile)
%   -> rgbfile
% USAGE: [clrs]=modify_clrfile('colorfix',rgbfile)
%   -> fixed rgbfile
% 
% USAGE: [clrs]=modify_clrfile('singlecolumn','clr.txt',2)
% USAGE: [clrs]=modify_clrfile('trim','clr.txt',1,3)
% USAGE: [clrs]=modify_clrfile('simplify','clr.txt',20)
% USAGE: [clrs]=modify_clrfile('sort','clr.txt','modbox.prm',[5:15])
%       * this would sort using columns 5-15 of the data, sequentially
% USAGE: [clrs]=modify_clrfile('color','clr.txt')
% USAGE: [clrs]=modify_clrfile('colorfix',rgbfile)
% 
% NOTES: 11/9/10 set 'simplify' junk module to -1 rather than 0

close;

% get the basics of the clrfile
clrs=load(clrfile);
[a b]=size(clrs);
[pth,fname,ext]=filenamefinder(clrfile,'dotsout');

switch type
    %%%%%%%%%%%%%%%%%%%
    case 'singlecolumn' % excise a single column
    %%%%%%%%%%%%%%%%%%%
    
        column=varargin{1,1};
        fname=[ pth '/' fname '_col' num2str(column) '.txt' ];
        clrs=clrs(:,column);
        dlmwrite(fname,clrs,'\t');
    
    %%%%%%%%%%%
    case 'trim' % excise a set of columns
    %%%%%%%%%%%
        
        columnA=varargin{1,1};
        columnZ=varargin{1,2};
        fname=[ pth '/' fname '_col' num2str(columnA) 'to' num2str(columnZ) '.txt' ];
        clrs=clrs(:,columnA:columnZ);
        dlmwrite(fname,clrs,'\t');
    
    %%%%%%%%%%%%%%%
    case 'simplify' % at each threshold, set modules with <X members to 0
    %%%%%%%%%%%%%%%    
        
        minsize=varargin{1,1};
        subplot(1,2,1);
        imagesc(clrs);
        olduniques=unique(clrs);
        for j=1:b
            clear uniques;
            uniques(:,1)=unique(clrs(:,j));
            for i=1:size(uniques,1)
                if nnz(clrs(:,j)==uniques(i,1)) < minsize
                    clrs((clrs(:,j)==uniques(i,1)),j)=-1;
                end
            end
        end
        subplot(1,2,2);
        imagesc(clrs);
        fname=[pth '/' fname '_minsize' num2str(minsize) '.txt' ];
        dlmwrite(fname,clrs,'\t');
        newuniques=unique(clrs);
        fprintf('Modules reduced from %d to %d\n',size(olduniques,1),size(newuniques,1));
    
    %%%%%%%%%%%
    case 'sort'
    %%%%%%%%%%%
    
        prmfile=varargin{1,1};
        sortparams=varargin{1,2};
        subplot(1,2,1);
        imagesc(clrs);
                
        % sort the clrs using the sortparams
        [clrs index] =sortrows(clrs,sortparams);
        
        % write the sorted colors
        dlmwrite([pth '/' fname '_sorted.txt' ],clrs,'\t');
        
        % write the index for the sorting (use with the 'resort' switch)
        dlmwrite([pth '/' fname '_sorted_index.txt' ],index,'\t');
        
        % show the clrs
        subplot(1,2,2);
        imagesc(clrs);
        
        % create new roifile, matfile, and prmfile based on the sorting
        reindex_prmfile(prmfile,index,'_sorted');
        
    %%%%%%%%%%%%
    case 'resort'
    %%%%%%%%%%%%
        
        % assign the variables
        prmfile=varargin{1,1};
        index=varargin{1,2};
        
        % if the index isn't passed in as an array, presume it's a text file
        % that should be read in (e.g. the output of 'sort')
        if ~isnumeric(index)
            index=load(index);
        end
        
        reindex_prmfile(prmfile,index,'_resorted',clrfile);
        
 
    %%%%%%%%%%%%    
    case 'color'
    %%%%%%%%%%%%
    
        
        
        % show the original file
        close;
        subplot(1,3,1);
        imagesc(clrs);
        
        % how many colors to assign
        uniques=unique(clrs);
        numuniques=size(uniques,1);
        rgb=zeros(a,b,3);
        rgblist=[0 0 0];
        
        %for each module
        for i=1:numuniques
            repeating=1;
            fprintf('Filling module %d (%d/%d)\n',uniques(i),i,numuniques);
            while repeating
                holdrgb=rgb; % store the prior rgb in case they dislike the new one
                
                % display the current rgb
                fprintf('rgblist is currently:\n');
                rgblist
                
                % show them the module they're targeting
                mask=ismember(clrs,uniques(i));
                subplot(1,3,2);
                imagesc(mask);
                
                % enter a color for the module
                rgblist(i,1)=input('Enter r (0-255) ');
                rgblist(i,2)=input('Enter r (0-255) ');
                rgblist(i,3)=input('Enter r (0-255) ');
                
                % show those colors on the module
                rgb(:,:,1)=rgb(:,:,1)+mask.*rgblist(i,1);
                rgb(:,:,2)=rgb(:,:,2)+mask.*rgblist(i,2);
                rgb(:,:,3)=rgb(:,:,3)+mask.*rgblist(i,3);
                subplot(1,3,3);
                imshow(rgb/255);
                
                requestinginput=1;
                while requestinginput
                    repeat=input('Look ok? (y/n) ','s');
                    switch repeat
                        case 'n'
                            repeating=1;
                            requestinginput=0;
                            rgb=holdrgb; % restore prior rgb
                        case 'y'
                            repeating=0;
                            requestinginput=0;
                        otherwise
                            disp('Enter y or n');
                    end
                end
            end
        end
        
        fileend=input('Enter a string to tag this rgbfile with: ','s');
        [pth,fname,ext]=filenamefinder(clrfile,'dotsout');
        rgbfile=[ pth '/' fname '_' fileend '_RGB.txt' ];
        dlmwrite(rgbfile,rgblist,'\t');
        
        jpgfile=[ pth '/' fname '_' fileend 'RGB.tiff' ];
        imwrite(rgb/255,jpgfile);
    
    %%%%%%%%%%%%%%%
    case 'colorfix'
    %%%%%%%%%%%%%%%
    
        % load the rgbfile and show people their options
        a=load(clrfile); % here clrfile is an rgbmat from color or the rgbmapper scripts
        fprintf('There are %d colors here\n',size(a,1));
        b=reshape(a,[size(a,1) 1 size(a,2)]); % shift it into the image format
        subplot(1,2,1);
        imshow(b/255);
        
        keepfixing=1;
        while keepfixing
            
            % alter one of the colors and show the fix
            line=input('Which line should we fix? ');
            a(line,1)=input('Enter r (0-255) ');
            a(line,2)=input('Enter r (0-255) ');
            a(line,3)=input('Enter r (0-255) ');
            subplot(1,2,2);
            b=reshape(a,[size(a,1) 1 size(a,2)]);
            imshow(b/255);
            
            % make sure they're happy with the fix
            doneyet=input('Is that all? (y/n) ','s');
            switch doneyet
                case 'y'
                    keepfixing=0;
            end
        end
        
        % save the alterations
        fileend=input('Enter a string to tag this rgbfile with: ','s');
        [pth,fname,ext]=filenamefinder(clrfile,'dotsout');
        rgbfile=[ pth '/' fname fileend '.txt' ];
        dlmwrite(rgbfile,a,'\t');
    
    %%%%%%%%%%%%
    case 'split'
    %%%%%%%%%%%%
    
        % this splits a network by module, making little roifiles and
        % matrices for the modules.
        
        % set the variables
        column=varargin{1,1};
        prmfile=varargin{1,2};
                
        % read in the prmfile, etc
        [matfile roifile stemname subjectA subjectZ loend step hiend writepath threshold boxcarsize boxcarstep thresholdtype] = prmfilereader(prmfile);
        
        % find which analysis we're interested in
        [subjectarray thresholdarray numanalyses xarray] = matrix_parameter_setter(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,'thr');
        
        % read in the roifile
        [xyz nodenames] = roifilereader(roifile);
        origindx=[1:size(xyz,1)]';
        
        % load the matrix for splitting
        [matrix] = matfile_loader(matfile);

        % make directory for output
        outdir=[ pth '/' fname '_col' num2str(column) ];
        if ~exist(outdir)
            mkdir(outdir);
        end
        
        % extract the column of interest, find the unique values
        clrs=clrs(:,column);
        uclrs=unique(clrs);
        numuclrs=size(uclrs,1);
        
        % split the data by module
        for i=1:numuclrs
            outbase=[ 'module' num2str(uclrs(i)) ];
            indx=ismember(clrs,uclrs(i));
            quickroifile(xyz(indx,:),[outdir '/' outbase '.roi' ]);
            dlmwrite([outdir '/' outbase '_originindex.txt' ],origindx(indx));
            dlmwrite([outdir '/' outbase '_dummyclr.clr' ],ones(nnz(indx),1));
            clear modmat; 
            modmat=matrix(indx,indx,:);
            save([outdir '/' outbase '.mat' ],'modmat');
            quickprmfile([outdir '/' outbase '.prm' ],[outbase '.mat' ],[outbase '.roi' ],outbase,subjectarray(column,1),subjectarray(column,2),thresholdarray(column),thresholdarray(column),thresholdarray(column),outdir,threshold,boxcarsize,boxcarstep,thresholdtype);
            
        end
        
        fprintf('Output prmfiles are set up assuming you were using ''thr''. If you''re using ''box'' you will need to reconfigure them.\n');

    otherwise
        error('Type help modify_clrfile');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reindex_prmfile(prmfile,index,filetag,varargin)

% This function takes and index, and reorders matrices, roifiles, and
% clrfiles to match the indexing. This is useful for ordering cohorts, etc.
% 
% If the user passes in an extra variable, this is a clrfile that needs
% reordering as well.

% read the parameters
[matfile roifile stemname subjectA subjectZ loend step hiend writepath threshold boxcarsize boxcarstep thresholdtype] = prmfilereader(prmfile);

% resort the matrix
[pth,fname,ext]=filenamefinder(matfile,'dotsout');
newmatfile = [ fname filetag '.mat' ];
[oldmat] = matfile_loader(matfile);
[a b c]=size(oldmat);
newmat=zeros(a,b,c);
for k=1:c
    for i=1:a
        for j=1:a
            newmat(i,j,k)=oldmat(index(i),index(j),k);
        end
    end
end
save([pth '/' newmatfile],'newmat');

% read in the roifile
fid=fopen(roifile,'r');
[C]=textscan(fid,'%d%d%d%d%s%s%s%s%s','HeaderLines',1);
fclose(fid);
num=single(C{1,1});
X=single(C{1,2});
Y=single(C{1,3});
Z=single(C{1,4});
name=char(C{1,5});
dummy1=char(C{1,6});
dummy2=char(C{1,7});
dummy3=char(C{1,8});
dummy4=char(C{1,9});

% write a newly sorted ROIfile
[pth,fname,ext]=filenamefinder(roifile,'dotsout');
newroifile = [ fname filetag '.roi' ];
fid=fopen([ pth '/' newroifile ],'w');
fprintf(fid,'ROI\tX\tY\tZ\tName\tAnatomy\tAcolor\tBname\tBcolor\n');
for i=1:a
    fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n',i,X(index(i)),Y(index(i)),Z(index(i)),name(index(i),:),dummy1(index(i),:),dummy2(index(i),:),dummy3(index(i),:),dummy4(index(i),:));
end
fclose(fid);

% make new prmfile
[pth,fname,ext]=filenamefinder(prmfile,'dotsout');
newprmfile = [ fname filetag ext ];
quickprmfile(newprmfile,newmatfile,newroifile,[stemname filetag],subjectA,subjectZ,loend,step,hiend,writepath,threshold,boxcarsize,boxcarstep,thresholdtype);

% if a clrfile needs to be reordered also (e.g. with 'resort')
if ~isempty(varargin)
    clrfile=varargin{1,1};
    
    % use can pass in as txtfile or array. should almost always be txtfile
    if ~isnumeric(clrfile)
        clrs=load(clrfile);
        [pth,fname,ext]=filenamefinder(clrfile,'dotsout');
    else
        clrs=clrfile;
        fname='sort';
    end
    
    % create the indexed clrs
    newclrs=clrs(index,:);
    
    % write the indexed clrs
    dlmwrite([pth '/' fname filetag ext],newclrs,'\t');
end

    
