function [clrs] =  regularize(clrs,varargin)
%[clrs] =  regularize(clrs,varargin)
% 
% This script takes a matrix (typically of clustering assignments) and
% organizes it such that values in each column represent the same assignment.
% The output is returned as a variable, and can be written to a file as
% well if a second argument is passed in.
% 
% USAGE: [newclrs] = regularize(clrfile,*outputfile*)
% USAGE: [newclrs] = regularize(clrmatrix)
%
% EMG 08/31/15, adapted from 'rawoutput2clr', by jdp 10/10/10

% assume matrix passed in, otherwise load the matrix
if ~isnumeric(clrs)
    clrs=load(clrs);
end

% if zeros are in here, give then a new number (b/c zeros are special for other scripts in graphtools)
[c d]=size(clrs);
maxclr=max(max(clrs));
if nnz(ismember(clrs,0))
    disp 'Setting zeros to something else';
    clrs(clrs==0)=maxclr+1;
    maxclr=maxclr+1;
end

% zero out the new color matrix and assign its first column
newclrs=zeros(c,d);
newclrs(:,1)=clrs(:,1);

% now cycle through the other columns
for j=2:d
    %j
    % get unique colors of adjacent columns
    clear firstuns seconduns
    firstuns=unique(newclrs(:,j-1)); % from newly assigned colors
    seconduns=unique(clrs(:,j)); % from existing colors
    
    % make a matrix of the numbers of overlaps between pairs of unique
    % values in adjacent columns
    oerlaps=zeros(size(firstuns,1),size(seconduns,1));
    for x=1:size(firstuns,1)
        for y=1:size(seconduns,1)
            oerlaps(x,y)=nnz((ismember(newclrs(:,j-1),firstuns(x))) & (ismember(clrs(:,j),seconduns(y))));
        end
    end
    
    % get sorted overlaps
    clear olap
    [olap(:,1) olap(:,2) olap(:,3)]=find(oerlaps);
    [cc dd] = sortrows(olap,-3);

    
    % assign row2 to row1 values when row1 unused
    firstused=zeros(size(firstuns,1),1);
    secondused=zeros(size(seconduns,1),1);
    while ~isempty(cc)

        clear exitnum insertnum;
        
        % if both old and new colors in question haven't been assigned
        if (firstused(cc(1,1))==0 && secondused(cc(1,2))==0)
            
            % now they're assigned
            firstused(cc(1,1))=1;
            secondused(cc(1,2))=1;
            
            % oldclrs with this value will be reassigned
            exitnum=seconduns(cc(1,2));
            % to this value from the newclrs
            insertnum=firstuns(cc(1,1));
            % and here they are assigned
            newclrs((clrs(:,j)==exitnum),j)=insertnum;
            
            % this deletes all entries of this oldclr from the unique list
            cc=cc(~(cc(:,2)==cc(1,2)),:);

        
        % if an oldclr can't find parthers in newclrs (the newclr has
        % already been given away, probably)
        elseif (firstused(cc(1,1))==1 && secondused(cc(1,2))==0)
            
            % now this color is used from oldclrs
            secondused(cc(1,2))=1;
            % and we will be replacing oldclrs with this value
            exitnum=seconduns(cc(1,2));
            % with the highest number possible
            newclrs((clrs(:,j)==exitnum),j)=maxclr;
            % and then we'll bump the highest number possible up one more
            maxclr=maxclr+1;
            % and now we'll delete entries of this oldclr from the list
            cc=cc(~(cc(:,2)==cc(1,2)),:);
        
        % the case before this case should take care of all scenarios, but
        % this is a catch in case somehow a newclr didn't find a partner in
        % oldclrs, and the oldclr has already been used
        elseif (firstused(cc(1,1))==0 && secondused(cc(1,2))==1)
            disp newold! how did i get in here?
            cc=cc(~(cc(:,2)==cc(1,2)),:);
        end
    end
end

% now drop all values to the minimum possible
vals=unique(newclrs);
clrs=zeros(c,d);
for i=1:size(vals,1)
    clrs(newclrs==vals(i))=i;
end
clrs=single(clrs);
clrs(clrs <= 1) = 0;
clrs = clrs -1;



% set the output file name if provided
if ~isempty(varargin)
    outfile=varargin{1,1};
    dlmwrite(outfile,clrs,'delimiter',' ');
end
