function [Ci] = run_infomap_on_pajekfile(pajekfilename,reps)
%[Ci] = run_infomap_on_pajekfile(pajekfilename,reps)
%
%
% This script runs infomap on a pajekfile with some number of
% repetitions. It then returns the community assignments found.

infomapfolder = '/data/cn/data1/scripts/CIFTI_RELATED/Infomap';

% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);


% run infomap
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
[failed, message] = system([infomapfolder '/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr]);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));
if logical(failed)
    disp(message)
end



% So parfor doesn't crap out
isclufile = exist(clufile);
while isclufile == 0
    pause(60)
    isclufile = exist(clufile);
end

Ci = textread(clufile,'%d','headerlines',1);
end
