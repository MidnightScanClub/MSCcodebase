function [pth fbase ext] = filenamefinder(filename,dots);

% this takes a file, and returns the path, name, and extension.
% you choose dotsin or dotsout to keep or remove dots in filenames

% jdp 6/2010

[pth,fbase,ext]=fileparts(filename);
if isempty(pth)
    pth=pwd;
end

switch dots
    case 'dotsin'
    case 'dotsout'
        fbase=regexprep(fbase,'\.','');
        pth=regexprep(pth,'\.','');
    otherwise
        error('Use dotsin or dotsout');
end