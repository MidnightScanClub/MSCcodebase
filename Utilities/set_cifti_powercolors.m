function set_cifti_powercolors(filename)


%insertstring = sprintf('\n\t<MetaData>\n\t\t<MD>\n\t\t\t<Name>PaletteColorMapping</Name>\n\t\t\t<Value>&lt;PaletteColorMapping Version=&quot;1&quot;&gt;\n\t&lt;ScaleMode&gt;MODE_USER_SCALE&lt;/ScaleMode&gt;\n\t&lt;AutoScalePercentageValues&gt;98.000000 2.000000 2.000000 98.000000&lt;/AutoScalePercentageValues&gt;\n\t&lt;AutoScaleAbsolutePercentageValues&gt;2.000000 98.000000&lt;/AutoScaleAbsolutePercentageValues&gt;\n\t&lt;UserScaleValues&gt;-100.000000 0.000000 1.000000 18.000000&lt;/UserScaleValues&gt;\n\t&lt;PaletteName&gt;power_surf&lt;/PaletteName&gt;\n\t&lt;InterpolatePalette&gt;true&lt;/InterpolatePalette&gt;\n\t&lt;DisplayPositiveData&gt;true&lt;/DisplayPositiveData&gt;\n\t&lt;DisplayZeroData&gt;false&lt;/DisplayZeroData&gt;\n\t&lt;DisplayNegativeData&gt;false&lt;/DisplayNegativeData&gt;\n\t&lt;ThresholdTest&gt;THRESHOLD_TEST_SHOW_OUTSIDE&lt;/ThresholdTest&gt;\n\t&lt;ThresholdType&gt;THRESHOLD_TYPE_OFF&lt;/ThresholdType&gt;\n\t&lt;ThresholdFailureInGreen&gt;false&lt;/ThresholdFailureInGreen&gt;\n\t&lt;ThresholdNormalValues&gt;-1.000000 1.000000&lt;/ThresholdNormalValues&gt;\n\t&lt;ThresholdMappedValues&gt;-1.000000 1.000000&lt;/ThresholdMappedValues&gt;\n\t&lt;ThresholdMappedAvgAreaValues&gt;-1.000000 1.000000&lt;/ThresholdMappedAvgAreaValues&gt;\n\t&lt;ThresholdDataName&gt;&lt;/ThresholdDataName&gt;\n\t&lt;ThresholdRangeMode&gt;PALETTE_THRESHOLD_RANGE_MODE_MAP&lt;/ThresholdRangeMode&gt;\n\t&lt;ThresholdLowHighLinked&gt;false&lt;/ThresholdLowHighLinked&gt;\n&lt;/PaletteColorMapping&gt;\n</Value>\n\t\t</MD>\n\t</MetaData>');
insertstring = sprintf('\n\t<MetaData>\n\t\t<MD>\n\t\t\t<Name>PaletteColorMapping</Name>\n\t\t\t<Value>&lt;PaletteColorMapping Version=&quot;1&quot;&gt;\n\t&lt;ScaleMode&gt;MODE_USER_SCALE&lt;/ScaleMode&gt;\n\t&lt;AutoScalePercentageValues&gt;98.000000 2.000000 2.000000 98.000000&lt;/AutoScalePercentageValues&gt;\n\t&lt;UserScaleValues&gt;-100.000000 0.000000 1.000000 18.000000&lt;/UserScaleValues&gt;\n\t&lt;PaletteName&gt;power_surf&lt;/PaletteName&gt;\n\t&lt;InterpolatePalette&gt;true&lt;/InterpolatePalette&gt;\n\t&lt;DisplayPositiveData&gt;true&lt;/DisplayPositiveData&gt;\n\t&lt;DisplayZeroData&gt;false&lt;/DisplayZeroData&gt;\n\t&lt;DisplayNegativeData&gt;false&lt;/DisplayNegativeData&gt;\n\t&lt;ThresholdTest&gt;THRESHOLD_TEST_SHOW_OUTSIDE&lt;/ThresholdTest&gt;\n\t&lt;ThresholdType&gt;THRESHOLD_TYPE_OFF&lt;/ThresholdType&gt;\n\t&lt;ThresholdFailureInGreen&gt;false&lt;/ThresholdFailureInGreen&gt;\n\t&lt;ThresholdNormalValues&gt;-1.000000 1.000000&lt;/ThresholdNormalValues&gt;\n\t&lt;ThresholdMappedValues&gt;-1.000000 1.000000&lt;/ThresholdMappedValues&gt;\n\t&lt;ThresholdMappedAvgAreaValues&gt;-1.000000 1.000000&lt;/ThresholdMappedAvgAreaValues&gt;\n\t&lt;ThresholdDataName&gt;&lt;/ThresholdDataName&gt;\n\t&lt;ThresholdRangeMode&gt;PALETTE_THRESHOLD_RANGE_MODE_MAP&lt;/ThresholdRangeMode&gt;\n&lt;/PaletteColorMapping&gt;\n</Value>\n\t\t</MD>\n\t</MetaData>');
rightafterstrings = {'<Matrix>','</MapName>'};


data = ft_read_cifti_mod(filename);

dat    = data.data;
dimord = data.dimord;

clear data

switch dimord
  case {'pos' 'pos_scalar' 'scalar_pos'}
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS
    extension   = '.dscalar.nii';
    intent_code = 3006;
    intent_name = 'ConnDenseScalar';
    dat = transpose(dat);
    dimord = 'scalar_pos';
%   case {'scalar_pos'}
%     % NIFTI_INTENT_CONNECTIVITY_DENSE_SCALARS
%     extension   = '.dscalar.nii';
%     intent_code = 3006;
%     intent_name = 'ConnDenseScalar';
  case 'pos_pos'
    % NIFTI_INTENT_CONNECTIVITY_DENSE
    extension = '.dconn.nii';
    intent_code = 3001;
    intent_name = 'ConnDense';
  case 'pos_time'
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES
    extension = '.dtseries.nii';
    intent_code = 3002;
    intent_name = 'ConnDenseSeries';
    dat = transpose(dat);
    dimord = 'time_pos';
  case 'pos_freq'
    % NIFTI_INTENT_CONNECTIVITY_DENSE_SERIES
    extension = '.dtseries.nii';
    intent_code = 3002;
    intent_name = 'ConnDenseSeries';
    dat = transpose(dat);
    dimord = 'freq_pos';
    
  case {'chan' 'chan_scalar'}
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SCALARS
    extension   = '.pscalar.nii';
    intent_code = 3006;
    intent_name = 'ConnParcelScalr'; % due to length constraints of the NIfTI header field, the last ?a? is removed
    dat = transpose(dat);
    dimord = 'scalar_chan';
  case 'chan_chan'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED
    extension = '.pconn.nii';
    intent_code = 3003;
    intent_name = 'ConnParcels';
  case 'chan_time'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SERIES
    extension = '.ptseries.nii';
    intent_code = 3000;
    intent_name = 'ConnParcelSries'; % due to length constraints of the NIfTI header field, the first "e" is removed
    dat = transpose(dat);
    dimord = 'time_chan';
  case 'chan_freq'
    % NIFTI_INTENT_CONNECTIVITY_PARCELLATED_SERIES
    extension = '.ptseries.nii';
    intent_code = 3000;
    intent_name = 'ConnParcelSries'; % due to length constraints of the NIfTI header field, the first "e" is removed
    dat = transpose(dat);
    dimord = 'freq_chan';
    
  otherwise
    error('unsupported dimord "%s"', dimord);
end % switch



% read the header section
hdr = read_nifti2_hdr(filename);

% xml_offset = 540+12;
% xml_size   = hdr.vox_offset-xml_offset-8;

fid = fopen(filename, 'rb', hdr.endian);

% determine the file size, this is used to catch endian errors
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

% set the default for readdata
%if isempty(readdata)
%   if filesize>1e9
%     warning('Not reading data by default in case filesize>1GB. Please specify the ''readdata'' option.');
%     readdata = false;
%   else
%    readdata = true;
%  end
%end

fseek(fid, 540, 'bof');
hdrext = fread(fid, [1 4], 'int8');
if hdrext(1)~=1
  error('cifti requires a header extension');
end

% determine the size of the header extension
esize = double(fread(fid, 1, 'int32=>int32'));
etype = fread(fid, 1, 'int32=>int32');

hdrsize = 540;
voxsize = filesize-double(hdr.vox_offset);
if esize>(filesize-hdrsize-voxsize)
  warning('the endianness of the header extension is inconsistent with the nifti-2 header');
  esize = swapbytes(esize);
  etype = swapbytes(etype);
end

if etype~=32 && etype~=swapbytes(int32(32)) % FIXME there is an endian problem
  error('the header extension type is not cifti');
end

% read the extension content, subtract the 8 bytes from esize and etype
xmldat = fread(fid, [1 esize-8], 'uint8=>char');

% the size of the extension must be an integer multiple of 16 bytes according to http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/extension.html
% consequently there might be some zero-padding at the end of the XML section
if any(xmldat==0)
  xmldat = xmldat(xmldat>0);
end

xmldat_orig = xmldat;

for i = 1:length(rightafterstrings)
    xmldat = strrep(xmldat, rightafterstrings{i}, [rightafterstrings{i} insertstring]);
end

whitespace = false(size(xmldat));

gt = int8('>');
lt = int8('<');
ws = int8(sprintf(' \t\r\n'));

b = find(xmldat==gt);
e = find(xmldat==lt);
b = b(1:end-1); % the XML section ends with ">", this is not of relevance
e = e(2:end);   % the XML section starts with "<", this is not of relevance
b = b+1;
e = e-1;

for i=1:length(b)
  for j=b(i):1:e(i)
    if any(ws==xmldat(j))
      whitespace(j) = true;
    else
      break
    end
  end
end

for i=1:length(b)
  for j=e(i):-1:b(i)
    if any(ws==xmldat(j))
      whitespace(j) = true;
    else
      break
    end
  end
end

% keep it if there is _only_ whitespace between the  ">" and "<"
for i=1:length(b)
  if all(whitespace(b(i):e(i)))
    whitespace(b(i):e(i)) = false;
  end
end

% remove the padding whitespace
xmldat  = xmldat(~whitespace);

xmlsize = length(xmldat);
xmlpad  = ceil((xmlsize+8)/16)*16 - (xmlsize+8);


hdr.intent_p1       = 0;
hdr.intent_p2       = 0;
hdr.intent_p3       = 0;
hdr.pixdim          = [0 1 1 1 1 1 1 1];
hdr.vox_offset      = 4+540+8+xmlsize+xmlpad;
hdr.scl_slope       = 1; % WorkBench sets scl_slope/scl_inter to 1 and 0, although 0 and 0 would also be fine - both mean the same thing according to the nifti spec
hdr.scl_inter       = 0;
hdr.cal_max         = 0;
hdr.cal_min         = 0;
hdr.slice_duration  = 0;
hdr.toffset         = 0;
hdr.slice_start     = 0;
hdr.slice_end       = 0;
hdr.descrip         = char(zeros(1,80));
hdr.aux_file        = char(zeros(1,24));
hdr.qform_code      = 0;
hdr.sform_code      = 0;
hdr.quatern_b       = 0;
hdr.quatern_c       = 0;
hdr.quatern_d       = 0;
hdr.qoffset_x       = 0;
hdr.qoffset_y       = 0;
hdr.qOffset_z       = 0;
hdr.srow_x          = [0 0 0 0];
hdr.srow_y          = [0 0 0 0];
hdr.srow_z          = [0 0 0 0];
hdr.slice_code      = 0;
hdr.xyzt_units      = 0;
hdr.intent_code     = intent_code;
hdr.intent_name     = cat(2, intent_name, zeros(1, 16-length(intent_name))); % zero-pad up to 16 characters
hdr.dim_info        = 0;
hdr.unused_str      = char(zeros(1,15));


% open the file
fid = fopen(filename, 'wb');

% write the header, this is 4+540 bytes
write_nifti2_hdr(fid, hdr);


% write the cifti header extension
fwrite(fid, [1 0 0 0], 'uint8');
fwrite(fid, 8+xmlsize+xmlpad, 'int32');   % esize
fwrite(fid, 32, 'int32');                 % etype
fwrite(fid, xmldat, 'char');              % write the ascii XML section
fwrite(fid, zeros(1,xmlpad), 'uint8');    % zero-pad to the next 16 byte boundary



% write the actual data
fwrite(fid, dat, 'single');

fclose(fid);



end


function [hdr] = read_nifti2_hdr(filename)

% READ_NIFTI2_HDR
%
% Use as
%   [hdr] = read_nifti2_hdr(filename)
% where
%   filename   = string
%   
% This implements the format as described at
%   http://www.nitrc.org/forum/forum.php?thread_id=2148&forum_id=1941
%
% Please note that it is different from the suggested format described here
%   http://www.nitrc.org/forum/forum.php?thread_id=2070&forum_id=1941
% and
%   https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2011-February/017482.html
% Notably, the unused fields have been removed and the size has been
% reduced from 560 to 540 bytes.
%
% See also WRITE_NIFTI_HDR, READ_CIFTI, WRITE_CIFTI

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

hdr.endian = 'l';
fid = fopen(filename, 'rb', hdr.endian);
hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'); % 0

if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
  % try opening as big endian
  fclose(fid);
  hdr.endian = 'b';
  fid = fopen(filename, 'r', 'b');
  hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'); % 0
end

if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
  fclose(fid);
  error('cannot open %s as nifti file, hdr size = %d, should be 348 or 540\n', filename, hdr.sizeof_hdr);
else
  % the file is now open with the appropriate little or big-endianness
end

if hdr.sizeof_hdr==348
  % the remainder of the code is for nifti-2 files
  error('%s seems to be a nifti-1 file', filename)
end

hdr.magic           = fread(fid, [1 8 ], 'int8=>int8'     ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
hdr.datatype        = fread(fid, [1 1 ], 'int16=>int16'   ); % 12      See file formats
hdr.bitpix          = fread(fid, [1 1 ], 'int16=>int16'   ); % 14      See file formats
hdr.dim             = fread(fid, [1 8 ], 'int64=>double'  ); % 16      See file formats

if hdr.dim(1)<1 || hdr.dim(1)>7
  % see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim.html
  error('inconsistent endianness in the header');
end

hdr.intent_p1       = fread(fid, [1 1 ], 'double=>double' ); % 80      0
hdr.intent_p2       = fread(fid, [1 1 ], 'double=>double' ); % 88      0
hdr.intent_p3       = fread(fid, [1 1 ], 'double=>double' ); % 96      0
hdr.pixdim          = fread(fid, [1 8 ], 'double=>double' ); % 104     0,1,1,1,1,1,1,1
hdr.vox_offset      = fread(fid, [1 1 ], 'int64=>int64'   ); % 168     Offset of data, minimum=544
hdr.scl_slope       = fread(fid, [1 1 ], 'double=>double' ); % 176     1
hdr.scl_inter       = fread(fid, [1 1 ], 'double=>double' ); % 184     0
hdr.cal_max         = fread(fid, [1 1 ], 'double=>double' ); % 192     0
hdr.cal_min         = fread(fid, [1 1 ], 'double=>double' ); % 200     0
hdr.slice_duration  = fread(fid, [1 1 ], 'double=>double' ); % 208     0
hdr.toffset         = fread(fid, [1 1 ], 'double=>double' ); % 216     0
hdr.slice_start     = fread(fid, [1 1 ], 'int64=>int64'   ); % 224     0
hdr.slice_end       = fread(fid, [1 1 ], 'int64=>int64'   ); % 232     0
hdr.descrip         = fread(fid, [1 80], 'int8=>char'     ); % 240     All zeros
hdr.aux_file        = fread(fid, [1 24], 'int8=>char'     ); % 320     All zeros
hdr.qform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 344     NIFTI_XFORM_UNKNOWN (0)
hdr.sform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 348     NIFTI_XFORM_UNKNOWN (0)
hdr.quatern_b       = fread(fid, [1 1 ], 'double=>double' ); % 352     0
hdr.quatern_c       = fread(fid, [1 1 ], 'double=>double' ); % 360     0
hdr.quatern_d       = fread(fid, [1 1 ], 'double=>double' ); % 368     0
hdr.qoffset_x       = fread(fid, [1 1 ], 'double=>double' ); % 376     0
hdr.qoffset_y       = fread(fid, [1 1 ], 'double=>double' ); % 384     0
hdr.qOffset_z       = fread(fid, [1 1 ], 'double=>double' ); % 392     0
hdr.srow_x          = fread(fid, [1 4 ], 'double=>double' ); % 400     0,0,0,0
hdr.srow_y          = fread(fid, [1 4 ], 'double=>double' ); % 432     0,0,0,0
hdr.srow_z          = fread(fid, [1 4 ], 'double=>double' ); % 464     0,0,0,0
hdr.slice_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 496     0
hdr.xyzt_units      = fread(fid, [1 1 ], 'int32=>int32'   ); % 500     0xC (seconds, millimeters)
hdr.intent_code     = fread(fid, [1 1 ], 'int32=>int32'   ); % 504     See file formats
hdr.intent_name     = fread(fid, [1 16], 'int8=>char'     ); % 508     See file formats
hdr.dim_info        = fread(fid, [1 1 ], 'int8=>int8'     ); % 524     0
hdr.unused_str      = fread(fid, [1 15], 'int8=>char'     ); % 525     All zeros
% disp(ftell(fid));                                          % 540     End of the header

fclose(fid);

end

function write_nifti2_hdr(filename, hdr)

% WRITE_NIFTI2_HDR
%
% Use as
%   write_nifti2_hdr(filename, hdr)
% where
%   filename   = string
%   hdr        = structure with nifti-2 header information
%
% See also READ_NIFTI_HDR, READ_CIFTI, WRITE_CIFTI

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if ischar(filename)
  fid = fopen(filename, 'w');
else
  fid = filename;
end

fwrite(fid, 540, 'int32');

assert(fwrite(fid, hdr.magic          , 'int8'   )==8 ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
assert(fwrite(fid, hdr.datatype       , 'int16'  )==1 ); % 12      See file formats
assert(fwrite(fid, hdr.bitpix         , 'int16'  )==1 ); % 14      See file formats
assert(fwrite(fid, hdr.dim            , 'int64'  )==8 ); % 16      See file formats
assert(fwrite(fid, hdr.intent_p1      , 'double' )==1 ); % 80      0
assert(fwrite(fid, hdr.intent_p2      , 'double' )==1 ); % 88      0
assert(fwrite(fid, hdr.intent_p3      , 'double' )==1 ); % 96      0
assert(fwrite(fid, hdr.pixdim         , 'double' )==8 ); % 104     0,1,1,1,1,1,1,1
assert(fwrite(fid, hdr.vox_offset     , 'int64'  )==1 ); % 168     Offset of data, minimum=544
assert(fwrite(fid, hdr.scl_slope      , 'double' )==1 ); % 176     1
assert(fwrite(fid, hdr.scl_inter      , 'double' )==1 ); % 184     0
assert(fwrite(fid, hdr.cal_max        , 'double' )==1 ); % 192     0
assert(fwrite(fid, hdr.cal_min        , 'double' )==1 ); % 200     0
assert(fwrite(fid, hdr.slice_duration , 'double' )==1 ); % 208     0
assert(fwrite(fid, hdr.toffset        , 'double' )==1 ); % 216     0
assert(fwrite(fid, hdr.slice_start    , 'int64'  )==1 ); % 224     0
assert(fwrite(fid, hdr.slice_end      , 'int64'  )==1 ); % 232     0
assert(fwrite(fid, hdr.descrip        , 'int8'   )==80); % 240     All zeros
assert(fwrite(fid, hdr.aux_file       , 'int8'   )==24); % 320     All zeros
assert(fwrite(fid, hdr.qform_code     , 'int32'  )==1 ); % 344     NIFTI_XFORM_UNKNOWN (0)
assert(fwrite(fid, hdr.sform_code     , 'int32'  )==1 ); % 348     NIFTI_XFORM_UNKNOWN (0)
assert(fwrite(fid, hdr.quatern_b      , 'double' )==1 ); % 352     0
assert(fwrite(fid, hdr.quatern_c      , 'double' )==1 ); % 360     0
assert(fwrite(fid, hdr.quatern_d      , 'double' )==1 ); % 368     0
assert(fwrite(fid, hdr.qoffset_x      , 'double' )==1 ); % 376     0
assert(fwrite(fid, hdr.qoffset_y      , 'double' )==1 ); % 384     0
assert(fwrite(fid, hdr.qOffset_z      , 'double' )==1 ); % 392     0
assert(fwrite(fid, hdr.srow_x         , 'double' )==4 ); % 400     0,0,0,0
assert(fwrite(fid, hdr.srow_y         , 'double' )==4 ); % 432     0,0,0,0
assert(fwrite(fid, hdr.srow_z         , 'double' )==4 ); % 464     0,0,0,0
assert(fwrite(fid, hdr.slice_code     , 'int32'  )==1 ); % 496     0
assert(fwrite(fid, hdr.xyzt_units     , 'int32'  )==1 ); % 500     0xC (seconds, millimeters)
assert(fwrite(fid, hdr.intent_code    , 'int32'  )==1 ); % 504     See file formats
assert(fwrite(fid, hdr.intent_name    , 'int8'   )==16); % 508     See file formats
assert(fwrite(fid, hdr.dim_info       , 'int8'   )==1 ); % 524     0
assert(fwrite(fid, hdr.unused_str     , 'int8'   )==15); % 525     All zeros
% disp(ftell(fid));                                     % 540     End of the header

if ischar(filename)
  fclose(fid);
end

end