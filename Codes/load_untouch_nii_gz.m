% load nii.gz 
% edit by Yurui Gao

function nii = load_untouch_nii_gz(fname)

if strcmp(fname(end-6:end), '.nii.gz') %compares the sorted chars are '.nii.gz'
    fname = fname(1:end-7);  %'273_CT'
elseif strcmp(fname(end-3:end), '.nii')
    fname = fname(1:end-4);
end

gunzip([fname '.nii.gz']); 
nii = load_untouch_nii([fname '.nii']);

if exist([fname '.nii.gz'], 'file')
    delete([fname '.nii'])
end