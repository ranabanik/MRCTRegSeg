% load nii.gz 

% edit by Yurui Gao


function nii = load_nii_gz(fname)

if strcmp(fname(end-6:end), '.nii.gz')
    fname = fname(1:end-7);
elseif strcmp(fname(end-3:end), '.nii')
    fname = fname(1:end-4);
end

gunzip([fname '.nii.gz']); 
nii = load_nii([fname '.nii']);

if exist([fname '.nii.gz'], 'file')
    delete([fname '.nii'])
end