% save nii.gz

% a wrapper of save_nii
% edit by Yurui Gao

function save_nii_gz(nii, fname)

if strcmp(fname(end-6:end), '.nii.gz')
    fname = fname(1:end-7);
elseif strcmp(fname(end-3:end), '.nii')
    fname = fname(1:end-4);
end

save_nii(nii, [fname '.nii']);
gzip([fname '.nii']);
delete([fname '.nii']);


