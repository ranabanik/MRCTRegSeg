%% nearest neighbor - x
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819/',...
    '/media/banikr2/DATA/emModels/MatlabTestRun/')
rootDir = '/media/banikr2/DATA/emModels/CT-MR_batch2/CTonMRDone/';
cd([rootDir,subID]);
lab= lab_skin_done;
numel(unique(lab))
% zz : sagittal(R->L); xx : coronal(P->A); yy : axial(I->S) (170, 256, 256)
[czz, cxx, cyy] = ind2sub(sz, find(lab==10|lab==4)); % candidate points from either skin or bones
[gzz, gxx, gyy] = ind2sub(sz, find(lab==20)); % given poinst from electrodes
nn = compute_nearest_neighbour([czz,cxx,cyy], [gzz,gxx,gyy]); % nearest points
[nzz, nxx, nyy] = ind2sub(sz, nn);
% nn=zeros(size([gzz,gxx,gyy]));
% for i=1:numel(gzz)
%     nn(i) = ind2sub(sz, compute_nearest_neighbour([czz,cxx,cyy],
%     [gzz,gxx,gyy])); %doesn't work
% end
% lab([gzz,gxx,gyy]) = lab([nzz,nxx,nyy]);
% for i=1:numel(gxx) %worked
%     lab(gzz(i),gxx(i),gyy(i))=lab(nzz(i),nxx(i),nyy(i));
% end
nii.img = dem_lab;
save_untouch_nii_gz(nii, [subID '_just_electrodes.nii.gz']);
dem_lab = zeros(size(lab));
dem_lab([gzz,gxx,gyy])=1;
for i=1:numel(gxx)
    dem_lab(nzz(i),nxx(i),nyy(i))=1; %zz,xx,yy sequence perfectly matches with image. 
end
volumeViewer(dem_lab)
%%
lab = lab_superior;
elec = zeros(size(lab));
elec(lab==20) = 1;
I1 = lab;
% I1(lab~=0)=1; % im2bw, regions are ones background zeros
% [BB,nb] = bwlabeln((lab~=20&lab~=0)&~elec);    % label regions?
% [BB,nb] = bwlabeln((lab~=0)&~elec);            % region without labels and without electrodes?
[BB,nb] = bwlabeln(I1&~elec);                    % labels different regions with different numbers
[EE,ne] = bwlabeln(elec);                        % label all electrode objects. 
II = BB*0;
for i=1:ne
    E1 = EE == i;                                % select one electrode
    E2 = imdilate(E1,strel('sphere', 4)); 
    for j = 1:nb
        B1 = BB == j;                            % select one region
        % imshowpair(E2,E1|B1)                   % 3D volumes won't work
        % pause(1)
        tmp = B1 & E2;                           % compare electrode and region
        if any(tmp(:))                             % if region and electrode are close
            II = II+j*E1;
            break
        end        
    end
end
A = lab(lab~=20)+II;
%% from darova
lab = lab_superior;     % the data shared
elec = zeros(size(lab));
elec(lab==20)=1;        % select electrodes
BB = lab;
BB(BB==20) = 0;         % remove electrodes
nb = unique(BB);        % unique regions
nb(nb==0)=[];
nb = flip(nb);
[EE,ne] = bwlabeln(elec);  
II = BB*0;
for i = 1:ne
    E1 = EE == i;   
    E2 = imdilate(E1,strel('sphere', 1)); 
    for j = nb(:)'
        B1 = BB == j;                                       % select one region
        tmp = B1 & E2;                                      % compare electrode and region
        if any(tmp(:))                                      % if region and electrode are close
            II = II+j*E1;
            break
        end        
    end
end
A = II + BB;
lab_sup = zeros(size(mr));
lab_sup(:,:,1:101) = A;
nii.img = lab_sup;
save_untouch_nii_gz(nii, 'darova1.nii.gz');
%%
nii = load_untouch_nii_gz('CT-MRI-5634_MR_norm_lab_2_.nii.gz');
norm_lab_5634 = nii.img;
sz_norm = size(norm_lab_5634);
[gzz, gxx, gyy] = ind2sub(sz_norm, find(norm_lab_5634==20)); 
[czz, cxx, cyy] = ind2sub(sz_norm, find(norm_lab_5634~=20&norm_lab_5634~=0)); % candidate points from either skin or bones
% kills the process
nn = compute_nearest_neighbour([czz, cxx, cyy], [gzz,gxx,gyy]); % nearest points
[nzz, nxx, nyy] = ind2sub(sz_norm, nn);
dem_lab = zeros(size(norm_lab_5634));
for i=1:numel(gxx) %worked
    dem_lab(gzz(i),gxx(i),gyy(i))=1;
end
% dem_lab([gzz,gxx,gyy]) = 1;   %dem_lab([nzz,nxx,nyy]);
% max([gzz, gxx, gyy])
%%
msLabCT= ((ct>-160&ct<40)&(lab==0&mr>0));
% uses mr intensity which is not standardized
stat = regionprops(msLabCT, 'Area', 'PixelIdxList');
lab_skin = zeros(sz);
for nn=1:length(stat)
    s = stat(nn);
    if s.Area<100
        continue;
    end % remove small area
    bi = zeros(sz);
    bi(s.PixelIdxList) = 1;
    bi = imfill(bi,'holes'); % fill holes
    bi = imopen(bi, strel('sphere', 1));
    bi_ring = bi - imerode(bi, strel('sphere',1)); bi(bi_ring>0&mr>50) = 0; %mr intensity is not standardized
    bi = imopen(bi, strel('sphere', 1));
    bi_ring = bi - imerode(bi, strel('sphere',1)); bi(bi_ring>0&mr>50) = 0;
    bi = imclose(bi, strel('sphere', 1)); bi = imclose(bi, strel('sphere', 1));
    lab_skin(bi==1) = 10;
end
%%
% msLabCT((ct>-160&ct<40)&(lab==0&mr>0))=1;
msLabCT = imfill(msLabCT,'holes');
msLabCT_dilated2 = imdilate(msLabCT, strel('sphere', 2));
lab_b4_skin(msLabCT_dilated2==1&lab_b4_skin==0)=10;
% nii.img = lab_b4_skin;
% save_untouch_nii_gz(nii, [subID '_half_lab_good_skin.nii.gz']); % save the brain
msLabCT_dilated2(msLabCT_dilated2==1&lab~=0)=0; % take only the outer most skins. 
msLabCT_dilated3 = imfill(msLabCT_dilated2,'holes');
msLabCT_dilated4 = imdilate(msLabCT_dilated2, strel('sphere', 1));
%% creating an eye mask?
nii = load_untouch_nii_gz([rootDir, 'average305_t1_tal_lin.nii.gz']);
mni305 = nii.img;
[~, bi_r] = regionGrowing(mni305, [121, 180, 26], 10, 20);
[~, bi_l] = regionGrowing(mni305, [55, 180, 26], 10, 20);
eyeMask = zeros(size(mni305));
eyeMask(bi_r==1|bi_l==1) = 1;
eyeMask = imdilate(eyeMask, strel('sphere', 1));
volumeViewer(eyeMask)
nii2.img = eyeMask;
save_untouch_nii_gz(nii2, 'MNI305_eyeMask.nii.gz');
%% 
cmd = ['flirt'...
' -in ' subID '_mr_norm.nii.gz'... % input MR is in MNI preregistered.
' -ref ' rootDir 'average305_t1_tal_lin.nii'...
' -refweight ' rootDir 'weightedMNI305Mask2std.nii.gz'... % weighted ref mask for eyes
' -out ' subID '_MRstd_2_MNI_eye3.nii.gz'...  % needed for eyeball seeding in MNI space
' -omat ' subID '_MR_2_MNI_eyemask3.mat'...
' -bins 256'...
' -cost normmi'...
' -searchrx -90 90'...
' -searchry -90 90'...
' -searchrz -90 90'...
' -dof 12'...
' -interp trilinear'];
system(cmd);
%% fnirt 
%fnirt --ref=<some template> --in=<some image> --infwhm=8,4,2 --subsamp=4,2,1 --warpres=8,8,8
cmd = ['fnirt'...
    ' --ref=' rootDir 'average305_t1_tal_lin.nii'...
    ' --in=' subID '_MRstd_2_MNI.nii.gz'...
    ' --refmask=' rootDir 'MNI305_eyeMask.nii.gz'... %not inmask
    ' --aff=' subID '_MR_2_MNI_eyemask.mat'... %output of affine registration
    ' --iout=' subID '_MR_NL_MNI.nii.gz'
%     ' --subsamp=' 4,2,1,1 ... 
%     ' --applyrefmask=' 1,1,1,1
    ];
system(cmd)
%% from Taylor Hanayik suggestion
weightedMNImask = zeros(size(mr_norm));
nii = load_untouch_nii_gz('MNIHeadMask.nii.gz');
headmaskMNI = nii.img;
nii = load_untouch_nii_gz('MNI305_eyeMask.nii.gz');
eyeMaskMNI = nii.img;
weightedMNImask(eyeMaskMNI==1)= 0.9;
weightedMNImask(headmaskMNI==1&eyeMaskMNI~=1)= 0.2;
unique(weightedMNImask)
nii.img = weightedMNImask;
save_untouch_nii_gz(nii, 'weightedMNImask.nii.gz'); %orientation ok with MR2std
%%
cmd = ['fslreorient2std '  '273_MRstd_2_MNI.nii.gz '  '273_MRstd_2_MNI.nii.gz'];
system(cmd);
cmd = ['fslreorient2std '  'average305_t1_tal_lin.nii.gz '  'Okay.nii.gz'];
%%
nii = load_untouch_nii_gz('weightedMNI305Mask2std.nii.gz');
weightedMNImask = nii.img;
%% worked. 
cmd = ['fslmaths MNI305_eyeMask.nii.gz -mul 0.9 eye_weight.nii.gz'];
system(cmd)

cmd = ['flirt'...
    ' -in ' subID '_MRstd_2_MNI.nii.gz'... %_MR2std %use normalized MR
    ' -ref ' rootDir 'average305_t1_tal_lin.nii.gz'...
    ' -refweight ' rootDir 'refweight.nii.gz'...
    ' -out ' subID '_Eye_found.nii.gz']; %'_taylor_mr_from_stdspace.nii.gz'];
system(cmd)
%% adding skin or removing skin gaps
lab_superior = niftiread('darova.nii.gz');
h = strel('disk', 5);
bone1 = imdilate((lab_superior==4), h);     % bone1 is binary
boneExt = bone1&(lab_superior~=4);          % select only dilated parts of the bones
% boneExt = bone1 - (lab_superior==4);      % same as previous line
skin1 = boneExt&(lab_superior==10);
volumeViewer(skin1);
% volumeViewer(lab_superior==10)
volumeViewer(boneExt)
% for i = 1:size()
%% September Eyeball
% works for 288, 282 
subID = '302'
cmd = ['fslmaths MNI305_eyeMask.nii.gz -mul 0.9 eye_weight.nii.gz'];
system(cmd)
cmd = ['fslmaths average305_t1_tal_lin.nii.gz -bin -mul 0.1 -add eye_weight.nii.gz refweight'];
system(cmd)
% fsl oriented is not required in 282 
% cmd = ['fslreorient2std '  'average305_t1_tal_lin.nii.gz '  'Okay.nii.gz'];
cd([rootDir,subID]);
cmd = ['flirt'...
    ' -in ' subID '_MR2std.nii.gz'... %_MR2std %use normalized MR
    ' -ref ' rootDir 'average305_t1_tal_lin.nii.gz'...
    ' -refweight ' rootDir 'refweight.nii.gz'...
    ' -out ' subID '_Eye_MNI.nii.gz'...%]; %'_taylor_mr_from_stdspace.nii.gz'];
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
system(cmd)
%% Is it possible to make muscle and skin label make more uniform little fuzzy? 
% after segmenting fat in lab file:
msLabCT = zeros(size(ct));
% msLabCT((ct>-160&ct<40)&mr>0)=1; % lab == 0 excluded
msLabCT((ct>-160&ct<40)&(lab==0&mr>0))=1;
% msLabCT = imfill(msLabCT,'holes'); % morphological operation not now
v = zeros(size(ct));
v(msLabCT==1) = mr(msLabCT==1);
v_flat = v(v~=0);
Y = prctile(v_flat,[2.5, 97.5], 'all');
% msLabMR = zeros(size(mr));
msLabMR = (mr>min(Y)&mr<max(Y)); %=1;
% volumeViewer(msLabMR);
% stat = regionprops(msLabMR, 'Area', 'PixelIdxList');
% for nn=1:length(stat)
%     s = stat(nn);
%     if s.Area<50 
%         continue; 
%     end % remove small objects
%     bi = zeros(sz); bi(s.PixelIdxList) = 1;
%     bi_ring = imdilate(bi, strel('sphere', 1)) - bi; 
%     bi(bi_ring==1&mr>40) = 1; % fill gaps
%     lab(bi==1) = 8;
% end
cc = bwconncomp(msLabMR);
numPixelClusters = cellfun(@numel, cc.PixelIdxList);
msLabFinal = zeros(size(mr));
for nn=1:cc.NumObjects
    if length(cc.PixelIdxList{nn})<10
        continue;
    end
%     msLabFinal = zeros(size(mr));
    msLabFinal(cc.PixelIdxList{nn})=1;
%     ms_ring = imdilate(msLabFinal, strel('sphere', 1)) - msLabFinal;
%     msLabFinal(ms_ring==1) = 1;
%     lab(msLabFinal==1)=10;
end
% ms_ring = imdilate(msLabFinal, strel('sphere', 1)) - msLabFinal;