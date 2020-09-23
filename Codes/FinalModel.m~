close all;
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819/',...
    '/media/banikr2/DATA/emModels/MatlabTestRun/')
rootDir = '/media/banikr2/DATA/emModels/CT-MR_batch2/SeptemberModels/'; %CTonMRDone/';
ls = dir(rootDir);
ls = ls(3:end);
ls.name;
subID = '307';
cd([rootDir,subID]);
% MAS file is not also in std orientation
cmd =['fslreorient2std ' 'Multi_Atlas/orig_target_seg.nii.gz ' subID '_MAS_MR_seg2std.nii.gz'];
system(cmd);
nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
mr = nii.img;
nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
ct = nii.img;
% extract brain mask from CT 
cmd = ['/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/./ct_bet ' subID '_CTonMR.nii.gz']; %creates the bet mask 
system(cmd);
nii = load_untouch_nii_gz([subID '_CTonMR_brain_mask.nii.gz']);
brainMask = nii.img;
%    registration from MR standard to MNI
% September: is it needed? 
% cmd = ['flirt'...
%     ' -in ' subID '_MR2std.nii.gz'... 
%     ' -ref ' rootDir 'average305_t1_tal_lin.nii'...%subID '_MR_norm.nii.gz '...
%     ' -out ' subID '_MRstd_2_MNI.nii.gz'...  % needed for eyeball seeding in MNI space
%     ' -omat ' subID '_MR_2_MNI.mat '...
%     ' -bins 256'...
%     ' -cost normmi'...
%     ' -searchrx -90 90'...
%     ' -searchry -90 90'...
%     ' -searchrz -90 90'...
%     ' -dof 12'...
%     ' -interp trilinear'];
% system(cmd);
mr(brainMask==0) = 0; 
nii.img = mr;
save_untouch_nii_gz(nii, [subID '_MR_brain.nii.gz']);
cmd = ['flirt'... % to get the standard2input matrix
    ' -in ' rootDir 'average305_t1_tal_lin.nii'...
    ' -ref ' subID '_MR2std.nii.gz'... 
    ' -out ' subID '_MNI2MR.nii.gz'... 
    ' -omat ' subID '_MNI2MR.mat '...  % needed in FAST
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
system(cmd);
cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o -a ' subID '_MNI2MR.mat ' subID '_MR_brain '  subID '_MR_brain'];
system(cmd);
nii = load_untouch_nii_gz('-a_pveseg.nii.gz');  %subID '_MR_brain
fst_lab = nii.img;  % >> the output of fsl fast 
nii = load_untouch_nii_gz([subID '_MAS_MR_seg2std.nii.gz']); 
atl_lab = nii.img;
[gx,gy,gz] = imgradientxyz(ct,'prewitt'); 
ct_grad = sqrt(gx.^2+gy.^2+gz.^2); 
nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']); 
mr = nii.img;    
[gx,gy,gz] = imgradientxyz(mr,'prewitt'); 
mr_grad = sqrt(gx.^2+gy.^2+gz.^2); 
sz = size(mr);
lab = zeros(sz);
lab(fst_lab==1) = 1; lab(fst_lab==2) = 2; lab(fst_lab==3) = 3;
stat = regionprops(fst_lab==2, 'Area', 'PixelIdxList');
for nn = 1:length(stat)
    s = stat(nn);
    if s.Area>=2 % 
        continue;  % do nothing, next iteration
    end
    if fst_lab(s.PixelIdxList+1)==3
        lab(s.PixelIdxList) = 3;  % increases the number of WM pixels and reduces GM pixels
    end
end
lab(atl_lab==55|atl_lab==56) = 56;
lab(atl_lab==59|atl_lab==60) = 56; 
%% Electrodes are cancelled
stat = regionprops(ct>2000&mr<10, 'Area', 'PixelIdxList'); % 1600; ct -> CT2MR mr<10 % 6000 for high intensity CT volumes
for nn = 1:length(stat)
    s = stat(nn);
    if s.Area<10
        continue;
    end % remove small area 
    [~,ii] = max(ct(s.PixelIdxList)); % find coords of max ct voxel                #todo: shall we run it on normed? 
    [xx,yy,zz] = ind2sub(sz, s.PixelIdxList(ii));
    if zz > size(lab,3)/2 % avoid taking teeths
        [~, bi] = regionGrowing(ct, [xx, yy, zz], 1600, 10); % region grow (1600 could be changed based on whether the result bi include skull)
    % Threshold = 1600; maxDistance = 50 or inf
    % the seed point comes from max intensity which should not be norm
    % specific
        bi_ring = imdilate(bi, strel('sphere', 2)) - bi; % ACCRE: disk->ball
%         lab(bi==1|(bi_ring==1&ct>1400)|(bi_ring==1&ct_grad>30000)) = 20;
        lab(lab==0& bi==1|(bi_ring==1&ct>2000)|(bi_ring==1&ct_grad>30000)) = 20;
    end
end
%% Bone/Skull 
bi = (lab==0&ct>400&ct<2000&ct_grad<30000&mr~=0); % ct<2000; 300-1500 % value 400 gets rid of peripheral ct components of machine
bi_superior = bi(:,:,156:end);
bi_backQ = bi_superior(:,1:128,:); %posterior half of superior volume
bi_backQ = imfill(bi_backQ,'holes');
bi_superior(:,1:128,:) = bi_backQ;
R = 1.5; % radius of structuring element
[x,y,z] = meshgrid(-R:R); % X Y Z 3d matrices
el = x.^2+y.^2+z.^2 <= R^2; % sphere of radius 'R'
bi_eroded = imerode(bi_superior, el);
bi_dilated = imdilate(bi_eroded, el); 
%     bi(:,:,166:end) = bi_dilated(:,:,11:end); % skipping the lower 10 slices as they are eroded too much
cc = bwconncomp(bi_dilated);
numPixels = cellfun(@numel, cc.PixelIdxList);
[~,idx] = max(numPixels);
bi_chunk = zeros(size(bi_dilated));
bi_chunk(cc.PixelIdxList{idx})=1;
bi(:,:,156:end) = bi_chunk;
stat = regionprops(bi, 'Area', 'PixelIdxList');
for nn=1:length(stat)
    s = stat(nn);
    if s.Area<1000 
        continue; 
    end % remove small area
    bi = zeros(sz); 
    bi(s.PixelIdxList) = 1;
    bi = imfill(bi,'holes'); 
    bi = imdilate(bi, strel('sphere', 1));
    lab(bi==1) = 4; %&lab==0
end
%% Fat
fatCTLabel = zeros(size(ct));
fatCTLabel(ct>-150&ct<-20&mr~=0) = 1; % [-150,-30]
mr(fatCTLabel~=1) = 0;
J = mr(mr(:)~=0);
Y = prctile(J,[2.5, 97.5], 'all');
fatMRLabel = zeros(size(mr));
fatMRLabel(mr>min(Y) & mr<max(Y))=1;
mr(fatMRLabel~=1) = 0;
mrFatFlat = mr(:);
J = mrFatFlat(mrFatFlat~=0);
GMModel = fitgmdist(double(J),2);
sFatRng = max(GMModel.mu)-120; bFatRng = max(GMModel.mu)+120; % no explanation for 100
nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
mr = nii.img;
fatMRLabel = zeros(size(mr));
fatMRLabel(mr>sFatRng&mr<bFatRng) = 1;
fat_ring = imdilate(fatMRLabel, strel('sphere', 3)) - fatMRLabel; 
fatMRLabel(fat_ring==1&(mr>sFatRng&mr<bFatRng)) = 1;
lab(fatMRLabel==1&lab==0) = 8;
% fatMRLabel = fatMRLabel==1&lab==0;
%% Eyeballs
% check September update in checkutility.m
% cmd = ['flirt'...
%     ' -in ' subID '_MRstd_2_MNI.nii.gz'... 
%     ' -ref ' rootDir 'average305_t1_tal_lin.nii.gz'...
%     ' -refweight ' rootDir 'refweight.nii.gz'... 
%     ' -out ' subID '_MRstd2MNI4eye.nii.gz'];
% system(cmd)
% nii2 = load_untouch_nii_gz([subID, '_Eye_MNI.nii.gz']); % different sizes of norm so nii2 instead of nii
nii2 = load_untouch_nii_gz([subID, '_MR2std.nii.gz']);
eye_norm = nii2.img;
[~, bi_r] = regionGrowing(eye_norm, [60, 194, 133], 10, 20); % [121, 178, 30], 10, 20  % have to provide mr brain? [124,177,29] 
[~, bi_l] = regionGrowing(eye_norm, [133, 195, 133], 10, 20); %[55, 178, 30], 10, 20);
eye2MRnorm = zeros(size(eye_norm));
eye2MRnorm(bi_r==1|bi_l==1) = 1;
v = zeros(size(eye_norm));
v(eye2MRnorm==1) = eye_norm(eye2MRnorm==1);
v_flat = v(v~=0);
Y = prctile(v_flat,[2.5, 97.5], 'all');
eyeball = zeros(size(mr));
eyeball(mr>min(Y)&mr<max(Y))=1;
eyeball_copy = eyeball;
eyeball_copy(:,:,1:ceil(sz(3)/3)) = 0;
eyeball_copy(:,1:ceil(sz(2)/2)+5,:) = 0;
eyeball_copy(lab~=0) = 0;
erodedEyeball = imerode(eyeball_copy, strel('sphere', 1));
cc = bwconncomp(erodedEyeball);
numPixels = cellfun(@numel, cc.PixelIdxList);
finalEye = zeros(size(mr));
[~,idx] = max(numPixels); % idx = index of maximum volume of number of voxels
finalEye(cc.PixelIdxList{idx})= 1;
erodedEyeball(cc.PixelIdxList{idx})= 0;
cc = bwconncomp(erodedEyeball);
numPixels = cellfun(@numel, cc.PixelIdxList);
[~,idx] = max(numPixels);
finalEye(cc.PixelIdxList{idx})= 1; % finalEye2 in cases where eyeball is missing
EyeLab=imdilate(finalEye, strel('sphere', 3)); % 2
lab(EyeLab==1)=11; % 07.09.2020 eyelabels should be separate. 
%% Muscle, Skin
% msLabCT = zeros(size(ct));
msLabCT = ((ct>-160&ct<40)&(lab==0&mr>0));
cc = bwconncomp(msLabCT);
numPixelClusters = cellfun(@numel, cc.PixelIdxList);
msLabFinal = zeros(size(mr));
for nn=1:cc.NumObjects
    if length(cc.PixelIdxList{nn})<10
        continue;
    end
    msLabFinal(cc.PixelIdxList{nn})=1;
end
msLabCT = imfill(msLabFinal,'holes'); % msLabCT before
lab(msLabCT==1) = 10;
%% binary MR mask
bi_mask = (mr>10);
z_mask = zeros(size(bi_mask));
stat = regionprops(bi_mask, 'Area', 'PixelIdxList');
for nn=1:length(stat)
    s = stat(nn);
    if s.Area<10 %1000 
        continue; 
    end % remove small objects
%     bi_mask = zeros(sz); 
    z_mask(s.PixelIdxList) = 1;
end
% volumeViewer(z_mask);
lab_jaw = lab(:,:,1:ceil(sz(3)/3));
bi_mask = (mr>10);
z_mask_jaw = z_mask(:,:,1:ceil(sz(3)/3));
lab_jaw(z_mask_jaw==1&lab_jaw==0)=10;
lab(:,:,1:ceil(sz(3)/3)) = lab_jaw; 
%% to convert electrodes into fat tissue
lab(lab==20)= 8;
%%
nii.img = lab;
save_untouch_nii_gz(nii, [subID '_label.nii.gz']);