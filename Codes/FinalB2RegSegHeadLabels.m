%rana.banik@vanderbilt.edu
close all;
% clear all;
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819/',...
    '/media/banikr2/DATA/emModels/MatlabTestRun/')
rootDir = '/media/banikr2/DATA/emModels/CT-MR_batch2/CTonMRDone/';
ls = dir(rootDir);
ls = ls(3:end);
ls.name; % takes the subfolders right under rootDir
% for nSub = 1:length(ls) %number of subjects
    subID = '271'; %ls(nSub).name;  
%     display(['Subject ' num2str(nSub) ': ' subID]);
    cd([rootDir,subID]);
    % reorient the native MR and CT images:
%     cmd = ['fslreorient2std '  subID '_CT.nii.gz '  subID '_CT2std.nii.gz'];
%     system(cmd);
    cmd = ['fslreorient2std '  subID '_MR.nii.gz '  subID '_MR2std.nii.gz']; %second one MR2std? was just MR.nii.gz before
    system(cmd);
%     MAS file is not also in std orientation 04.29.2020
    cmd =['fslreorient2std ' 'Multi_Atlas/orig_target_seg.nii.gz ' subID '_MAS_MR_seg2std.nii.gz'];
    system(cmd);
%     load the CT and MR 04.29.2020
%     nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
%     ct = nii.img; 
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
    mr = nii.img;
%     register the native CT to native MR
%     cmd =['flirt'...
%     ' -in ' subID '_CT2std.nii.gz' ... 
%     ' -ref ' subID '_MR2std.nii.gz'... %ref or fixed image: MR
%     ' -out ' subID '_CTonMR.nii.gz'... %output native CT registered on native MR
%     ' -omat ' subID '_CTonMR.mat'... % 4x4 affine transformation matrix(ascii format)
%     ' -bins 256'...
%     ' -cost normmi'...
%     ' -searchrx -90 90'... %search angle -90 90
%     ' -searchry -90 90'...
%     ' -searchrz -90 90'...
%     ' -dof 12'... %3 rotation, 3 translation, 3 shear, 3 scale
%     ' -interp trilinear'];
%     system(cmd);
    nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
    ct = nii.img; %CT file
    % extract brain mask from CT 04.29.2020
    cmd = ['/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/./ct_bet ' subID '_CTonMR.nii.gz']; %creates the bet mask 
    system(cmd);
    nii = load_untouch_nii_gz([subID '_CTonMR_brain_mask.nii.gz']);
    brainMask = nii.img;
%    registration from MR standard to MNI
% Not performed in september
    cmd = ['flirt'...
    ' -in ' subID '_MR2std.nii.gz'... 
    ' -ref ' rootDir 'average305_t1_tal_lin.nii'...%subID '_MR_norm.nii.gz '...
    ' -out ' subID '_MRstd_2_MNI.nii.gz'...  % needed for eyeball seeding in MNI space
    ' -omat ' subID '_MR_2_MNI.mat '...
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
    system(cmd);
    % FAST segmentation on MR brain 04.29.2020
    mr(brainMask==0) = 0; 
%     mr(mr<0)=0.5; % there is no mr < 0 in batch 2
    nii.img = mr;
    save_untouch_nii_gz(nii, [subID '_MR_brain.nii.gz']); %orientation ok with MR2std
    % 04.29.2020
    % registration from MNI to MR standard space
    %     ' -in ' '/media/banikr2/DATA/emModels/CT-MR_batch2/average305_t1_tal_lin.nii'...
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
    disp('done')
%     %% MR to MNI 
%     cmd = ['flirt'... % to get the standard2input matrix
%     ' -in ' subID '_MR.nii.gz'...
%     ' -ref ' 'average305_t1_tal_lin.nii'... %subID '_MR_norm.nii.gz'... % ... 
%     ' -out ' subID '_MR2std.nii.gz'... 
%     ' -omat ' subID '_MR2std.mat '...
%     ' -bins 256'...
%     ' -cost normmi'...
%     ' -searchrx -90 90'...
%     ' -searchry -90 90'...
%     ' -searchrz -90 90'...
%     ' -dof 12'...
%     ' -interp trilinear'];
%     system(cmd);
%     %% apply the xformation matrix to MR_brain
%     % this part is not necessary I suppose coz we apply standard2input
%     % xformation matrix on MR_brain through FAST
%     cmd = ['flirt'... 
%     ' -in ' subID '_MR_brain.nii.gz'... 
%     ' -ref ' subID '_MR_norm.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
%     ' -out ' subID '_MR_brain_norm.nii.gz'... 
%     ' -applyxfm -init ' subID '_MR2norm.mat'...%'_MR_TO_MNI.mat'... 
%     ' -interp nearestneighbour'];
%     system(cmd);
    % FAST 04.29.2020
    % with probability mapping
    cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o -a ' subID '_MNI2MR.mat ' subID '_MR_brain '  subID '_MR_brain'];
    system(cmd);
    % fsl2reorient here ??? Not done
    nii = load_untouch_nii_gz('-a_pveseg.nii.gz');  %subID '_MR_brain
    fst_lab = nii.img;  % >> the output of fsl fast 
    nii = load_untouch_nii_gz([subID '_MAS_MR_seg2std.nii.gz']); 
    atl_lab = nii.img;
    % Merge labels 04.29.2020
    [gx,gy,gz] = imgradientxyz(ct,'prewitt'); 
    ct_grad = sqrt(gx.^2+gy.^2+gz.^2); 
%     nii.img = ct_grad; 
%     save_untouch_nii_gz(nii, [subID '_CT_gradmag.nii.gz']); % no fsl2re..
%     nii = load_untouch_nii_gz([subID '_CT_gradmag.nii.gz']);
%     ct_grad = nii.img;
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']); 
    mr = nii.img;    
    [gx,gy,gz] = imgradientxyz(mr,'prewitt'); 
    mr_grad = sqrt(gx.^2+gy.^2+gz.^2); 
%     nii.img = mr_grad; 
%     save_untouch_nii_gz(nii, [subID '_MR_gradmag.nii.gz']);  % no fsl2re
    % 04.29.2020
    sz = size(mr);
    lab = zeros(sz);
    % label 3 - white matter
    % label 2 - grey matter
    % label 1 - csf
    lab(fst_lab==1) = 1; lab(fst_lab==2) = 2; lab(fst_lab==3) = 3;
    stat = regionprops(fst_lab==2, 'Area', 'PixelIdxList');
    for nn = 1:length(stat)
        s = stat(nn);
        if s.Area>=2
            continue;  % do nothing, next iteration
        end
        if fst_lab(s.PixelIdxList+1)==3
            lab(s.PixelIdxList) = 3;  % increases the number of WM pixels and reduces GM pixels
        end
    end
    % label thalamus and pallidus
    lab(atl_lab==55|atl_lab==56) = 56;
    lab(atl_lab==59|atl_lab==60) = 56; % 60; cheanging both thalamus and pallidus to the same label
    % label electrodes : ct>2000
    stat = regionprops(ct>1600&mr<30, 'Area', 'PixelIdxList'); % 1600; ct -> CT2MR mr<10
    for nn = 1:length(stat)
        s = stat(nn);
        if s.Area>10  %<10
            continue;
        end % remove small area 
        [~,ii] = max(ct(s.PixelIdxList)); % find coords of max ct voxel                #todo: shall we run it on normed? 
        [xx,yy,zz] = ind2sub(sz, s.PixelIdxList(ii));
        if zz > size(lab,3)/2 % avoid taking teeths
            [~, bi] = regionGrowing(ct, [xx, yy, zz], 1600, 50); % region grow (1600 could be changed based on whether the result bi include skull)
        % Threshold = 1600; maxDistance = 50 or inf
        % the seed point comes from max intensity which should not be norm
        % specific
            bi_ring = imdilate(bi, strel('sphere', 2)) - bi; % ACCRE: disk->ball
            lab(bi==1|(bi_ring==1&ct>1400)|(bi_ring==1&ct_grad>30000)) = 20;
        end
    end % looks good
    % label bones : ct=[300,2000]
    bi = (lab==0&ct>370&ct<2000&ct_grad<30000&mr~=0); % ct<2000; 300-1500 % value 400 gets rid of peripheral ct components of machine
    % mr > 0 added. 
    % this creates a binary mask; excluded mr from the thresholding 
    % remove fake electrodes/bones present in skull exosurface:
    % how about other subjects ?
    stat = regionprops(bi, 'Area', 'PixelIdxList');
    for nn=1:length(stat)
        s = stat(nn);
        if s.Area<1000 
            continue; 
        end % remove small area
        bi = zeros(sz); 
        bi(s.PixelIdxList) = 1;
        bi = imfill(bi,'holes'); % fill holes
        % if you take label = 0 pixels then the following loop is not...
        % ...necessary I suppose. 
%                         for sl = (sz(3)/2):sz(3) % check every axial slice to remove electrode
%                             ss = regionprops(bi(:,:,sl), 'Area', 'PixelIdxList');
%                             if length(ss)==1 
%                                 continue; 
%                             end
%                         end
    %         bi = imopen(bi, strel('sphere', 1)); 
    %         bi_ring = bi - imerode(bi, strel('sphere',1)); 
    %         bi(bi_ring>0&mr>50) = 0;
    %         bi = imopen(bi, strel('sphere', 1));
    %         bi_ring = bi - imerode(bi, strel('sphere',1)); 
    %         bi(bi_ring>0&mr>50) = 0;
    %         bi = imclose(bi, strel('sphere', 1)); 
    %         bi = imclose(bi, strel('sphere', 1));
        bi = imdilate(bi, strel('sphere', 1));
        lab(bi==1) = 4; %&lab==0
    end
    %     bonewtelectrode = zeros(sz);
    %     bonewtelectrode(bi==1 & lab==0) = 1;
    %     volumeViewer(lab)
    % GMMFat
    fatCTLabel = zeros(size(ct));
    fatCTLabel(ct>-150&ct<-20&mr~=0) = 1; % [-150,-30]
    mr(fatCTLabel~=1) = 0;
    J = mr(mr(:)~=0);
    % newly added on 30th april
    Y = prctile(J,[2.5, 97.5], 'all');   % take it? 
    %     H = histogram(J);
    fatMRLabel = zeros(size(mr));
    fatMRLabel(mr>min(Y) & mr<max(Y))=1;
    mr(fatMRLabel~=1) = 0;
    mrFatFlat = mr(:);
    J = mrFatFlat(mrFatFlat~=0);
    % 30th april 
    GMModel = fitgmdist(double(J),2);
    sFatRng = max(GMModel.mu)-100; bFatRng = max(GMModel.mu)+100; % no explanation for 100
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
    mr = nii.img;
    fatMRLabel = zeros(size(mr));
    fatMRLabel(mr>sFatRng&mr<bFatRng) = 1;
    % new 
    %     stat = regionprops(fatMRLabel, 'Area', 'PixelIdxList');
    fat_ring = imdilate(fatMRLabel, strel('sphere', 1)) - fatMRLabel; 
    fatMRLabel(fat_ring==1&(mr>sFatRng&mr<bFatRng)) = 1; % fill gaps
    %     lab(fatCTLabel==1) = 8;
    %     volumeViewer(lab)
    lab(fatMRLabel==1&lab==0) = 8;
    % 05.07.2020 changed eyeball registration
    cmd = ['flirt'...
    ' -in ' subID '_MR2std.nii.gz'... 
    ' -ref ' rootDir 'average305_t1_tal_lin.nii.gz'...
    ' -refweight ' rootDir 'refweight.nii.gz'... 
    ' -out ' subID '_MRstd2MNI4eye.nii.gz'];
    system(cmd)
    nii2 = load_untouch_nii_gz([subID, '_MRstd2MNI4eye.nii.gz']); % different sizes of norm so nii2 instead of nii
    eye_norm = nii2.img;
    % clear nii
    [~, bi_r] = regionGrowing(eye_norm, [121, 178, 30], 10, 20); %have to provide mr brain? [124,177,29] 
    [~, bi_l] = regionGrowing(eye_norm, [55, 178, 30], 10, 20);
    eye2MRnorm = zeros(size(eye_norm));
    eye2MRnorm(bi_r==1|bi_l==1) = 1;
    v = zeros(size(eye_norm));
    v(eye2MRnorm==1) = eye_norm(eye2MRnorm==1);
    v_flat = v(v~=0);
    Y = prctile(v_flat,[2.5, 97.5], 'all');
    eyeball = zeros(size(mr));
    eyeball(mr>min(Y)&mr<max(Y))=1;
    %     eyeball(lab~=0) = 0;
    eyeball_copy = eyeball;
    eyeball_copy(:,:,1:ceil(sz(3)/3)) = 0; % removing 1/3 of the voxels from inferior 
    %     volumeViewer(eyeball_copy)
    eyeball_copy(:,1:ceil(sz(2)/2)+5,:) = 0;
    %     volumeViewer(eyeball_copy)
    eyeball_copy(lab~=0) = 0;
    erodedEyeball = imerode(eyeball_copy, strel('sphere', 1));
    %     volumeViewer(erodedEyeball)
    cc = bwconncomp(erodedEyeball);
    %     S = regionprops(cc,'Centroid');
    numPixels = cellfun(@numel, cc.PixelIdxList);
    %     numel(numPixels)
    %     length(numPixels)
    finalEye = zeros(size(mr));
    [~,idx] = max(numPixels);
    finalEye(cc.PixelIdxList{idx})= 1;
    %     volumeViewer(finalEye)
    erodedEyeball(cc.PixelIdxList{idx})= 0;
    cc = bwconncomp(erodedEyeball);
    %     S = regionprops(cc,'Centroid');
    numPixels = cellfun(@numel, cc.PixelIdxList);
    [~,idx] = max(numPixels);
    finalEye(cc.PixelIdxList{idx})= 1;
    %     volumeViewer(finalEye)
    EyeLab=imdilate(finalEye, strel('sphere', 2));
    %     volumeViewer(EyeLab)
    lab(EyeLab==1)=11;
    % nii.img = lab;
    % save_untouch_nii_gz(nii,[subID, '_afterEyelab2.nii.gz']);
    % muscle and skin
    msLabCT = zeros(size(ct)); % or mr
    % msLabCT((ct>-160 & ct<200) & lab==0 & (mr>0&mr<sFatRng)) = 1;
    msLabCT((ct>-160&ct<40)&(lab==0&mr>0))=1;
    msLabCT = imfill(msLabCT,'holes');
    lab(msLabCT==1) = 10;
    % nii.img = lab;
    % save_untouch_nii_gz(nii, [subID '_final_label.nii.gz']);
    %
    % take the MRstd to MNI space, it doesn't change the intensity dist
    % this is already done.
%     cmd = ['flirt'...
%     ' -in ' subID '_MR2std.nii.gz'... 
%     ' -ref ' rootDir 'average305_t1_tal_lin.nii'...%subID '_MR_norm.nii.gz '...
%     ' -out ' subID '_MRstd_2_MNI.nii.gz'... 
%     ' -omat ' subID '_MR_2_MNI.mat '...
%     ' -bins 256'...
%     ' -cost normmi'...
%     ' -searchrx -90 90'...
%     ' -searchry -90 90'...
%     ' -searchrz -90 90'...
%     ' -dof 12'...
%     ' -interp trilinear'];
%     system(cmd);

%     %% label fat: ct=[-150,-20]
% not implemented because of mr_grad is not valid
%     bi = (lab==0 & ((ct>-150&ct<-20&mr>70)|mr>50&mr_grad>1500|mr>80)); % fat
%     stat = regionprops(bi, 'Area', 'PixelIdxList');
%     for nn=1:length(stat)
%         s = stat(nn);
%         if s.Area<500 
%             continue; 
%         end % remove small objects
%         bi = zeros(sz); 
%         bi(s.PixelIdxList) = 1;
%         bi_ring = imdilate(bi, strel('sphere', 1)) - bi; 
%         bi(bi_ring==1&mr>40) = 1; % fill gaps
%         lab(bi==1) = 8;
%     end
%     % label muscle & skin: mr=[30,100]
%     lab(lab==0&(mr>30&mr<100)) = 10;
%     bi = lab==0&(mr>30&mr<100);
    % label water in eyeball
% Eyeball
%     cmd = ['flirt'... 
%     ' -in ' subID '_MR_brain.nii.gz'... 
%     ' -ref ' subID '_MR_norm.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
%     ' -out ' subID '_MR_brain_norm.nii.gz'... 
%     ' -applyxfm -init ' subID '_MR2norm.mat'... 
%     ' -interp nearestneighbour'];
%     system(cmd);

%     nii = load_untouch_nii_gz([subID '_MR_brain_norm.nii.gz']); 
%     mr_brain_norm = nii.img;
%     sz_norm = size(mr_brain_norm);
%     [~, bi_r] = regionGrowing(mr_brain_norm, [124, 180, 36], 10, 15); %have to provide mr brain? 
%     [~, bi_l] = regionGrowing(mr_brain_norm, [51, 180, 36], 10, 15);
%     lab_norm = zeros(sz_norm);
%     
%     lab_norm(bi_r==1|bi_l==1) = 1;
%     nii.img = lab_norm;
%     save_untouch_nii_gz(nii, [subID '_eyeball_norm_lab.nii.gz']);
%     % xform the label from norm to MR
%     cmd = ['flirt'... 
%     ' -in ' subID '_eyeball_norm_lab.nii.gz'... 
%     ' -ref ' subID '_MR2std.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
%     ' -out ' subID '_eyeball_lab.nii.gz'... 
%     ' -applyxfm -init ' subID '_MNI2MR.mat'... 
%     ' -interp nearestneighbour'];
%     system(cmd);
%     cmd = ['fslreorient2std '  subID '_eyeball_lab.nii.gz '  subID '_eyeball_lab.nii.gz']; 
%     system(cmd);
%     nii = load_untouch_nii_gz([subID '_eyeball_lab.nii.gz']); 
%     eyeball_lab = nii.img;
%     lab(eyeball_lab==1) = 1;
%     nii.img = lab;
%     save_untouch_nii_gz(nii, [subID '_final_label_reduced_eyeball.nii.gz']);
%     break;
% end
% %% remove electrodes from bones.
% clc
% % need further finetuning in phase 2 
% lab_wo_electrode = zeros(sz(1), sz(2), sz(3)-156); %100
% ind = 1;
% for i = sz(3)-100:sz(3) %128:256
%     % ax_slice = bi(:,:,i);
%     cc = bwconncomp(bi(:,:,i)); % slice #203 has 5 objects.
%     disp(cc)
%     S = regionprops(cc, 'Area', 'PixelIdxList');%'Centroid');
%     ax_slice = zeros(sz(1), sz(2));
%     numPixels = cellfun(@numel, cc.PixelIdxList);
%     % disp(numPixels)
%     [biggest,idx] = max(numPixels);
%     % ax_slice(cc.PixelIdxList{idx})= 1;
%     lab_wo_electrode(:,:,ind) = ax_slice;
%     ind=ind+1;
% end
%% binary MR mask
bi_mask = (mr>10);
z_mask = zeros(size(bi_mask));
stat = regionprops(bi_mask, 'Area', 'PixelIdxList');
for nn=1:length(stat)
    s = stat(nn);
    if s.Area<100 %1000 
        continue; 
    end % remove small objects
%     bi_mask = zeros(sz); 
    z_mask(s.PixelIdxList) = 1;
end
% volumeViewer(z_mask);
%
lab_jaw = lab(:,:,1:ceil(sz(3)/3));
bi_mask = (mr>10);
z_mask_jaw = z_mask(:,:,1:ceil(sz(3)/3));
lab_jaw(z_mask_jaw==1&lab_jaw==0)=10;
lab(:,:,1:ceil(sz(3)/3)) = lab_jaw; 
nii.img = lab;
save_untouch_nii_gz(nii, [subID '_final_label_check.nii.gz']);