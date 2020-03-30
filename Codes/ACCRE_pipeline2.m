close all;
clear all;
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes', '/media/banikr2/DATA/emModels/test', ...
    '/media/banikr2/DATA/emModels/TwentySubjectAll/5767');
rootDir = '/media/banikr2/DATA/emModels/test/';
ls = dir(rootDir); ls = ls(3:end);
ls.name;
%%
for nSub = 1:length(ls) %number of subjects
    subID = ls(nSub).name;  
    display(['Subject ' num2str(nSub) ': ' subID]);
    cd([rootDir subID]);  
end

%%
nii = load_untouch_nii_gz([subID '_CT.nii.gz']); ct = nii.img; %CT file
% nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_seg.nii.gz']); msk = nii.img; %  

%%
cmd =['flirt'...
' -in /media/banikr/DATA/emModels/test/271_MR_bet_norm.nii.gz' ... %input or moving image: CT
' -ref /media/banikr/DATA/emModels/test/271/271_MR.nii.gz'... %ref or fixed image: MR
' -out /media/banikr/DATA/emModels/test/271_MR_bet2.nii.gz'... %output native CT registered on native MR
' -omat /media/banikr/DATA/emModels/test/271_MR_bet2.mat'... % 4x4 affine transformation matrix(ascii format)
' -bins 256'...
' -cost normmi'...
' -searchrx -90 90'... %search angle -90 90
' -searchry -90 90'...
' -searchrz -90 90'...
' -dof 12'... %3 rotation, 3 translation, 3 shear, 3 scale
' -interp trilinear'];
system(cmd)
%%
Data = load('CT-MRI-273_MR_TO_MR_norm_0.mat','-ASCII');
%%
nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']); ct = nii.img; %CT file  %%%%%%%%CT_on_MR
nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_seg.nii.gz']); msk = nii.img; %     
ring = imdilate(msk==1, strel('sphere', 2))-imerode(msk==1, strel('sphere', 2)); 
ring(ring>0&ct>200) = 2; msk(ring>0) = ring(ring>0); 
stat = regionprops(msk==1, 'Area', 'PixelIdxList'); 
for ii=1:length(stat)
    if stat(ii).Area>10000; continue; end;
    msk(stat(ii).PixelIdxList) = 2;
end
nii.img = msk;
save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_seg1.nii.gz']);
% Manual correct mask on bottom spinal part
% Brain extraction
nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_seg2.nii.gz']); mask = nii.img==1; 
nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm.nii.gz']); mr = nii.img; 
mr(mask==0) = 0; mr(mr<0)=0.5;
nii.img = mr;
save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_brain.nii.gz']);   
%% bet on norm: bet doesn't work if there is .nii file of the same .gz
cmd = ['bet /media/banikr/DATA/emModels/TwentySubjectAll/5767/CT-MRI-5767_MR_norm.nii.gz /media/banikr/DATA/emModels/test/CT-MRI-5767_MR_norm_robust_smallf_ng.nii.gz -m -n -f 0.15 -R -S -t'];
system(cmd)
% if head segmentation required CT -> best result. 
% -f 0.2259 on CT both norm and native gets the entire head. 
% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET/UserGuide
% https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro2/index.html#bet_optional
% bet doesn't perform as the CT-MRI-273_MR_norm_seg
% nii = load_untouch_nii_gz('/media/banikr/DATA/emModels/test/CT-MRI-5767_CT_norm_robust_smallf_ng_mask.nii.gz'); 
% volumeViewer(nii.img)
%% bet 
cmd = ['bet /media/banikr/DATA/emModels/TwentySubjectAll/5767/CT-MRI-5767_MR_norm_seg.nii.gz /media/banikr/DATA/emModels/TwentySubjectAll/5767_seg_mask_bet.nii.gz -m -f 0.15 -R'];
system(cmd)
%% Checking fslreorient2std
cmd = ['fslreorient2std /media/banikr/DATA/emModels/test/271_MR_bet.nii.gz /media/banikr/DATA/emModels/test/271_MR_bet.nii.gz']; %have to write twice for some strange syntax
system(cmd)
cmd = ['fslreorient2std /media/banikr/DATA/emModels/test/271/271_MR.nii.gz /media/banikr/DATA/emModels/test/271/271_MR.nii.gz'];
system(cmd)
%%
% Correct brain mask
% subID = '5767'
% nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']); ct = nii.img; % CT file  %CT_on_MR
% nii = load_untouch_nii_gz(['CT-MRI-',subID,'_MR_norm_seg.nii.gz']); msk = nii.img; %     
% ring = imdilate(msk==1, strel('sphere', 2))-imerode(msk==1, strel('sphere', 2)); 
% ring(ring>0&ct>200) = 2; msk(ring>0) = ring(ring>0); 
% stat = regionprops(msk==1, 'Area', 'PixelIdxList'); 
% for ii=1:length(stat)
%     if stat(ii).Area>10000; continue;
%     end;
%     msk(stat(ii).PixelIdxList) = 2;
% end
% nii.img = msk;
% save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_seg1.nii.gz']);
%% trial on bet
subID = '5767'
nii = load_untouch_nii_gz([subID '_MR_norm_bet_mask.nii.gz']); 
mrBetMask = nii.img; 
% 271_norm_CT.nii.gz
nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']); 
ct = nii.img;
dilatedBrain = imdilate(mrBetMask, strel('sphere', 2));
erodedBrain = imerode(mrBetMask, strel('sphere', 2));
brainResidue = dilatedBrain - erodedBrain; % >> is the ring

brainResidue(brainResidue>0&ct>200) = 2; % >> ct > 200 is bones/skull and ring>0 is only 1; 2 means skull? 
mrBetMask(brainResidue>0) = brainResidue(brainResidue>0); % similar as before!

stat = regionprops(mrBetMask==1, 'Area', 'PixelIdxList'); 
for ii=1:length(stat)
    if stat(ii).Area>10000; 
        continue; 
    end;
    mrBetMask(stat(ii).PixelIdxList) = 2; % >> enriching skull/bone values
end
% nii.img = mask012; % >>detailed skull mask012 
%% betsurf

%% bet on binary mask
cmd = ['bet /media/banikr/DATA/emModels/test/CT-MRI-5767_brainCT.nii.gz /media/banikr/DATA/emModels/test/CT-MRI-5767_brainCT.nii.gz -m -f 0.15'];
system(cmd)
%% creating seg label
subID = '5767';
nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']);
ct = nii.img;
nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm.nii.gz']);
mr = nii.img;
mr(mr<0)=0.5;
sz = size(ct);
seg012 = zeros(sz);
seg012(ct>0 & ct<81)=1; % CSF, WM, GM, Acute blood
seg012(mr == 0) = 0;
seg012(ct>190) = 2;
msk = imdilate(seg012==1, strel('sphere',2));
% volumeViewer(seg012)
% for i = 1:156
% seg012 = regionfill(seg012==1);
% end
volumeViewer(seg012)
% seg012 = imdilate(seg012==1, strel('sphere',2));
% volumeViewer(seg012)
% nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_seg2.nii.gz']); mask = nii.img;%==1;  %% Do not have this segmentation also
% volumeViewer(mask)
%% all about masks
subID = '5767';
nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']);
ct = nii.img;
sz = size(ct);
% csf + wm + gm
brainCT = zeros(sz);
brainCT(ct>-1&ct<81)=1;
volumeViewer(brainCT)
% only csf 
csf = zeros(sz);
csf(ct>-1&ct<11)=1;
volumeViewer(csf)
csf(ct>-1&ct<11)=1;
nii = load_untouch_nii_gz(['CT-MRI-5767_CT_norm_robust_smallf_ng_mask.nii.gz']);
skullmask = nii.img;
ct(skullmask==0) = 0;
% skullCTmask(skullmask==1) = ct;
volumeViewer(ct)
nii.img = ct;
% save_untouch_nii_gz(nii, ['CT-MRI-5767_skull_CT.nii.gz']);
nii = load_untouch_nii_gz(['CT-MRI-5767_CT_norm_robust_smallf_ng_mask.nii.gz']);
skmaskCT = nii.img;
%% ct_bet_master
nii = load_untouch_nii_gz(['CT-MRI-5767_CT_norm_brain_mask.nii.gz']);
betCTmsk =  nii.img;

% check performance with fsl_bet
nii = load_untouch_nii_gz(['CT-MRI-5767_MR_norm_brain.nii.gz']);
mrBrain = nii.img;

nii = load_untouch_nii_gz(['CT-MRI-5767_MR_norm.nii.gz']);
mr = nii.img;

% brainCropBet = betCTmsk.*mr;
brainCropBet = zeros(sz);
brainCropBet(betCTmsk==1) = mr(betCTmsk==1); % it cropped out the brian from the mr
nii.img=brainCropBet;
save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_brainCTbet.nii.gz']);
%% Check FSL fast on ct_bet
% label 3 - white matter
% label 2 - grey matter
% label 1 - csf
cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o  CT-MRI-' subID '_MR_norm_brainCTbet CT-MRI-' subID '_MR_norm_brainCTbet'];
system(cmd);
%% after FAST
lab = zeros(sz);
nii = load_untouch_nii_gz(['CT-MRI-5767_MR_norm_brainCTbet_seg.nii.gz']);
fst_lab = nii.img; 
volumeViewer(fst_lab)
lab(fst_lab==1) = 1; 
lab(fst_lab==2) = 2; 
lab(fst_lab==3) = 3;
stat = regionprops(fst_lab==2, 'Area', 'PixelIdxList'); % props of grey matter
for nn = 1:length(stat) % the number of grey areas found
    s = stat(nn);
    if s.Area>=2  % ig they are greater than 2 voxels then surely grey 
        continue; % do nothing, next iteration
    end
    if fst_lab(s.PixelIdxList+1)==2 % smaller areas less than 2 if white
        lab(s.PixelIdxList) = 3; % convert it in grey in new label
    end
end
volumeViewer(lab)
%% Check if the previous MAS and current MAS has same/similar labels
% how many labels in the image: 
finalLabel = niftiread('/media/banikr/DATA/emModels/MRCT_result/5639/CT-MRI-5639_MR_norm_lab_2_.nii.gz');
numel(unique(finalLabel)) % 8 
batch2MASLabel = niftiread('/media/banikr/DATA/emModels/test/271/Multi_Atlas/orig_target_seg.nii.gz');
numel(unique(batch2MASLabel)) % 133
batch1MASLabel = niftiread('/media/banikr/DATA/emModels/TwentySubjectAll/5767/Multi_Atlas/SEG/orig_target_seg.nii.gz');
numel(unique(batch1MASLabel)) % 133
% check thalamus pallidus
isThalPall = zeros(size(batch1MASLabel));
isThalPall(batch1MASLabel==55|batch1MASLabel==56|batch1MASLabel==59|batch1MASLabel==60) = 1;
volumeViewer(isThalPall)
% use separate labels for thalamus and pallidas. . . 
%% Regiongrowing
subID = '5767';
nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']);
ct = nii.img;
nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm.nii.gz']);
mr = nii.img;
lab = zeros(sz);
stat = regionprops(ct>2000&mr<10, 'Area', 'PixelIdxList');
for nn = 1:length(stat)
    s = stat(nn);
    if s.Area<10
        continue;
    end % remove small area      
    [~,ii] = max(ct(s.PixelIdxList)); % find coords of max ct voxel
    [xx,yy,zz] = ind2sub(sz, s.PixelIdxList(ii));

    [~, bi] = regionGrowing(ct, [xx, yy, zz], 1600, 50); % region grow (1600 could be changed based on whether the result bi include skull)
    bi_ring = imdilate(bi, strel('sphere', 1)) - bi; % ACCRE: disk->ball
    lab(bi==1|(bi_ring==1&ct>1000)) = 20; % |(bi_ring==1&ct_grad>30000)
end
volumeViewer(lab)
%% The main code here

% load native MR
% load native CT
% 1. fslorient2std : 

% 2. flirt CT on MR : CTonMR, xform.mat
% 3. ct_bet_python on CTonMR : brainMask (CSF, WM, GM) + also can include the
%                            robust fslbet in the pipeline. 


% 4. Mask out brain from MR native
% 5. FSL FAST
%% check the pipeline
close all;
clear all;
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819')
rootDir = '/media/banikr2/DATA/emModels/MatlabTestRun/';
ls = dir(rootDir);
ls = ls(3:end);
ls.name; % takes the subfolders right under rootDir
for nSub = 1:length(ls) %number of subjects
    subID = ls(nSub).name;  
    display(['Subject ' num2str(nSub) ': ' subID]);
    cd([rootDir,subID]);
    
    cmd = ['fslreorient2std '  subID '_CT.nii.gz '  subID '_CT.nii.gz']; %have to write twice for some strange syntax
    system(cmd)
    cmd = ['fslreorient2std '  subID '_MR.nii.gz '  subID '_MR.nii.gz'];
    system(cmd)
%     nii = load_untouch_nii_gz([subID '_CT.nii.gz']); 
%     ct = nii.img; %CT file  %native CT has no variable use I guess
%     nii = load_untouch_nii_gz([subID '_MR.nii.gz']);
%     mr = nii.img;
%     cmd =['flirt'...
%     ' -in ' subID '_CT.nii.gz' ... %input or moving image: CT
%     ' -ref ' subID '_MR.nii.gz'... %ref or fixed image: MR
%     ' -out ' subID '_CTonMR.nii.gz'... %output native CT registered on native MR
%     ' -omat ' subID '_CTonMR.mat'... % 4x4 affine transformation matrix(ascii format)
%     ' -bins 256'...
%     ' -cost normmi'...
%     ' -searchrx -90 90'... %search angle -90 90
%     ' -searchry -90 90'...
%     ' -searchrz -90 90'...
%     ' -dof 12'... %3 rotation, 3 translation, 3 shear, 3 scale
%     ' -interp trilinear'];
%     system(cmd)
% %     
%     nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
%     ct = nii.img; %CT file
%     cmd = ['./ct_bet ' rootDir subID '/' subID '_CTonMR.nii.gz']; %creates the bet mask 
%     system(cmd)
%     
%     nii = load_untouch_nii_gz([rootDir subID '/' subID '_CTonMR_brain_mask.nii.gz']);
%     mask = nii.img; 
%     % shall the mask be finetuned by morphological operations? 
%     mr(mask==0) = 0; 
%     mr(mr<0)= 0.5;
%     nii.img = mr;
%     save_untouch_nii_gz(nii, [subID '_MR_brain.nii.gz']); % save the brain
%     
%     % have to perform FSL FAST in norm space of brain. We need MR to norm
%     % xformation matrix for that. 
    
    cmd = ['flirt'...
    ' -in ' subID '_MR.nii.gz'... 
    ' -ref ' subID '_MR_norm.nii.gz '...
    ' -out ' subID '_MR_TO_MR_norm.nii.gz'... 
    ' -omat ' subID '_MR_TO_MR_norm.mat '...
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
    system(cmd);
    
%     cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o -a ' subID '_MR_TO_MR_norm.mat ' subID '_MR_brain '  subID '_MR_brain'];
    cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o  ' subID '_MR_brain '  subID '_MR_brain'];
    % syntax: fast without .nii.gz and twice.
    % generated seg file. 
    system(cmd);
    break;
end
%%
% cmd = ['./ct_bet ' rootDir subID '/' subID '_CTonMR.nii.gz'];
% system(cmd)
cmd = ['flirt'...
    ' -in ' subID '_MR.nii.gz'... 
    ' -ref ' rootDir 'average305_t1_tal_lin.nii'...%subID '_MR_norm.nii.gz '...
    ' -out ' subID '_MR_2_MNI.nii.gz'... 
    ' -omat ' subID '_MR_2_MNI.mat '...
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
    system(cmd);
%% region growing on MR norm head
nii = load_untouch_nii_gz('309_MR_norm.nii.gz');
mr_norm = nii.img;
[~, bi_r] = regionGrowing(mr_norm, [124, 176, 25], 10, 20); %have to provide mr brain? 
[~, bi_l] = regionGrowing(mr_norm, [55, 175, 25], 10, 20);
eye2MRnorm = zeros(size(mr_norm));
eye2MRnorm(bi_r==1|bi_l==1) = 1;
% volumeViewer(mr_norm)
nii.img = eye2MRnorm;
save_untouch_nii_gz(nii, ['309_eyeball_norm_lab.nii.gz']);
%% accelerated eyeball
