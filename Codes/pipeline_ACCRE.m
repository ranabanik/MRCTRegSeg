% % MR-CT pair

% Data?

% MR space
% CT-MRI-####_MR.nii.gz                     (MR: ~170x256x256, 1x1x1mm)
% /Multi_Atlas/SEG/orig_target_seg.nii.gz   (Lab:~170x256x256, 1x1x1mm)

% CT space
% CT-MRI-####_CT.nii.gz                     (MR: 512x512x~272, 0.47x0.48x0.62mm)

% Common space (MNI 305)
% CT-MRI-####_CT_norm.nii.gz                (CT: 172x220x156, 1x1x1mm)
% CT-MRI-####_MR_norm.nii.gz                (MR: 172x220x156, 1x1x1mm, old name: average305_t1_tal_lin_via_CT-MRI-5652.nii.gz)
% CT-MRI-####_MR_norm_seg.nii.gz            (Lab:172x220x156,1x1x1mm, old name: average305_t1_tal_lin_via_CT-MRI-5652_seg.nii.gz 
%                                                "1"-brain; "2"-CSF )

% Normal CT Brain
% Tissue            HU
% Ca+               > 1000
% Acute blood       60 ~ 80
% Gray matter       32 ~ 42
% White matter      22 ~ 32
% CSF               0  ~ 10
% Fat              -50 ~ -80
% Air               -1000

close all; 
clear all;
addpath('/home/gaoy3/matlab/NIFTI_20100819/', '/scratch/gaoy3/CT-MR/'); % scratch is in ACCRE
rootDir = '/scratch/gaoy3/CT-MR/data/';

% list all subjects
ls = dir(rootDir); ls = ls(3:end); 

% cmd = 'export FSLDIR=/scratch/gaoy3/fsl'; system(cmd); 
% cmd = '. /scratch/gaoy3/fsl/etc/fslconf/fsl.sh'; system(cmd); 
% cmd = 'export FSL_DIR=$FSLDIR'; system(cmd); 
% cmd = 'export PATH=$PATH:/scratch/gaoy3/fsl:/scratch/gaoy3/fsl/bin'; system(cmd); 

for nSub = 1:length(ls)
    subID = ls(nSub).name;  display(['Subject ' num2str(nSub) ': ' subID]);
    
    cd([rootDir subID]); 
    
    % 0. FSL FLIRT (MR/CT -> MNI 305) 
    
    
    % 1. Lab WM/GM/CSF on MR (FSL FAST)
    % Correct brain mask
    nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']); ct = nii.img;
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
    % FSL FAST
    cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o  CT-MRI-' subID '_MR_norm_brain CT-MRI-' subID '_MR_norm_brain'];
    system(cmd);
   
    % 2. Lab subcortical on MR (multi atlas)
    % FSL FLIRT: MR.nii.gz -> MR_norm.nii.gz (MR_norm_seg.nii.gz)
    cmd = ['flirt -in CT-MRI-' subID '_MR.nii.gz -ref CT-MRI-' subID '_MR_norm.nii.gz '...
        '-out CT-MRI-' subID '_MR_TO_MR_norm.nii.gz -omat CT-MRI-' subID '_MR_TO_MR_norm.mat '...
        '-bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear'];
    system(cmd);
    cmd = ['flirt -in Multi_Atlas/SEG/orig_target_seg_improv.nii.gz -ref CT-MRI-' subID '_MR_norm.nii.gz '...
        '-out CT-MRI-' subID '_MR_TO_MR_norm_shadowreg_multiatlas_seg.nii.gz '...
        '-applyxfm -init CT-MRI-' subID '_MR_TO_MR_norm.mat -interp nearestneighbour'];
    system(cmd);

    % 3. Merge labs
    nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_TO_MR_norm_shadowreg_multiatlas_seg.nii.gz']); atl_lab = nii.img;
    nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_brain_seg.nii.gz']); fst_lab = nii.img;    
    nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm.nii.gz']); ct = nii.img; 
    [gx,gy,gz] = imgradientxyz(ct,'prewitt'); ct_grad = sqrt(gx.^2+gy.^2+gz.^2); 
    nii.img = ct_grad; save_untouch_nii_gz(nii, ['CT-MRI-' subID '_CT_norm_gradmag.nii.gz']); % ACCRE only
    nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm.nii.gz']); mr = nii.img;    
    [gx,gy,gz] = imgradientxyz(mr,'prewitt'); mr_grad = sqrt(gx.^2+gy.^2+gz.^2); 
    nii.img = mr_grad; save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_gradmag.nii.gz']); % ACCRE only
%     nii = load_untouch_nii_gz(['CT-MRI-' subID '_CT_norm_gradmag.nii.gz']); ct_grad = nii.img;
%     nii = load_untouch_nii_gz(['CT-MRI-' subID '_MR_norm_gradmag.nii.gz']); mr_grad = nii.img;
    sz = size(mr);

    lab = zeros(sz);
    
    % label CSF / GM / WM
    lab(fst_lab==1) = 1; lab(fst_lab==2) = 2; lab(fst_lab==3) = 3;
    stat = regionprops(fst_lab==2, 'Area', 'PixelIdxList');
    for nn = 1:length(stat)
        s = stat(nn);
        if s.Area>=2, continue; end;
        if fst_lab(s.PixelIdxList+1)==3, lab(s.PixelIdxList) = 3; end;
    end
    
    % label thalmus and pallidas
    lab(atl_lab==55|atl_lab==56|atl_lab==59|atl_lab==60) = 2; % ACCRE
    
    % label electrodes : ct>2000
    stat = regionprops(ct>2000&mr<10, 'Area', 'PixelIdxList');
    for nn = 1:length(stat)
        s = stat(nn);
        if s.Area<10, continue; end; % remove small area
        
        [~,ii] = max(ct(s.PixelIdxList)); % find coords of max ct voxel
        [xx,yy,zz] = ind2sub(sz, s.PixelIdxList(ii));
        
        [~, bi] = regionGrowing(ct, [xx, yy, zz], 1600, 50); % region grow (1600 could be changed based on whether the result bi include skull)
        bi_ring = imdilate(bi, strel('sphere', 1)) - bi; % ACCRE: disk->ball
        lab(bi==1|(bi_ring==1&ct>1000)|(bi_ring==1&ct_grad>30000)) = 20;
    end
    
    % label bone : ct=[300,2000]
    bi = (lab==0&ct>300&ct<2000&mr<150&ct_grad<30000); % 300-1500
    stat = regionprops(bi, 'Area', 'PixelIdxList');
    for nn=1:length(stat)
        s = stat(nn);
        if s.Area<1000, continue; end; % remove small area
        bi = zeros(sz); bi(s.PixelIdxList) = 1;
        bi = imfill(bi,'holes'); % fill holes
        %             for sl = (sz(3)/2):sz(3) % check everay axial slice to remove electrode
        %                 ss = regionprops(bi(:,:,sl), 'Area', 'PixelIdxList');
        %                 if length(ss)==1, continue; end;
        %
        %             end
        bi = imopen(bi, strel('sphere', 1)); 
        bi_ring = bi - imerode(bi, strel('sphere',1)); bi(bi_ring>0&mr>50) = 0;
        bi = imopen(bi, strel('sphere', 1));
        bi_ring = bi - imerode(bi, strel('sphere',1)); bi(bi_ring>0&mr>50) = 0;
        bi = imclose(bi, strel('sphere', 1)); bi = imclose(bi, strel('sphere', 1));
        lab(bi==1) = 4;
    end
    
    % label fat: ct=[-150,-20]
    bi = (lab==0&((ct>-150&ct<-20&mr>70)|mr>50&mr_grad>1500|mr>80)); % fat
    stat = regionprops(bi, 'Area', 'PixelIdxList');
    for nn=1:length(stat)
        s = stat(nn);
        if s.Area<500, continue; end; % remove small objects
        bi = zeros(sz); bi(s.PixelIdxList) = 1;
        bi_ring = imdilate(bi, strel('sphere', 1)) - bi; bi(bi_ring==1&mr>40) = 1; % fill gaps
        lab(bi==1) = 8;
    end
    
    % label muscle & skin: mr=[30,100]
    lab(lab==0&(mr>30&mr<100)) = 10;
    
    % label water in eyeball
    [~, bi_r] = regionGrowing(mr, [124, 180, 36], 10, 20);
    [~, bi_l] = regionGrowing(mr, [51, 180, 36], 10, 20);
    lab(bi_r==1|bi_l==1) = 1;
    
    nii.img = lab;
    save_untouch_nii_gz(nii, ['CT-MRI-' subID '_MR_norm_lab.nii.gz']);

%     ct_lab = zeros(size(ct)); 
%     ct_lab(ct>300) = 4; % bone
%     bi = ct_lab ==4; bi = imclose(bi, strel('sphere', 3)); ct_lab(bi>0) = 4;
%     ct_lab(ct>-20   & ct<=300) = 6; % white/gray matter/csf/musle
%     ct_lab(ct>-150 & ct<=-20) = 8; % fat
%     ct_lab(ct<=-150) = 0; % air   

end

cd('/scratch/gaoy3/CT-MR/');