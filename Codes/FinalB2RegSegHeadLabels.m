%rana.banik@vanderbilt.edu
close all;
% clear all;
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819/')

rootDir = '/media/banikr2/DATA/emModels/MatlabTestRun/';
ls = dir(rootDir);
ls = ls(3:end);
ls.name; % takes the subfolders right under rootDir
for nSub = 1:length(ls) %number of subjects
    subID = ls(nSub).name;  
    display(['Subject ' num2str(nSub) ': ' subID]);
    cd([rootDir,subID]);
    %% reorient the native MR and CT images:
    cmd = ['fslreorient2std '  subID '_CT.nii.gz '  subID '_CT2std.nii.gz'];
    system(cmd);
    cmd = ['fslreorient2std '  subID '_MR.nii.gz '  subID '_MR2std.nii.gz'];
    system(cmd);
    %% MAS file is not also in std orientation
    cmd =['fslreorient2std ' 'Multi_Atlas/orig_target_seg.nii.gz ' subID '_MAS_MR_seg2std.nii.gz'];
    system(cmd);
    %% load the CT and MR
    nii = load_untouch_nii_gz([subID '_CT2std.nii.gz']); 
    ct = nii.img; 
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
    mr = nii.img;
    %% register the native CT to native MR
    cmd =['flirt'...
    ' -in ' subID '_CT2std.nii.gz' ... 
    ' -ref ' subID '_MR2std.nii.gz'... %ref or fixed image: MR
    ' -out ' subID '_CTonMR.nii.gz'... %output native CT registered on native MR
    ' -omat ' subID '_CTonMR.mat'... % 4x4 affine transformation matrix(ascii format)
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'... %search angle -90 90
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'... %3 rotation, 3 translation, 3 shear, 3 scale
    ' -interp trilinear'];
    system(cmd);
    
    nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
    ct = nii.img; %CT file
    %% extract brain mask from CT
    cmd = ['/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/./ct_bet ' subID '_CTonMR.nii.gz']; %creates the bet mask 
    system(cmd);
    nii = load_untouch_nii_gz([subID '_CTonMR_brain_mask.nii.gz']);
    brainMask = nii.img;
    %% FAST segmentation on MR brain
    mr(brainMask==0) = 0; 
    mr(mr<0)=0.5;
    nii.img = mr;
    save_untouch_nii_gz(nii, [subID '_MR_brain.nii.gz']); %orientation ok with MR2std
    %%
    cmd = ['flirt'... % to get the standard2input matrix
    ' -in ' rootDir 'average305_t1_tal_lin.nii'...
    ' -ref ' subID '_MR2std.nii.gz'... 
    ' -out ' subID '_MNI2MR.nii.gz'... 
    ' -omat ' subID '_MNI2MR.mat '...
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
    system(cmd);
    %% MR to MNI 
    cmd = ['flirt'... % to get the standard2input matrix
    ' -in ' subID '_MR2std.nii.gz'...
    ' -ref ' subID '_MR_norm.nii.gz'... % rootDir 'average305_t1_tal_lin.nii'... 
    ' -out ' subID '_MR2norm.nii.gz'... 
    ' -omat ' subID '_MR2norm.mat '...
    ' -bins 256'...
    ' -cost normmi'...
    ' -searchrx -90 90'...
    ' -searchry -90 90'...
    ' -searchrz -90 90'...
    ' -dof 12'...
    ' -interp trilinear'];
    system(cmd);
    %% apply the xformation matrix to MR_brain
    % this part is not necessary I suppose coz we apply standard2input
    % xformation matrix on MR_brain through FAST
    cmd = ['flirt'... 
    ' -in ' subID '_MR_brain.nii.gz'... 
    ' -ref ' subID '_MR_norm.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
    ' -out ' subID '_MR_brain_norm.nii.gz'... 
    ' -applyxfm -init ' subID '_MR2norm.mat'...%'_MR_TO_MNI.mat'... 
    ' -interp nearestneighbour'];
    system(cmd);
    %% FAST
    % with probability mapping
    cmd = ['fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o -a ' subID '_MNI2MR.mat ' subID '_MR_brain '  subID '_MR_brain'];
    system(cmd);
    % fsl2reorient here ??? Not done
    nii = load_untouch_nii_gz(['-a_pveseg.nii.gz']);  %subID '_MR_brain
    fst_lab = nii.img;  % >> the output of fsl fast 
    nii = load_untouch_nii_gz([subID '_MAS_MR_seg2std.nii.gz']); 
    atl_lab = nii.img;
    %% Merge labels
    [gx,gy,gz] = imgradientxyz(ct,'prewitt'); 
    ct_grad = sqrt(gx.^2+gy.^2+gz.^2); 
    nii.img = ct_grad; 
%     save_untouch_nii_gz(nii, [subID '_CT_gradmag.nii.gz']); % no fsl2re..
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']); 
    mr = nii.img;    
    [gx,gy,gz] = imgradientxyz(mr,'prewitt'); 
    mr_grad = sqrt(gx.^2+gy.^2+gz.^2); 
    nii.img = mr_grad; 
%     save_untouch_nii_gz(nii, [subID '_MR_gradmag.nii.gz']);  % no fsl2re
    
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
            lab(s.PixelIdxList) = 3; 
        end
    end
    % label thalmus and pallidas
    lab(atl_lab==55|atl_lab==56) = 56;
    lab(atl_lab==59|atl_lab==60) = 60;
    % label electrodes : ct>2000
    stat = regionprops(ct>2000&mr<10, 'Area', 'PixelIdxList'); % ct -> CT2MR
    for nn = 1:length(stat)
        s = stat(nn);
        if s.Area<10
            continue;
        end % remove small area 
        [~,ii] = max(ct(s.PixelIdxList)); % find coords of max ct voxel                #todo: shall we run it on normed? 
        [xx,yy,zz] = ind2sub(sz, s.PixelIdxList(ii));
        [~, bi] = regionGrowing(ct, [xx, yy, zz], 1600, 50); % region grow (1600 could be changed based on whether the result bi include skull)
        % Threshold = 1600; maxDistance = 50 or inf
        % the seed point comes from max intensity which should not be norm
        % specific
        bi_ring = imdilate(bi, strel('sphere', 1)) - bi; % ACCRE: disk->ball
        lab(bi==1|(bi_ring==1&ct>1600)|(bi_ring==1&ct_grad>30000)) = 20;
    end % looks good
    % label bones : ct=[300,2000]
    bi = (lab==0&ct>370&ct<2000&mr<150&ct_grad<30000); % 300-1500 % value 400 gets rid of peripheral ct components of machine
    stat = regionprops(bi, 'Area', 'PixelIdxList');
    for nn=1:length(stat)
        s = stat(nn);
        if s.Area<1000 
            continue; 
        end % remove small area
        bi = zeros(sz); bi(s.PixelIdxList) = 1;
        bi = imfill(bi,'holes'); % fill holes
        %             for sl = (sz(3)/2):sz(3) % check every axial slice to remove electrode
        %                 ss = regionprops(bi(:,:,sl), 'Area', 'PixelIdxList');
        %                 if length(ss)==1, continue; end;
        %
        %             end
%         bi = imopen(bi, strel('sphere', 1)); 
%         bi_ring = bi - imerode(bi, strel('sphere',1)); 
%         bi(bi_ring>0&mr>50) = 0;
%         bi = imopen(bi, strel('sphere', 1));
%         bi_ring = bi - imerode(bi, strel('sphere',1)); 
%         bi(bi_ring>0&mr>50) = 0;
%         bi = imclose(bi, strel('sphere', 1)); 
%         bi = imclose(bi, strel('sphere', 1));
        bi = imdilate(bi, strel('sphere', 1));
        lab(bi==1) = 4;
    end
    % label fat: ct=[-150,-20]
    bi = (lab==0&((ct>-150&ct<-20&mr>70)|mr>50&mr_grad>1500|mr>80)); % fat
    stat = regionprops(bi, 'Area', 'PixelIdxList');
    for nn=1:length(stat)
        s = stat(nn);
        if s.Area<500 
            continue; 
        end % remove small objects
        bi = zeros(sz); 
        bi(s.PixelIdxList) = 1;
        bi_ring = imdilate(bi, strel('sphere', 1)) - bi; 
        bi(bi_ring==1&mr>40) = 1; % fill gaps
        lab(bi==1) = 8;
    end
    % label muscle & skin: mr=[30,100]
    lab(lab==0&(mr>30&mr<100)) = 10;
    bi = lab==0&(mr>30&mr<100);
    % label water in eyeball
    
    %% Eyeball
    cmd = ['flirt'... 
    ' -in ' subID '_MR_brain.nii.gz'... 
    ' -ref ' subID '_MR_norm.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
    ' -out ' subID '_MR_brain_norm.nii.gz'... 
    ' -applyxfm -init ' subID '_MR2norm.mat'... 
    ' -interp nearestneighbour'];
    system(cmd);
    
    
    nii = load_untouch_nii_gz([subID '_MR_brain_norm.nii.gz']); 
    mr_brain_norm = nii.img;
    sz_norm = size(mr_brain_norm);
    [~, bi_r] = regionGrowing(mr_brain_norm, [124, 180, 36], 10, 15); %have to provide mr brain? 
    [~, bi_l] = regionGrowing(mr_brain_norm, [51, 180, 36], 10, 15);
    lab_norm = zeros(sz_norm);
    
    lab_norm(bi_r==1|bi_l==1) = 1;
    nii.img = lab_norm;
    save_untouch_nii_gz(nii, [subID '_eyeball_norm_lab.nii.gz']);
    % xform the label from norm to MR
    cmd = ['flirt'... 
    ' -in ' subID '_eyeball_norm_lab.nii.gz'... 
    ' -ref ' subID '_MR2std.nii.gz '... % if the ref is _MR_MNI.nii.gz output resolution turns bad. 
    ' -out ' subID '_eyeball_lab.nii.gz'... 
    ' -applyxfm -init ' subID '_MNI2MR.mat'... 
    ' -interp nearestneighbour'];
    system(cmd);
    cmd = ['fslreorient2std '  subID '_eyeball_lab.nii.gz '  subID '_eyeball_lab.nii.gz']; 
    system(cmd);
    
    nii = load_untouch_nii_gz([subID '_eyeball_lab.nii.gz']); 
    eyeball_lab = nii.img;
    
    lab(eyeball_lab==1) = 1;
    
    nii.img = lab;
    save_untouch_nii_gz(nii, [subID '_final_label_reduced_eyeball.nii.gz']);
%     break;
end