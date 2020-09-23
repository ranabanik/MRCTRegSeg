rootDir = '/media/banikr2/DATA/emModels/CT-MR_batch2/SeptemberModels2/';
ls = dir(rootDir);
ls = ls(3:end); % to remove . and .. folders
ls.name;
for nSub = 1:length(ls)
    subID = ls(nSub).name
    cd([rootDir,subID]);
    cmd = ['fslreorient2std '  subID '_CT.nii.gz '  subID '_CT2std.nii.gz'];
    system(cmd);
    cmd = ['fslreorient2std '  subID '_MR.nii.gz '  subID '_MR2std.nii.gz'];
    system(cmd);
    
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
    
    cd(rootDir);
    cmd = ['fslmaths MNI305_eyeMask.nii.gz -mul 0.9 eye_weight.nii.gz'];
    system(cmd)
    cmd = ['fslmaths average305_t1_tal_lin.nii.gz -bin -mul 0.1 -add eye_weight.nii.gz refweight'];
    system(cmd)
    
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
end
%% 
mainDir = '/media/banikr2/DATA/emModels/CT-MR_batch2/CTonMRDone/';
addpath('/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/',...
    '/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Codes/NIFTI_20100819/');
ls = dir(mainDir);
ls = ls(3:end-1);
rangelist = zeros(4,length(ls));
for nSub = 1:length(ls)
    subID = ls(nSub).name; 
    cd([mainDir,subID]);
    display(['Subject ' num2str(nSub) ': ' subID])
    nii = load_untouch_nii_gz([subID '_CTonMR.nii.gz']); 
    ct = nii.img; 
    nii = load_untouch_nii_gz([subID '_MR2std.nii.gz']);
    mr = nii.img;
    rangelist(3,nSub) = max(mr(:));
    rangelist(4,nSub) = min(mr(:));
    fatCTLabel = zeros(size(ct));
    fatCTLabel(ct>-150&ct<-20&mr~=0) = 1; %[-150,-30]
    mr(fatCTLabel~=1) = 0;
    % rangelist(3,nSub) = max(mr(:));
    % rangelist(4,nSub) = min(mr(:));
    J = mr(mr(:)~=0);
%     size(J)
    H = histogram(J);
    % break;
    % end
    % histfit(J, 10)
    GMModel = fitgmdist(double(J),2);
    rangelist(1,nSub) = max(GMModel.mu); % max(mr(:));
%     rangelist(2,nSub) = GMModel.mu(2); % min(mr(:));
    rangelist(2,nSub) = str2double(subID);
% xgrid = linspace(-4,12,1001)';
% hold on; plot(xgrid,pdf(GMModel,xgrid),'r-'); hold off
end
plot(rangelist(1,:),'r')
%% Soft window tissue normalization
L = 40; % window
W = 200; % level or range
maxThr = L+W;
minThr = l-W;
ct = (255*(ct - minThr)) / (maxThr-minThr);
%% About GMModel
% J = mr(mr(:)~=0);
% histogram(J)
% sum(fst_lab(:)==2)
x = [randn(4000,1)/2; 5+2*randn(6000,1)];
% plot(x)
histogram(x,'Normalization','pdf')
%% plot
plot(rangelist(1,:),'r');%legend1('fat gaussian mean')
hold on
plot(rangelist(3,:),'b');%legend2('max mr intensity')
hold on
plot(rangelist(4,:),'g');%legend3('min mr intensity')
xticklabels('fontsize',20)