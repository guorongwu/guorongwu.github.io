%% demo code for voxel-wise HRF deconvolution
%% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
%% Guo-Rong Wu, gronwu@gmail.com, UESTC, UGent, 2013.9.12
%% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D.. 
%% A blind deconvolution approach to recover effective connectivity brain networks 
%% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .
clc,clear
%% Mask file
brainmask = 'E:\spm8\toolbox\DPARSF_V2.1_120101\Templates\BrainMask_05_61x73x61.hdr';%change the mask as you like
brain = spm_read_vols(spm_vol(brainmask));
data_tmp = zeros(size(brain));
voxel_ind = find(brain>0); %% change as you like.
num_voxel = length(voxel_ind);
nobs = 240; % number of time points
bsig = zeros(nobs,num_voxel); 

%% open Matlab parallel computing, NumWorkers: set a reasonable number yourself. If you don't have parallel facilities no prob, but change "parfor" to normal "for"
try 
    myCluster = parcluster('local');
    myCluster.NumWorkers = 8; 
    saveAsProfile(myCluster,'local2'); 
    matlabpool open 'local2' 8
end
TR = 2; %in seconds 
thr = 1; % threshold, for example 1 SD.
event_lag_max = 5; % the (estimated) maximum lagged time from neural event to BOLD event, in points. 


main= 'G:\RS_data\';  % data directory
data_dir = fullfile(main,'FunImgNormalizedCovremovedDetrendedFiltered'); %% where data are stored after smoothing,regression, filtering, detrending or whatever preprocessing
sub = dir(data_dir); 
sub(1:2)=[]; %this removes the "." and ".." 
save_dir = fullfile(main,'FunImgNor_HRF'); %% save dir, change the name as you like.
save_dir2 = fullfile(main,'FunImgNor_Deconv'); %% save dir, change the name as you like.

% make directories for the HRF parameters
mkdir(fullfile(save_dir,'Height'));
mkdir(fullfile(save_dir,'T2P'));
mkdir(fullfile(save_dir,'FWHM'));
mkdir(save_dir2);

for isub=1:length(sub)
    disp('Reading data ...')
    sub_dir = fullfile(data_dir,sub(isub).name);
    disp(sub(isub).name)
    cd(sub_dir);
    clear imag
    % if your preprocessed data are not stored in an image, but in a
    % vector, you can call this vector rsig and skip the following lines
    imag = dir('*.img'); %% if *.nii, change it yourself.
    tic
    rsig = zeros(size(imag,1),num_voxel);
    parfor k = 1:length(imag)
        [data1] = spm_read_vols(spm_vol(imag(k).name));
        rsig(k,:) =  data1(voxel_ind);
    end
    toc
    disp('Done') 
    
%     tic
%     rsig = spm_detrend(rsig,3); % make sure stability
%     toc
%     disp('Finishing detrending')   
    disp('Retrieving HRF ...'); 
    tic
        [data_deconv onset hrf event_lag PARA] = wgr_deconv_canonhrf_par(rsig,thr,event_lag_max,TR)
    toc
    disp('Done');   
    
    save(fullfile(save_dir,[sub(isub).name,'_hrf.mat']),'event_lag_max','thr','TR','onset', 'hrf','event_lag', 'PARA','-v7.3');
    
    
    % Write HRF parameter
    % Height - h
    % Time to peak - p (in time units of TR seconds)
    % FWHM (at half peak) - w  

    v=spm_vol(brainmask);
    v.dt=[16,0]; 
    
    v.fname = fullfile(save_dir,'Height',[sub(isub).name,'_height.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(1,:);
    spm_write_vol(v,data);

    v.fname = fullfile(save_dir,'T2P',[sub(isub).name,'_Time2peak.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(2,:);
    spm_write_vol(v,data);

    v.fname = fullfile(save_dir,'FWHM',[sub(isub).name,'_FWHM.nii']);
    data = data_tmp;
    data(voxel_ind)=PARA(3,:);
    spm_write_vol(v,data);
    
    
    sub_save_dir = fullfile(save_dir2,sub(isub).name);
    mkdir(sub_save_dir)
    % writting back into nifti files
    for k = 1:length(imag)
        v.fname = fullfile(sub_save_dir,imag(k).name);
        data = data_tmp;
        data(voxel_ind) = data_deconv(k,:);
        spm_write_vol(v,data);
    end   
    
end
