%% preprocessing pipeline for rs and movie data
% created by kaden on 10-24-2019 for the purposes of MSc thesis project

%% ========================= Create Directories ===========================
clear all
close all
clc

% create paths to preprocessing software locations
fsleyes = '/usr/local/fsl/bin/fsleyes';
fsl = '/usr/local/fsl/bin/';
afni = '/Users/kaden_shearer/abin/';

% create control subject directory
subj_path = '/Users/kaden_shearer/Dropbox/kaden_thesis/movie_data/varsity_data/controls/';
subj_dir = dir([subj_path,'n*']);

%% ========================= RS preprocessing ===========================

for ii = 1%:length(subj_dir)
    
    tmp_subj = subj_dir(ii).name;
    disp(['============== Working on ',tmp_subj,' RS data... =================='])
    rs_path = [subj_path,tmp_subj,'/rs/norm/'];
    
    %% slice timing correction
    % TR = 4
    
    disp('>>>Running slice timing correction..')
    cd(rs_path)
    eval(['!',fsl,'slicetimer',' ','-i',' ','201*',' ','-o',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],...
        ' ','-v',' ','-r',' ','4',' ','--odd'])
  
    %% split ASL and BOLD echoes
    % echo 1 = ASL, echo 2 = BOLD
    
    disp('>>>Separating data into BOLD and ASL echoes...')

    cd(rs_path)
    % ASL echo
    eval(['!',fsl,'fslroi',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],' ',[tmp_subj,'_rs_ASL.nii.gz'],' ','0',...
    ' ','140']);
    % BOLD echo
    eval(['!',fsl,'fslroi',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],' ',[tmp_subj,'_rs_BOLD.nii.gz'],' ','140',...
    ' ','140']);
    
    %% delete first 2 volumes
    % allow MR signal to stabilize
    
    disp('>>>Stabilizing MR signal...')

    cd(rs_path)
    eval(['!',fsl,'fslroi',' ','*ASL.nii.gz',' ',[tmp_subj,'_ASL_corr.nii.gz'],' ','2',...
    ' ','138']);
    eval(['!',fsl,'fslroi',' ','*BOLD.nii.gz',' ',[tmp_subj,'_BOLD_corr.nii.gz'],' ','2',...
    ' ','138']);

    %% signal despiking - both echoes
    % use afni????
    
    disp('>>>Running signal despiking sequence...')
    
    eval(['!',afni,'3dDespike',' ','-nomask',' ',[tmp_subj,'_ASL_corr.nii.gz']]);
    % QUESTION: use the -nomask option or default???
    
    eval(['!',afni,'3dAFNItoNIFTI',' ','-prefix',' ',[tmp_subj,'_ASL_despike.nii.gz'],' ',...
        'despike+orig']);
    
    % remove default files
    eval(['!','rm ',' ','despike+orig.BRIK']);
    eval(['!','rm',' ','despike+orig.HEAD']);
    
    eval(['!',afni,'3dDespike',' ','-nomask',' ',[tmp_subj,'_BOLD_corr.nii.gz']]);
    
    eval(['!',afni,'3dAFNItoNIFTI',' ','-prefix',' ',[tmp_subj,'_BOLD_despike.nii.gz'],' ',...
        'despike+orig']);
    
    % remove default files
    eval(['!','rm',' ','despike+orig.BRIK']);
    eval(['!','rm',' ','despike+orig.HEAD']);
    
    %% motion correction
    % rigid body transformations using mcflirt (use volume 1 as template) -
    % least squares approach and 6 parameter spatial transformation.
    
    disp('>>>Running motion correction sequence...')
    
    % asl ts motion correction
    eval(['!',fsl,'mcflirt',' ','-in',' ',[tmp_subj,'_ASL_despike.nii.gz'],' ','-refvol',' ',...
        1,' ','-out',' ',[tmp_subj,'_ASL_mcf.nii.gz']])
    
    % bold ts motion correction
    eval(['!',fsl,'mcflirt',' ','-in',' ',[tmp_subj,'_BOLD_despike.nii.gz'],' ','-refvol',' ',...
        1,' ','-out',' ',[tmp_subj,'_BOLD_mcf'],' ','-plots'])
    
    %% plot the motion parameters
    % parameters will be used later in the pipeline
    
%     data = load([rs_path,[tmp_subj,'_BOLD_mcf.par']]);
% 
%     set(0,'DefaultFigureWindowStyle','docked');
%     
%     for ii = 2:length(data)
%     % calculate mean displacement
%     mean_displacement(ii,:) = (((data(ii,1)-data(ii-1,1))^2)+...
%         ((data(ii,2)-data(ii-1,2))^2)+((data(ii,3)-data(ii-1,3))^2)...
%         +((data(ii,4)-data(ii-1,4))^2)+((data(ii,5)-data(ii-1,5))^2)+...
%         ((data(ii,6)-data(ii-1,6))^2))^(1/2);  
%     end
% 
%     motion_parameters = [data,mean_displacement(:,1)];
%     motion_parameters = num2cell(motion_parameters);
% 
%     rotation_data = cell2mat(motion_parameters(:,1:3));
%     translation_data = cell2mat(motion_parameters(:,4:6));
%     mean_disp = cell2mat(motion_parameters(:,7));
% 
%     figure('Name','MCFLIRT Estimated Rotations','NumberTitle','off')
%     plot(rotation_data,'linewidth',2);
%     title('MCFLIRT Estimated Rotations');
%     xlabel('Volume');
%     ylabel('Rotation (radians)');
%     yline(0,'--');
%     legend({'x','y','z','ref'},'location','northwest');
%     legend boxoff
%     set(gca,'fontsize',20)
% 
%     figure('Name','MCFLIRT Estimated Translations','NumberTitle','off')
%     plot(translation_data,'linewidth',2);
%     title('MCFLIRT Estimated Translation');
%     xlabel('Volume');
%     ylabel('Translation (mm)');
%     yline(0,'--');
%     legend({'x','y','z','ref'},'location','northwest');
%     legend boxoff
%     set(gca,'fontsize',20)
% 
%     figure('Name','MCFLIRT Mean Displacement','NumberTitle','off')
%     plot(mean_disp,'linewidth',2);
%     title('MCFLIRT Mean Displacement');
%     xlabel('Volume');
%     ylabel('Mean Displacement (mm)');
%     yline(0,'--');
%     legend({'mean displacement'},'location','northwest');
%     legend boxoff
%     set(gca,'fontsize',20)
    
    %% create mean image for bold and asl echoes
    % align volumes over time to create mean image for each ts (bold and
    % asl)
    
    disp('>>>Creating mean images from time series...')
    
    % mean asl image
    eval(['!',fsl,'fslmaths',' ',[tmp_subj,'_ASL_mcf.nii.gz'],' ','-Tmean',' ',...
        [tmp_subj,'_ASL_meanvol.nii.gz']])
    
    % mean bold image
    eval(['!',fsl,'fslmaths',' ',[tmp_subj,'_BOLD_mcf.nii.gz'],' ','-Tmean',' ',...
        [tmp_subj,'_BOLD_meanvol.nii.gz']])
    
    %% brain exactraction
    % bet brain extration to remove non-brain tissue in each of the mean
    % images
    
    disp('>>>Performing brain extraction...')
    disp('>>>Creating binary brain masks...')
    
    % mean bold image bet
    eval(['!',fsl,'bet',' ',[tmp_subj,'_BOLD_meanvol.nii.gz'],' ',...
        [tmp_subj,'_BOLD_meanBrain.nii.gz'],' ','-m'])
    
    % mean asl image bet
    eval(['!',fsl,'bet',' ',[tmp_subj,'_ASL_meanvol.nii.gz'],' ',...
        [tmp_subj,'_ASL_meanBrain.nii.gz'],' ','-m']) 
    
    %% =================== BOLD signal reconstruction =====================
    
    %% brain extract the 4D bold ts
    % multiply the 3D binary brain mask over the ts to extract the rest of
    % the volumes
    
    disp('>>>Brain extracting BOLD ts volumes...')
    
    cd(rs_path)
    bold_mcf = load_untouch_nii([tmp_subj,'_BOLD_mcf.nii.gz']);
    bold_img = bold_mcf.img;
    
    bold_mask = load_untouch_nii([tmp_subj,'_BOLD_meanBrain_mask.nii.gz']);
    mask_img = bold_mask.img;
    mask_img = single(mask_img);
    
    clear bold_ts
    
    for jj = 1:size(bold_img,4)
        tmp_vol = bold_img(:,:,:,jj);
        bet_vol = tmp_vol.*mask_img;
        bold_ts(:,:,:,jj) = bet_vol;
    end
    
    dyn = size(bold_ts,4);
    
    % save it
    nii_hdr = bold_mcf.hdr; % header isn't changed
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img = bold_ts; % 4D data
    save_nii(img_nii.nii,[rs_path,[tmp_subj,'_BOLD_ts.nii.gz']]);
    
    %% window average the bold signal data
    % use a sliding window average technique
    
    disp('>>>Window averaging BOLD ts...')
    
    bold_ts_nii = load_untouch_nii([rs_path,[tmp_subj,'_BOLD_ts.nii.gz']]);
    dyn = size(bold_ts,4); 

    clear echo_ctr echo_lab
    for ii = 1:dyn/2
        echo_ctr(:,:,:,ii) = bold_ts(:,:,:,2*ii-1);
        echo_lab(:,:,:,ii) = bold_ts(:,:,:,2*ii);
    end

    % BOLD - sliding window average - addition/subtraction schemes can be
    % modified.

    clear BOLD_e2ts
    BOLD_e2ts = zeros(size(bold_ts));
    BOLD_e2ts(:,:,:,1) = (echo_ctr(:,:,:,1)+echo_lab(:,:,:,1))/2; % find mean between echo_ctr and echo_lab
    BOLD_e2ts(:,:,:,end) = (echo_ctr(:,:,:,end)+echo_lab(:,:,:,end))/2;

    for ii = 1:dyn/2-1
        BOLD_e2ts(:,:,:,2*ii+1) = (echo_ctr(:,:,:,ii)+echo_lab(:,:,:,ii))/2;
        BOLD_e2ts(:,:,:,2*ii) = (echo_lab(:,:,:,ii)+echo_ctr(:,:,:,ii+1))/2;
    end

    % save it
    nii_hdr = bold_ts_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=BOLD_e2ts; % BOLD_e2ts is the 
    save_nii(img_nii.nii,[rs_path,[tmp_subj,'_BOLD_window_avg.nii.gz']]);
    
    %% brain extract T1 anat image
    
    disp('>>>Running brain extraction for MPRAGE anat image...')
    
    anat_path = [subj_path,tmp_subj,'/anat/'];
    anat_file = dir([anat_path,'co*']).name;
    
    cd(anat_path)

    eval(['!',fsl,'bet',' ',anat_file,' ',[tmp_subj,'_anat_brain.nii.gz'],' ',...
        '-f',' ','0.65',' ','-g',' ','-0.5'])
    
    
    %% co-registration with anatomical MPRAGE image
    % use epi_reg and boundary-based registration in FSL
    % corrects for motion differences between the anatomical and functional
    % scans
    
    anat_wholehead = [anat_path,anat_file];
    anat_brain = [anat_path,[tmp_subj,'_anat_brain.nii.gz']];
    wm_seg = [anat_path,[tmp_subj,'_anat_brain_seg_2.nii.gz']];
    cd(rs_path)
    copyfile(anat_wholehead)
    copyfile(anat_brain)
    copyfile(wm_seg)
    
    disp('>>>Registering functional data to MNI152 space (via hi-res structural scan)...')
    
    % linear registration of func data to structural image
    eval('!',fsl,'flirt',' ','-v',' ','-ref',' ',[tmp_subj,'_anat_brain.nii.gz'],' ',...
    '-in',' ',[tmp_subj,'_BOLD_window_avg.nii.gz'],' ','-dof 6',' ','-omat',' ',...
    [tmp_subj,'_func2struct.mat'])
    
    % linear registration of structural image to MNI space
    eval(['!',fsl,'flirt',' ','-ref $FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '-in',' ',[tmp_subj,'_anat_brain.nii.gz'],' ','-omat',' ',[tmp_subj,'_T1_2_MNI.mat']])
    
    % non-linear registration of structural image to MNI space
    eval(['!',fsl,'fnirt',' ','--in=',anat_file,' ','--aff=',[tmp_subj,'_T1_2_MNI.mat'],...
    ' ','--cout=',[tmp_subj,'_T1_2_MNI_nonlinear'],' ','--config=T1_2_MNI152_2mm']])
    
    % apply transformation to 4d func data
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_BOLD_window_avg.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--premat=',[tmp_subj,'_func2struct.mat'],' ','--out=',[tmp_subj,'_BOLD_MNI'])
    
    disp(' ')
    disp('registration complete.')
    
end



disp(' ')
disp('COMPLETE')
