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
    processdir = [subj_path,tmp_subj,'/rs/norm/'];
    
    %% slice timing correction
    % TR = 4
    
    disp('>>>Running slice timing correction..')
    cd(processdir)
    eval(['!',fsl,'slicetimer',' ','-i',' ','201*',' ','-o',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],...
        ' ','-v',' ','-r',' ','4',' ','--odd'])
  
    %% split ASL and BOLD echoes
    % echo 1 = ASL, echo 2 = BOLD
    
    disp('>>>Separating data into BOLD and ASL echoes...')

    cd(processdir)
    % ASL echo
    eval(['!',fsl,'fslroi',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],' ',[tmp_subj,'_rs_ASL.nii.gz'],' ','0',...
    ' ','140']);
    % BOLD echo
    eval(['!',fsl,'fslroi',' ',[tmp_subj,'_RS_slicetimer.nii.gz'],' ',[tmp_subj,'_rs_BOLD.nii.gz'],' ','140',...
    ' ','140']);
    
    %% delete first 2 volumes
    % allow MR signal to stabilize
    
    disp('>>>Stabilizing MR signal...')

    cd(processdir)
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
        1,' ','-out',' ',[tmp_subj,'_ASL_mcf'],' -plots'])
    
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
    
    cd(processdir)
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
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_BOLD_ts.nii.gz']]);
    
    %% window average the bold signal data, parse control and label volumes
    % use a sliding window average technique
    
    disp('>>>Window averaging BOLD ts...')
    
    % import bold to sum control and label volumes
    bold_ts_nii = load_untouch_nii([processdir,[tmp_subj,'_BOLD_ts.nii.gz']]);
    dyn = size(bold_ts,4);
    
    % parse the control and label volumes
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
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_BOLD_window_avg.nii.gz']]);
    
    %% brain extract T1 anat image
    
    disp('>>>Running brain extraction for MPRAGE anat image...')
    
    anat_path = [subj_path,tmp_subj,'/anat/'];
    anat_file = dir([anat_path,'co*']).name;
    
    cd(anat_path)

    eval(['!',fsl,'bet',' ',anat_file,' ',[tmp_subj,'_anat_brain.nii.gz'],' ',...
        '-f',' ','0.65',' ','-g',' ','-0.5'])
    
    
    %% co-registration and transformation into MNI space (2mm)
    
    anat_wholehead = [anat_path,anat_file];
    anat_brain = [anat_path,[tmp_subj,'_anat_brain.nii.gz']];
    wm_seg = [anat_path,[tmp_subj,'_anat_brain_seg_2.nii.gz']];
    cd(processdir)
    copyfile(anat_wholehead)
    copyfile(anat_brain)
    copyfile(wm_seg)
    
    %%
    
    disp('>>>Registering functional data to MNI152 space (via hi-res structural scan)...')
    
    disp('  >>>Running FLIRT...')

    % linear registration of func data to structural image
    eval(['!',fsl,'flirt',' ','-ref',' ',[tmp_subj,'_anat_brain.nii.gz'],' ',...
    '-in',' ',[tmp_subj,'_BOLD_window_avg.nii.gz'],' ','-dof 6',' ','-omat',' ',...
    [tmp_subj,'_func2struct.mat']])
    
    % linear registration of structural image to MNI space
    eval(['!',fsl,'flirt',' ','-ref $FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '-in',' ',[tmp_subj,'_anat_brain.nii.gz'],' ','-omat',' ',[tmp_subj,'_T1_2_MNI.mat']])

    disp('  >>>Running FNIRT...')
    
    % non-linear registration of structural image to MNI space
    eval(['!',fsl,'fnirt',' ','--in=',anat_file,' ','--aff=',[tmp_subj,'_T1_2_MNI.mat'],...
    ' ','--cout=',[tmp_subj,'_T1_2_MNI_nonlinear'],' ','--config=T1_2_MNI152_2mm'])


%%
    disp('  >>>Applying transformation to BOLD ts...')
    
    % apply transformation to 4d func data
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_BOLD_window_avg.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--premat=',[tmp_subj,'_func2struct.mat'],' ','--out=',[tmp_subj,'_BOLD_MNI']])


    %% =================== prepare 4D ASL dataset ========================
    
    disp('>>>Preparing ASL dataset')
    
    %% use asl_bet mean image to brain extract the 4d asl ts
    
    disp(' >>Brain extracting 4D ASL data...')

    working_dir = processdir;
    ts_img_name = [tmp_subj,'_ASL_mcf.nii.gz'];
    mask_img_name = [tmp_subj,'_ASL_meanBrain_mask.nii.gz'];
    save_name = [tmp_subj,'_ASL_ts.nii.gz'];
    
    asl_ts = ts_bet(working_dir,ts_img_name,mask_img_name,save_name);
    
    %% run surround subtraction on ASL data
    % computes difference volumes 
    
    disp(' >>Running surround substraction...')

    % import bold to sum control and label volumes
    asl_nii = load_untouch_nii([processdir,[tmp_subj,'_ASL_ts.nii.gz']]);
    asl_data = double(asl_nii.img);
    
    dyn = size(asl_data,4);
    
    % surround subtraction of control and labelled images:
    % difference between each image and the average of its 2 nearest
    % neighbours is formed - reduces transient artifacts of BOLD weighting
    
    for ii = 2:size(asl_data,4)-1
        if mod(ii,2) == 0 % even
           diff_vols(:,:,:,ii) = asl_data(:,:,:,ii) - ((asl_data(:,:,:,ii-1) + asl_data(:,:,:,ii+1))./2); % (control) - (avg of surrounding tags)
        else % odd
           diff_vols(:,:,:,ii) = -asl_data(:,:,:,ii) + ((asl_data(:,:,:,ii-1) + asl_data(:,:,:,ii+1))./2); % -(tag) + (avg of surrounding controls)
        end
    end
    
    % fill in dummy data on either end
    diff_vols(:,:,:,1) = diff_vols(:,:,:,2);
    diff_vols(:,:,:,end+1) = diff_vols(:,:,:,end);
    
    % save it
    nii_hdr = asl_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=diff_vols; % 4d data
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_ASL_sub.nii.gz']]);
    
    %% compute voxelwise CBF0 map
    % index of basal vascular tension, with partial volume and T2*
    % corrections.
    
    cd(processdir)
    
    % path to M0 image
    m0_path = '/Users/kaden_shearer/Dropbox/kaden_thesis/movie_data/varsity_data/controls/nci.cq.lcf.8/m0/norm/';
    m0_file = [m0_path,'*.nii.gz'];
    
    %% FAST segmentation
    
    disp(' >>FAST segmenting anatomical image...')
    eval(['!',fsl,'fast -g -v ',[tmp_subj,'_anat_brain.nii.gz']])
    
   %% adjust for varying PLD due to 2D EPI readout
   % adjusted on a slice-by-slice basis using PLD = 1000 ms + (St)*(slice-1)
   % St = slice time correction factor in ms (53.8 ms)
   
%    asl_slices = transpose(1:138);
%    for xx = 1:length(asl_slices)
%        tmp_slice = asl_slices(xx);
%        slice_PLD(xx) = (1000 + (53.8*(tmp_slice - 1)))*0.001;
%    end
%    
%    cmd_line_PLD = '1';
%    
%    for zz = 2:length(slice_PLD)
%        tmp_PLD = slice_PLD(zz);
%        cmd_line_PLD = [cmd_line_PLD,',',num2str(tmp_PLD)];
%    end

% >>>  NOTE: Do you need to correct for PLD on slice-by-slice basis if
% slice timing  correction was already performed on the raw data?
   
   %% run oxford_asl with partial volume and T2* corrections
   
   cd(processdir)
   
   disp(' >>Running Oxford_ASL with pvcorr and T2*corr...')
   eval(['!',fsl,'oxford_asl -i ',[tmp_subj,'_ASL_sub.nii.gz'],' -o ',...
       [tmp_subj,'_output_cbf'],' -m ',[tmp_subj,'_ASL_meanBrain_mask.nii.gz'],...
       ' --tis=2.665 --bolus=1.665 --casl -r ',[tmp_subj,'_anat_brain.nii.gz'],...
       ' --te=10 -c ',m0_file,' --cmethod=voxel ','--alpha=',num2str(0.84),' ',...
       ' --t2star --echospacing=0.00047 --pvcorr --pvgm=',[tmp_subj,'_anat_brain_pve_1.nii.gz'],...
       ' --pvwm=',[tmp_subj,'_anat_brain_pve_2.nii.gz']])
   
   % perfusion_calib.nii.gz in native space dir is flow map in absolute
   % units (mL/100g/min)
   
   cd(processdir)
   native_space_path = '/Users/kaden_shearer/Dropbox/kaden_thesis/movie_data/varsity_data/controls/nci.cq.lcf.8/rs/norm/nci.cq.lcf.8_output_cbf/native_space/';
   copyfile([native_space_path,'perfusion_calib.nii.gz'])
   eval(['!','mv perfusion_calib.nii.gz ',[tmp_subj,'_ASL_perfusion_map.nii.gz']])
   
   %% resample 4d ASL data into 2mm MNI space
   % use linear and non-linear matricies calculated before for bold data
   
   disp(' >> Resampling 4D ASL data into 2 mm MNI space...')
   
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_ASL_sub.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--premat=',[tmp_subj,'_func2struct.mat'],' ','--out=',[tmp_subj,'_ASL_MNI']])

    %% resample ASL perfusion map into 2mm MNI space
    
    disp(' >> Resampling ASL perfusion map into 2 mm MNI space...')
   
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_ASL_perfusion_map.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--premat=',[tmp_subj,'_func2struct.mat'],' ','--out=',[tmp_subj,'_ASL_perfusion_map_MNI']])

    %% resample 3D segmented GM, WM, and CSF masks into 2mm MNI space
    % use concatenated linear and non-linear transformation mat generated
    % previously in the pipeline
    
    disp(' >> Resampling FAST segments into 2 mm MNI space...')
   
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_anat_brain_seg_0.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[tmp_subj,'_anat_brain_seg_0_MNI.nii.gz']])
    
    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_anat_brain_seg_1.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[tmp_subj,'_anat_brain_seg_1_MNI.nii.gz']])

    eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[tmp_subj,'_anat_brain_seg_2.nii.gz'],' ','--warp=',[tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[tmp_subj,'_anat_brain_seg_2_MNI.nii.gz']])
    
    %% ===== Denoise 4d bold data using fsl_regfilt =====
    % temporal regression of nuisance parameters such as cardiac/resp noise
    % and non-neuronally related signal
    
    disp(' >> Denoising BOLD data via temporal regression...')
    
    % load 4D bold data in registered in MNI space
    bold_nii = load_untouch_nii([processdir,[tmp_subj,'_BOLD_MNI.nii.gz']]);
    bold_img = double(bold_nii.img);
    
    dyn = size(bold_img,4);
    
    % extract global signal (Murphy and Fox, 2017)
    % mean of the voxel time series within the brain; 
    bold_avg_ts = make_avg_ts(bold_img);
    
    %% mean WM signal
    % WM and CSF signal is used as a physiological proxy for cardiac/resp noise
    % since they have minimal contribution to neural activation.
    % load WM segment
    
    disp(' >> Removing cardiac & respiratory artifacts')
    
    wm_seg_nii = load_untouch_nii([processdir,[tmp_subj,'_anat_brain_seg_2_MNI.nii.gz']]);
    wm_seg_img = double(wm_seg_nii.img);
        
    % apply segment to bold data
    for yy = 1:size(bold_img,4)
        tmp_vol = bold_img(:,:,:,yy);
        wm_sig(:,:,:,yy) = tmp_vol.*wm_seg_img;
    end
    
    % save it
    nii_hdr = bold_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=wm_sig; % 4d data
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_BOLD_WM_sig.nii.gz']]);
    
    wm_avg_ts = make_avg_ts(wm_sig);
    
    %% mean CSF signal
    % load CSF segment
    csf_seg_nii = load_untouch_nii([processdir,[tmp_subj,'_anat_brain_seg_0_MNI.nii.gz']]);
    csf_seg_img = double(csf_seg_nii.img);
    
    dyn = size(bold_img,4);
    
    % apply segment to bold data
    for yy = 1:size(bold_img,4)
        tmp_vol = bold_img(:,:,:,yy);
        csf_sig(:,:,:,yy) = tmp_vol.*csf_seg_img;
    end
    
    % save it
    nii_hdr = bold_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=csf_sig; % 4d data
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_BOLD_CSF_sig.nii.gz']]);
    
    csf_avg_ts = make_avg_ts(csf_sig);
    
    %% extract six rigid-body motion parameters 
    
    data = load([rs_path,[tmp_subj,'_BOLD_mcf.par']]);
    
    for ii = 2:length(data)
    % calculate mean displacement
    mean_displacement(ii,:) = (((data(ii,1)-data(ii-1,1))^2)+...
        ((data(ii,2)-data(ii-1,2))^2)+((data(ii,3)-data(ii-1,3))^2)...
        +((data(ii,4)-data(ii-1,4))^2)+((data(ii,5)-data(ii-1,5))^2)+...
        ((data(ii,6)-data(ii-1,6))^2))^(1/2);  
    end

    motion_parameters = [data,mean_displacement(:,1)];
    motion_parameters = num2cell(motion_parameters);

    rotation_data = cell2mat(motion_parameters(:,1:3));
    translation_data = cell2mat(motion_parameters(:,4:6));
    mean_disp = cell2mat(motion_parameters(:,7));
    
    %% create matrix with time courses of nuisance parameters
    noise_mat(:,1) = bold_avg_ts;
    noise_mat(:,2) = wm_avg_ts;
    noise_mat(:,3) = csf_avg_ts;
    noise_mat(:,4:6) = translation_data;
    noise_mat(:,7:9) = rotation_data;
    
    % create design matrix as input for fsl_regfilt
    dlmwrite('bold_regfilt_design.txt',noise_mat,'delimiter',' ');
    
    %% use fsl_regfilt to regress out extracted nuisance parameters
    eval(['!',fsl,'fsl_regfilt -i ',[tmp_subj,'_BOLD_MNI.nii.gz'],...
        ' -d bold_regfilt_design.txt -f "1,2,3,4,5,6,7,8,9" -o ',...
        [tmp_subj,'_BOLD_denoised.nii.gz'],' -v'])
    disp('fsl_regfilt COMPLETE')
    
    %% normalize denoised BOLD signal
    % convert absolute signal to %change in BOLD relative to baseline (avg
    % of first 10 volumes)
    
    image_path = [processdir,[tmp_subj,'_BOLD_denoised.nii.gz']];
    bsl_range = 1:10;
    save_dir = processdir;
    save_name = [tmp_subj,'_BOLD_normalized.nii.gz'];
    bold_normalized = normalize_signal(image_path,bsl_range,save_dir,save_name);
    
     %% check avg signal of denoised data
%     
%     denoised_ts = make_avg_ts(bold_normalized);
%     mni_ts = make_avg_ts(mni_normalized);
%     
%     plot(mni_ts,'linewidth',2)
%     hold on
%     plot(denoised_ts,'linewidth',2)
%     set(gca,'fontsize',14)
%     xlabel('ts volume')
%     ylabel('avg BOLD signal')
%     legend('raw BOLD','denoised BOLD')

    %% ===== Denoise 4d asl data using fsl_regfilt =====
    
    disp(' >> Denoising ASL data via temporal regression...')
    
    % load 4d asl data (resampled into MNI space)
    asl_nii = load_untouch_nii([processdir,[tmp_subj,'_ASL_MNI.nii.gz']]);
    asl_img = double(asl_nii.img);
    
    dyn = size(asl_img,4);
    
    % extract global signal from asl ts
    asl_avg_ts = make_avg_ts(asl_img);
    
    % extract avg wm signal from asl data (seg is the same from MNI)
    for yy = 1:size(asl_img,4)
        tmp_vol = asl_img(:,:,:,yy);
        wm_sig(:,:,:,yy) = tmp_vol.*wm_seg_img;
    end
    
    % save it
    nii_hdr = asl_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=wm_sig; % 4d data
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_ASL_WM_sig.nii.gz']]);
    
    asl_wm_avg_ts = make_avg_ts(wm_sig); % avg wm signal

    % extract avg csf signal from asl data
    for yy = 1:size(asl_img,4)
        tmp_vol = asl_img(:,:,:,yy);
        csf_sig(:,:,:,yy) = tmp_vol.*csf_seg_img;
    end
    
    % save it
    nii_hdr = asl_nii.hdr;
    nii_hdr.dime.dim(5)=dyn; % # of volumes
    img_nii.nii.hdr = nii_hdr;
    img_nii.nii.img=csf_sig; % 4d data
    save_nii(img_nii.nii,[processdir,[tmp_subj,'_ASL_CSF_sig.nii.gz']]);
    
    asl_csf_avg_ts = make_avg_ts(csf_sig); % avg csf signal 
    
    % rigid-body motion parameters same as bold data
    
    % create matrix with time courses of nuisance parameters
    asl_noise_mat(:,1) = asl_avg_ts;
    asl_noise_mat(:,2) = asl_wm_avg_ts;
    asl_noise_mat(:,3) = asl_csf_avg_ts;
    asl_noise_mat(:,4:6) = translation_data;
    asl_noise_mat(:,7:9) = rotation_data;
    
    % create design matrix as input for fsl_regfilt
    dlmwrite('asl_regfilt_design.txt',asl_noise_mat,'delimiter',' ');
    
    % use fsl_regfilt to regress out extracted nuisance parameters
    eval(['!',fsl,'fsl_regfilt -i ',[tmp_subj,'_ASL_MNI.nii.gz'],...
        ' -d asl_regfilt_design.txt -f "1,2,3,4,5,6,7,8,9" -o ',...
        [tmp_subj,'_ASL_denoised.nii.gz'],' -v'])
    disp('fsl_regfilt COMPLETE')
    
    %% normalize the asl signal
    
    image_path = [processdir,[tmp_subj,'_ASL_denoised.nii.gz']];
    bsl_range = 1:10;
    save_dir = processdir;
    save_name = [tmp_subj,'_ASL_normalized.nii.gz'];
    asl_normalized = normalize_signal(image_path,bsl_range,save_dir,save_name);
    
    %% calculate mean displacement for 4d data
    
    disp(' >>Calculating mean displacement of 4d data...')
    
    % calculate mean displacement for bold data
    bold_motion = load([rs_path,[tmp_subj,'_BOLD_mcf.par']]);
    for ii = 2:length(bold_motion)
    % calculate mean displacement
    bold_mean_displacement(ii,:) = (((bold_motion(ii,1)-bold_motion(ii-1,1))^2)+...
        ((bold_motion(ii,2)-bold_motion(ii-1,2))^2)+((bold_motion(ii,3)-bold_motion(ii-1,3))^2)...
        +((bold_motion(ii,4)-bold_motion(ii-1,4))^2)+((bold_motion(ii,5)-bold_motion(ii-1,5))^2)+...
        ((bold_motion(ii,6)-bold_motion(ii-1,6))^2))^(1/2);  
    end
    
    % calculate mean displacment for asl data
    asl_motion = load([rs_path,[tmp_subj,'_ASL_mcf.par']]);
    for ii = 2:length(asl_motion)
    % calculate mean displacement
    asl_mean_displacement(ii,:) = (((asl_motion(ii,1)-asl_motion(ii-1,1))^2)+...
        ((asl_motion(ii,2)-asl_motion(ii-1,2))^2)+((asl_motion(ii,3)-asl_motion(ii-1,3))^2)...
        +((asl_motion(ii,4)-asl_motion(ii-1,4))^2)+((asl_motion(ii,5)-asl_motion(ii-1,5))^2)+...
        ((asl_motion(ii,6)-asl_motion(ii-1,6))^2))^(1/2);
    end
    
%     figure('Name','BOLD Mean Displacement','NumberTitle','off')
%     plot(bold_mean_displacement,'linewidth',2);
%     title('BOLD Mean Displacement');
%     xlabel('Volume');
%     ylabel('Mean Displacement (mm)');
%     yline(0,'--');
%     legend({'mean displacement'},'location','northwest');
%     legend boxoff
%     set(gca,'fontsize',20)
%     
%     figure('Name','ASL Mean Displacement','NumberTitle','off')
%     plot(asl_mean_displacement,'linewidth',2);
%     title('ASL Mean Displacement');
%     xlabel('Volume');
%     ylabel('Mean Displacement (mm)');
%     yline(0,'--');
%     legend({'mean displacement'},'location','northwest');
%     legend boxoff
%     set(gca,'fontsize',20)
      
      %% === bold volume censoring (mean displacement > 0.2 mm) ===
      
      disp(' >>Volume censoring bold data (> 0.2mm displacement)')
      
      % censor out bold volumes with displacement > 0.2 mm
      bold_censor_locs = find(bold_mean_displacement > 0.2);
      
      counter = 1;
      
      for tt = 1:length(bold_normalized)
          tmp_vol = bold_normalized(:,:,:,tt);
          censor = sum(tt == bold_censor_locs);
          if censor == 0
              bold_censored_vols(:,:,:,counter) = tmp_vol;
              counter = counter + 1; 
          end
      end
      
      tmp_nii = load_untouch_nii([processdir,[tmp_subj,'_BOLD_normalized.nii.gz']]);
      dyn = size(bold_censored_vols,4);
      
      nii_hdr = tmp_nii.hdr;
      nii_hdr.dime.dim(5)=dyn; % # of volumes
      img_nii.nii.hdr = nii_hdr;
      img_nii.nii.img=bold_censored_vols; % 4d data
      save_nii(img_nii.nii,[processdir,[tmp_subj,'_BOLD_volcensored.nii.gz']]);
      
       %% === asl volume censoring (mean displacement > 0.2 mm) ===
      
      disp(' >>Volume censoring asl data (> 0.2mm displacement)')
      
      % censor out bold volumes with displacement > 0.2 mm
      asl_censor_locs = find(asl_mean_displacement > 0.2);
      
      counter = 1;
      
      for tt = 1:length(asl_normalized)
          tmp_vol = asl_normalized(:,:,:,tt);
          censor = sum(tt == asl_censor_locs);
          if censor == 0
              asl_censored_vols(:,:,:,counter) = tmp_vol;
              counter = counter + 1; 
          end
      end
      
      tmp_nii = load_untouch_nii([processdir,[tmp_subj,'_ASL_normalized.nii.gz']]);
      dyn = size(asl_censored_vols,4);
      
      nii_hdr = tmp_nii.hdr;
      nii_hdr.dime.dim(5)=dyn; % # of volumes
      img_nii.nii.hdr = nii_hdr;
      img_nii.nii.img=asl_censored_vols; % 4d data
      save_nii(img_nii.nii,[processdir,[tmp_subj,'_ASL_volcensored.nii.gz']]);
      
      %% smooth bold dataset
      % use 6mm full width at half-maximum (FWHM) Gaussian kernel.
      
      
    
end

disp(' ')
disp('COMPLETE')
