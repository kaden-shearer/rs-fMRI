%% Thesis Pipeline
% scripted by kaden on 11-12-2019
% used to process pCASL data and assess BOLD-CBF coupling between rs and
% movie datasets.

% >> Note: current pipeline only examines control subjects. Update
% directories for concussed subjects in future iterations.

%% ======================= Create Directories ============================

% paths to processing software
fsleyes = '/usr/local/fsl/bin/fsleyes';
fsl = '/usr/local/fsl/bin/';
afni = '/Users/kaden_shearer/abin/';
% control subject directory
subj_path = '/Users/kaden_shearer/Dropbox/kaden_thesis/data/varsity_data/controls/';
subj_dir = dir([subj_path,'n*']);

%% ========================= RS Pre-Processing =============================

for ii = 1%1:length(subj_dir)
    
close all
clc
clearvars -except fsleyes fsl afni subj_path subj_dir ii

tmp_subj = subj_dir(ii).name;
disp(['============== Working on ',tmp_subj,' RS data... =================='])
processdir = [subj_path,tmp_subj,'/rs/norm/'];
cd(processdir)

%% ========================= pCASL Processing =============================

disp(['>> Preprocessing pCASL data for ',tmp_subj])

%% Slice-timing correction
% TR = 4s

msg = sprintf('Running slice-timing correction...'); n=numel(msg); fprintf(msg)
eval(['!',fsl,'slicetimer',' ','-i',' ','201*',' ','-o',' ',[processdir,tmp_subj,'_pcasl_slicetimer.nii.gz'],...
    ' ','-r',' ','4',' ','--odd'])
fprintf(repmat('\b',1,n)); fprintf([msg,'complete']); disp(' ');

%% Separate BOLD and ASL echoes
% echo1 = ASL, echo2 = BOLD
fprintf('Separating BOLD and ASL echoes...'); disp(' ');

slicetimer_img = load_image([processdir,tmp_subj,'_pcasl_slicetimer.nii.gz']);

tmp_length = size(slicetimer_img,4);
thresh = tmp_length/2;


%%
% ASL data
eval(['!',fsl,'fslroi',' ',[processdir,tmp_subj,'_pcasl_slicetimer.nii.gz'],...
    ' ',[processdir,tmp_subj,'_pcasl_ASL.nii.gz'],' ','0',' ',num2str(thresh)]);
% BOLD data
eval(['!',fsl,'fslroi',' ',[processdir,tmp_subj,'_pcasl_slicetimer.nii.gz'],...
    ' ',[processdir,tmp_subj,'_pcasl_BOLD.nii.gz'],' ',num2str(thresh),' ',num2str(thresh)]);
fprintf('complete'); disp(' ');

%% ensure echoes are separated correctly

eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_pcasl_slicetimer.nii.gz',' &']])

choice = menu(['Assess Echo Separation (Slice ',num2str(thresh),'). Continue?'],'Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    disp('Echo Separation ERROR')
    return
end

%% Delete first 2 volumes in each echo
% ensures that the MR signal is stabilized and  subject is acclimitized to
% loud scanner noises

msg = sprintf('Ensuring MR signal stabilization...'); n=numel(msg); fprintf(msg)
eval(['!',fsl,'fslroi',' ',[processdir,tmp_subj,'_pcasl_ASL.nii.gz'],...
    ' ',[processdir,tmp_subj,'_ASL_corr.nii.gz'],' ','2',' ',num2str(thresh-2)]);
eval(['!',fsl,'fslroi',' ',[processdir,tmp_subj,'_pcasl_BOLD.nii.gz'],...
    ' ',[tmp_subj,'_BOLD_corr.nii.gz'],' ','2',' ',num2str(thresh-2)]);
fprintf(repmat('\b',1,n)); fprintf([msg,'complete']); disp(' ');

%% Signal despiking

msg = sprintf('Running signal despiking...'); n=numel(msg); fprintf(msg)

% despike ASL
eval(['!',afni,'3dDespike',' ','-nomask',' ',[processdir,tmp_subj,'_ASL_corr.nii.gz']]);
eval(['!',afni,'3dAFNItoNIFTI',' ','-prefix',' ',[processdir,tmp_subj,'_ASL_despike.nii.gz'],' ',...
        'despike+orig']);
eval(['!','rm ',' ','despike+orig.BRIK']);
eval(['!','rm',' ','despike+orig.HEAD']);

% despike BOLD
eval(['!',afni,'3dDespike',' ','-nomask',' ',[processdir,tmp_subj,'_BOLD_corr.nii.gz']]);
eval(['!',afni,'3dAFNItoNIFTI',' ','-prefix',' ',[processdir,tmp_subj,'_BOLD_despike.nii.gz'],' ',...
        'despike+orig']);
eval(['!','rm',' ','despike+orig.BRIK']);
eval(['!','rm',' ','despike+orig.HEAD']);

fprintf(repmat('\b',1,n)); fprintf([msg,'complete']); disp(' ');

%% check despiking

figure('name','BOLD despike','numbertitle','off')
plot_image_ts([processdir,tmp_subj,'_BOLD_corr.nii.gz']); hold on
plot_image_ts([processdir,tmp_subj,'_BOLD_despike.nii.gz']);
legend('raw','despiked');

figure('name','ASL despike','numbertitle','off')
plot_image_ts([processdir,tmp_subj,'_ASL_corr.nii.gz']); hold on
plot_image_ts([processdir,tmp_subj,'_ASL_despike.nii.gz']);
legend('raw','despiked');

choice = menu('Assess Signal Despiking. Continue?','Yes','No');
if choice==1 || choice==0
   close all
elseif choice==2
    close all
    disp('Signal Despiking ERROR')
    return
end


%% MCFLIRT motion correction
% use first volume as template, least square approach and 6 parameter
% spatial transformation (rigid-body)

msg = sprintf('Running MCFLIRT motion correction...'); n=numel(msg); fprintf(msg)

eval(['!',fsl,'mcflirt',' ','-in',' ',[processdir,tmp_subj,'_ASL_despike.nii.gz'],' ','-refvol',' ',...
        '1',' ','-out',' ',[processdir,tmp_subj,'_ASL_mcf'],' -plots'])
eval(['!',fsl,'mcflirt',' ','-in',' ',[processdir,tmp_subj,'_BOLD_despike.nii.gz'],' ','-refvol',' ',...
        '1',' ','-out',' ',[processdir,tmp_subj,'_BOLD_mcf'],' ','-plots'])

fprintf(repmat('\b',1,n)); fprintf([msg,'complete']); disp(' ');

%% check motion parameters

par_file = [processdir,tmp_subj,'_BOLD_mcf.par'];
figure('name','BOLD MCF','numbertitle','off')
[~,~,~] = plot_motion_parameters(par_file);

par_file = [processdir,tmp_subj,'_ASL_mcf.par'];
figure('name','ASL MCF','numbertitle','off')
[~,~,~] = plot_motion_parameters(par_file);

choice = menu('Assess Motion Parameters. Continue?','Yes','No');
if choice==1 || choice==0
   close all
elseif choice==2
    close all
    disp('Motion Correction ERROR')
    return
end

%% Brain extract ts
% create mean image, bet mean image, create brain mask, multiply mask
% through the ts to extract brain tissue

fprintf('Brain extracting 4D data...'); disp(' ');

% generate mean images
msg = sprintf('   Generating mean images...'); n=numel(msg); fprintf(msg)
eval(['!',fsl,'fslmaths',' ',[processdir,tmp_subj,'_ASL_mcf.nii.gz'],' ','-Tmean',' ',...
        [processdir,tmp_subj,'_ASL_meanvol.nii.gz']])
eval(['!',fsl,'fslmaths',' ',[processdir,tmp_subj,'_BOLD_mcf.nii.gz'],' ','-Tmean',' ',...
        [processdir,tmp_subj,'_BOLD_meanvol.nii.gz']])
fprintf('complete'); disp(' ');

% brain extract mean images, create binary masks
fprintf('   Running BET brain extraction...'); disp(' ');

%% BOLD BET
fprintf('     Running BOLD BET...')
% initial bet, then edit f and g parameters if required
eval(['!',fsl,'bet',' ',[processdir,tmp_subj,'_BOLD_meanvol.nii.gz'],' ',...
        [processdir,tmp_subj,'_BOLD_meanBrain.nii.gz'],' ','-m'])    
eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_BOLD_meanvol.nii.gz'],' ',...
    [processdir,tmp_subj,'_BOLD_meanBrain.nii.gz'],' --cmap blue-lightblue --alpha 50 &'])
choice = menu('Assess BET Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    defvalues = {'0','0.5'};
    while(1)
        prompt = {'Enter vertical gradient parameter (g; -1:1):','fractional intensity parameter (f; 0:1):'};
        dlgtitle = 'Bet Input';
        dims = [1,75];
        definput = defvalues;
        BetInput = inputdlg(prompt,dlgtitle,dims,definput);
        defvalues = {char(BetInput(1,1)),char(BetInput(2,1))};
        
        eval(['!',fsl,'bet',' ',[processdir,tmp_subj,'_BOLD_meanvol.nii.gz'],' ',...
        [processdir,tmp_subj,'_BOLD_meanBrain.nii.gz'],' ','-m -f ',char(BetInput(2,1)),...
        ' -g ',char(BetInput(1,1))])
        
        eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_BOLD_meanvol.nii.gz'],' ',...
            [processdir,tmp_subj,'_BOLD_meanBrain.nii.gz'],' --cmap blue-lightblue --alpha 50 &'])
        
        choice = menu('Assess BET Quality. Continue?','Yes','No');
        if choice==1 || choice==0
            eval(['!','pkill fsleyes']);
            break;
        end
        eval(['!','pkill fsleyes'])
    end
end

eval(['!','pkill fsleyes'])
fprintf('complete'); disp(' ')
   
%% ASL BET
fprintf('     Running ASL BET...')

eval(['!',fsl,'bet',' ',[processdir,tmp_subj,'_ASL_meanvol.nii.gz'],' ',...
        [processdir,tmp_subj,'_ASL_meanBrain.nii.gz'],' ','-m'])
    
eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_ASL_meanvol.nii.gz'],' ',...
    [processdir,tmp_subj,'_ASL_meanBrain.nii.gz'],' --cmap blue-lightblue --alpha 50 &'])

choice = menu('Assess BET Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    defvalues = {'0','0.5'};
    while(1)
        prompt = {'Enter vertical gradient parameter (g; -1:1):','fractional intensity parameter (f; 0:1):'};
        dlgtitle = 'Bet Input';
        dims = [1,75];
        definput = defvalues;
        BetInput = inputdlg(prompt,dlgtitle,dims,definput);
        defvalues = {char(BetInput(1,1)),char(BetInput(2,1))};
        
        eval(['!',fsl,'bet',' ',[processdir,tmp_subj,'_ASL_meanvol.nii.gz'],' ',...
        [processdir,tmp_subj,'_ASL_meanBrain.nii.gz'],' ','-m -f ',char(BetInput(2,1)),...
        ' -g ',char(BetInput(1,1))])
        
        eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_ASL_meanvol.nii.gz'],' ',...
            [processdir,tmp_subj,'_ASL_meanBrain.nii.gz'],' --cmap blue-lightblue --alpha 50 &'])
        
        choice = menu('Assess BET Quality. Continue?','Yes','No');
        if choice==1 || choice==0
            eval(['!','pkill fsleyes']);
            break;
        end
        eval(['!','pkill fsleyes'])
    end
end

eval(['!','pkill fsleyes'])
fprintf('complete'); disp(' ')
    
    
fprintf('Brain extraction complete'); disp(' ');

% extract 4D time series
fprintf('   Extracting 4D betted time series...');

% extract bold ts
bold_data = load_image([processdir,tmp_subj,'_BOLD_mcf.nii.gz']);
mask_data = single(load_image([processdir,tmp_subj,'_BOLD_meanBrain_mask.nii.gz']));
clear bold_ts
for jj = 1:size(bold_data,4)
    tmp_vol = bold_data(:,:,:,jj);
    bet_vol = tmp_vol.*mask_data;
    bold_ts(:,:,:,jj) = bet_vol;
end
save_image([processdir,tmp_subj,'_BOLD_mcf.nii.gz'],bold_ts,[processdir,tmp_subj,'_BOLD_ts.nii.gz'])

% extract asl ts
asl_data = load_image([processdir,tmp_subj,'_ASL_mcf.nii.gz']);
mask_data = single(load_image([processdir,tmp_subj,'_ASL_meanBrain_mask.nii.gz']));
clear asl_ts
for jj = 1:size(asl_data,4)
    tmp_vol = asl_data(:,:,:,jj);
    bet_vol = tmp_vol.*mask_data;
    asl_ts(:,:,:,jj) = bet_vol;
end
save_image([processdir,tmp_subj,'_ASL_mcf.nii.gz'],asl_ts,[processdir,tmp_subj,'_ASL_ts.nii.gz'])

fprintf('complete'); disp(' ');

%% =================== BOLD Signal Reconstruction ========================
disp('.........................................................')
disp(' '); disp(['Initiating BOLD signal reconstruction for ',tmp_subj]); disp(' ');

%% Window average BOLD data
% uses sliding window average - finds mean b/w control & label volumes

msg = sprintf('Window averaging BOLD ts...'); n=numel(msg); fprintf(msg)
bold_window_avg = window_average([processdir,tmp_subj,'_BOLD_ts.nii.gz'],...
    [processdir,tmp_subj,'_BOLD_window_avg.nii.gz']);
fprintf(repmat('\b',1,n)); fprintf([msg,'complete']); disp(' ');

figure('name','Window Averaging','numbertitle','off')
plot_image_ts([processdir,tmp_subj,'_BOLD_ts.nii.gz']); hold on;
plot_image_ts([processdir,tmp_subj,'_BOLD_window_avg.nii.gz']);
legend('raw','window averaged')

choice = menu('Assess Window Averaging. Continue?','Yes','No');
if choice==1 || choice==0
   close all
elseif choice==2
    close all
    disp('Window Averaging ERROR')
    return
end

%% Co-Registration and Transformation into 2mm MNI Space

disp('Transforming BOLD data into 2mm MNI space...');

% bet T1 anatomical image
fprintf('  Brain extracting structural image...'); disp(' ');
anatdir = [subj_path,tmp_subj,'/anat/'];
t1_img = dir([anatdir,'co*']).name;

eval(['!',fsl,'bet ',[anatdir,t1_img],' ',[processdir,tmp_subj,'_anat_bet.nii.gz']])

eval(['!',fsl,'fsleyes ',[anatdir,t1_img],' ',...
    [processdir,tmp_subj,'_anat_bet.nii.gz'],' --cmap blue-lightblue --alpha 50 &'])

choice = menu('Assess BET Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    defvalues = {'0','0.5'};
while(1)
prompt = {'Enter vertical gradient parameter (g; -1:1):','fractional intensity parameter (f; 0:1):'};
dlgtitle = 'Bet Input';
dims = [1,75];
definput = defvalues;
disp('f (0:1) - default=0.5; smaller values give larger brain outline estimates')
disp('g (-1:1) - default=0; positive values give larger brain outline at bottom, smaller at top')
BetInput = inputdlg(prompt,dlgtitle,dims,definput);
defvalues = {char(BetInput(1,1)),char(BetInput(2,1))};

fprintf('  Brain extracting structural image...')
anatdir = [subj_path,tmp_subj,'/anat/'];
t1_img = dir([anatdir,'co*']).name;
eval(['!',fsl,'bet ',[anatdir,t1_img],' ',[processdir,tmp_subj,'_anat_bet.nii.gz'],...
    ' -f ',char(BetInput(2,1)),' -g ',char(BetInput(1,1))])
fprintf('complete'); disp(' ');

eval(['!',fsl,'fsleyes ',[anatdir,t1_img],' ',[processdir,tmp_subj,'_anat_bet.nii.gz --cmap blue-lightblue --alpha 50 &']])

choice = menu('Assess BET Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
   break;
end

eval(['!','pkill fsleyes'])

end
end

%% fast segment structural image
fprintf('  Running FAST segmentation...')
eval(['!',fsl,'fast --segments ',[processdir,tmp_subj,'_anat_bet.nii.gz']])
fprintf('complete'); disp(' ');

% check segments
eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_anat_bet.nii.gz'],' ',[processdir,tmp_subj,'_anat_bet_seg_0.nii.gz'],...
    ' --cmap blue-lightblue --alpha 50 '...
    [processdir,tmp_subj,'_anat_bet_seg_1.nii.gz'],' --cmap green --alpha 50 ',...
    [processdir,tmp_subj,'_anat_bet_seg_2.nii.gz'],' --cmap red-yellow --alpha 50 ','&'])

choice = menu('Assess Segmentation Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    disp('FAST SEGMENTATION ERROR')
    return
end

%% co-registration of EPI data to structural image - use mean BOLD img for
% input
fprintf('  Co-registering BOLD EPI data to structural image...'); disp(' ');
eval(['!',fsl,'epi_reg --wmseg=',[processdir,tmp_subj,'_anat_bet_seg_2.nii.gz'],...
    ' --echospacing=0.00047 --pedir=-y --epi=',[processdir,tmp_subj,'_BOLD_meanBrain.nii.gz']...
    ' --t1=',[anatdir,t1_img],' --t1brain=',[processdir,tmp_subj,'_anat_bet.nii.gz'],...
    ' --out=',[processdir,tmp_subj,'_BOLD_epi2struct']])
fprintf('complete'); disp(' ');

if ~isfile([processdir,tmp_subj,'_BOLD_epi2struct.mat'])
    disp('ERROR: epi_reg file does not exist')
    return
end

% linear registration of T1 image to MNI space
fprintf('  Running FLIRT...');
eval(['!',fsl,'flirt -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain -in ',...
    [processdir,tmp_subj,'_anat_bet.nii.gz'],' -omat ',...
    [processdir,tmp_subj,'_T1_2_MNI.mat']])
fprintf('complete'); disp(' ');

if ~isfile([processdir,tmp_subj,'_T1_2_MNI.mat'])
    disp('ERROR: linear registration file does not exist')
    return
end

% non-linear registration of T1 image to MNI space
fprintf('  Running FNIRT...');
eval(['!',fsl,'fnirt --in=',[processdir,tmp_subj,'_anat_bet.nii.gz'],' --aff=',[processdir,tmp_subj,'_T1_2_MNI.mat'],...
    ' --cout=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],' --config=T1_2_MNI152_2mm'])
fprintf('complete'); disp(' ');
% edit: swapped [anatdir,t1_img] for [processdir,tmp_subj,'_anat_bet.nii.gz'] (12-12-2019)

if ~isfile([processdir,tmp_subj,'_T1_2_MNI_nonlinear.nii.gz'])
    disp('ERROR: nonlinear registration file does not exist')
    return
end

% apply transformation matricies to 4D func data
fprintf('  Applying transformation matricies to BOLD data...');
eval(['!',fsl,'applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --in=',...
    [processdir,tmp_subj,'_BOLD_window_avg.nii.gz'],' --warp=',...
    [processdir,tmp_subj,'_T1_2_MNI_nonlinear'],' --premat=',...
    [processdir,tmp_subj,'_BOLD_epi2struct.mat'],' --out=',...
    [processdir,tmp_subj,'_BOLD_MNI']])
fprintf('complete'); disp(' ');

% check registration in MNI space
eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_BOLD_MNI.nii.gz'],' $FSLDIR/data/standard/MNI152_T1_2mm_brain --cmap blue-lightblue --alpha 20 &'])
choice = menu('Assess Registration Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    disp('REGISTRATION ERROR')
    return
end


%% =================== ASL Signal Reconstruction ========================
disp('.........................................................')
disp(' '); disp(['Initiating ASL signal reconstruction for ',tmp_subj]); disp(' ');

%% Run surround subtraction on ASL ts
% difference between each image and the average of its 2 nearest
% neighbours is formed - reduces transient artifacts of BOLD weighting
% computes difference vols
 
image = [processdir,tmp_subj,'_ASL_ts.nii.gz'];
save_dir = processdir;
save_name = [tmp_subj,'_ASL_sub.nii.gz'];
asl_diff_vols = surround_subtraction(image,save_dir,save_name);

%% Compute voxelwise CBF0 map
%  index of vascular tension with partial volume and T2* corrections
%  (Chappell 2009/11)

% path to m0 image
m0dir = [subj_path,tmp_subj,'/m0/norm/'];
m0_file = dir([m0dir,'*nii.gz']).name;
m0_img = [m0dir,m0_file];

% run oxford_asl w/ partial volume  and T2* corrections - creates voxelwise
% cbf map in absolute units (mL/100g/min )
disp('Computing voxelwise CBF perfusion map')
fprintf('  Running Oxford_ASL...'); disp(' ');
eval(['!',fsl,'oxford_asl -i ',[processdir,tmp_subj,'_ASL_sub.nii.gz'],' -o ',...
       [processdir,tmp_subj,'_output_cbf'],' -m ',[processdir,tmp_subj,'_ASL_meanBrain_mask.nii.gz'],' --mc --iaf=diff ',...
       '--tis=2.665 --bolus=1.665 --casl -r ',[processdir,tmp_subj,'_anat_bet.nii.gz'],...
       ' --te=10 -c ',m0_img,' --cmethod=voxel ','--alpha=',num2str(0.84),' ',...
       ' --slicedt=0.0538 --t2star --echospacing=0.00047 --pvcorr --pvgm=',[processdir,tmp_subj,'_anat_bet_pve_1.nii.gz'],...
       ' --pvwm=',[processdir,tmp_subj,'_anat_bet_pve_2.nii.gz']])
fprintf('complete'); disp(' ');

% move absolute perfusion map into processdir
fprintf('  Generating ASL perfusion map...')
nativespacedir = [processdir,[tmp_subj,'_output_cbf'],'/native_space/'];
eval(['!','cp ',[nativespacedir,'perfusion_calib.nii.gz'],' ',processdir])
eval(['!','mv ',[processdir,'perfusion_calib.nii.gz'],' ',[processdir,tmp_subj,'_ASL_perfusion_map.nii.gz']])
fprintf('complete'); disp(' ');

%% check perfusion map

eval(['!',fsl,'fsleyes ',[processdir,tmp_subj,'_ASL_perfusion_map.nii.gz'],' &'])
choice = menu('Preview Perfusion Map. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    disp('Oxford_ASL ERROR')
    return
end

%% Resample ASL data into 2mm MNI space
% use linear and nonlinear matricies generated in previous steps for bold
% data

disp('Resampling ASL data into 2mm MNI space...')

% resample 4D ASL data into MNI
fprintf('  Resampling 4D ASL data...')
eval(['!',fsl,'applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --in=',...
    [processdir,tmp_subj,'_ASL_sub.nii.gz'],' --warp=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' --premat=',[processdir,tmp_subj,'_BOLD_epi2struct.mat'],' --out=',[processdir,tmp_subj,'_ASL_MNI']])
fprintf('complete'); disp(' ');

% resample perfusion map into MNI
fprintf('  Resampling ASL perfusion map...')
eval(['!',fsl,'applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --in=',...
    [processdir,tmp_subj,'_ASL_perfusion_map.nii.gz'],' --warp=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' --premat=',[processdir,tmp_subj,'_BOLD_epi2struct.mat'],' --out=',[processdir,tmp_subj,'_ASL_perfusion_map_MNI']])
fprintf('complete'); disp(' ');

%  resample segmented tissue masks into MNI
fprintf('  Resampling segmented tissue masks...')
eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[processdir,tmp_subj,'_anat_bet_seg_0.nii.gz'],' ','--warp=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[processdir,tmp_subj,'_anat_bet_seg_0_MNI.nii.gz']])
    
eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[processdir,tmp_subj,'_anat_bet_seg_1.nii.gz'],' ','--warp=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[processdir,tmp_subj,'_anat_bet_seg_1_MNI.nii.gz']])

eval(['!',fsl,'applywarp',' ','--ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain',' ',...
    '--in=',[processdir,tmp_subj,'_anat_bet_seg_2.nii.gz'],' ','--warp=',[processdir,tmp_subj,'_T1_2_MNI_nonlinear'],...
    ' ','--out=',[processdir,tmp_subj,'_anat_bet_seg_2_MNI.nii.gz']])
fprintf('complete'); disp(' ');

%% check resampled data in MNI space
eval(['!',fsl,'fsleyes $FSLDIR/data/standard/MNI152_T1_2mm_brain --cmap blue-lightblue --alpha 20 ',...
    [processdir,tmp_subj,'_ASL_MNI.nii.gz'],' ',[processdir,tmp_subj,'_BOLD_MNI.nii.gz',' &']])
choice = menu('Assess Registration Quality. Continue?','Yes','No');
if choice==1 || choice==0
   eval(['!','pkill fsleyes'])
elseif choice==2
    eval(['!','pkill fsleyes'])
    disp('REGISTRATION ERROR')
    return
end

disp('Resampling complete.')

%% =================== Signal Preprocessing ========================

disp('.........................................................'); disp(' '); 
disp(['Initiating Pre-processing for ',tmp_subj]); disp(' ');

%% Denoise BOLD data using temporal regression
% use fsl_regfilt to regress out nuisance parameters (i.e. cardiac/resp
% noise, non-neuronally related signal from motion, global signal).

disp('Denoising BOLD data via temporal regression...')

% extract WM signal
fprintf('  Extracting WM signal...')
image = [processdir,tmp_subj,'_BOLD_MNI.nii.gz'];
mask_image = [processdir,tmp_subj,'_anat_bet_seg_2_MNI.nii.gz'];
save_dir = processdir;
save_name = [tmp_subj,'_BOLD_WM_sig.nii.gz'];
[wm_sig,wm_ts] = apply_mask(image,mask_image,save_dir,save_name);
fprintf('complete'); disp(' ');

% extract CSF signal
fprintf('  Extracting CSF signal...')
image = [processdir,tmp_subj,'_BOLD_MNI.nii.gz'];
mask_image = [processdir,tmp_subj,'_anat_bet_seg_0_MNI.nii.gz'];
save_dir = processdir;
save_name = [tmp_subj,'_BOLD_CSF_sig.nii.gz'];
[csf_sig,csf_ts] = apply_mask(image,mask_image,save_dir,save_name);
fprintf('complete'); disp(' ');

% extract rigid-body motion parameters
fprintf('  Extracting rigid-body motion parameters...')
bold_motion = load([processdir,tmp_subj,'_BOLD_mcf.par']);
fprintf('complete'); disp(' ');

% create design matrix for fsl_regfilt
fprintf('  Creating GLM design matrix...')
bold_noise(:,1) = wm_ts;
bold_noise(:,2) = csf_ts;
bold_noise(:,3:8) = bold_motion;
dlmwrite('bold_regfilt_design.txt',bold_noise,'delimiter',' ');
fprintf('complete'); disp(' ');

% regress out nuisance parameters
fprintf('  Running fsl_regfilt...')
 eval(['!',fsl,'fsl_regfilt -i ',[processdir,tmp_subj,'_BOLD_MNI.nii.gz'],...
        ' -d bold_regfilt_design.txt -f "1,2,3,4,5,6,7,8" -o ',...
        [processdir,tmp_subj,'_BOLD_denoised.nii.gz']])
fprintf('complete'); disp(' ');

%% Denoise ASL data using temporal regression

disp('Denoising ASL data via temporal regression...')

% WM and CSF segments used as ROI for physiological noise, since neuronal
% activity in these regions is minimal

% extract WM signal
fprintf('  Extracting WM signal...')
image = [processdir,tmp_subj,'_ASL_MNI.nii.gz'];
mask_image = [processdir,tmp_subj,'_anat_bet_seg_2_MNI.nii.gz'];
save_dir = processdir;
save_name = [tmp_subj,'_ASL_WM_sig.nii.gz'];
[wm_sig,wm_ts] = apply_mask(image,mask_image,save_dir,save_name);
fprintf('complete'); disp(' ');

% extract CSF signal
fprintf('  Extracting CSF signal...')
image = [processdir,tmp_subj,'_ASL_MNI.nii.gz'];
mask_image = [processdir,tmp_subj,'_anat_bet_seg_0_MNI.nii.gz'];
save_dir = processdir;
save_name = [tmp_subj,'_ASL_CSF_sig.nii.gz'];
[csf_sig,csf_ts] = apply_mask(image,mask_image,save_dir,save_name);
fprintf('complete'); disp(' ');

% extract rigid-body motion parameters
fprintf('  Extracting rigid-body motion parameters...')
asl_motion = load([processdir,tmp_subj,'_ASL_mcf.par']);
fprintf('complete'); disp(' ');

% create design matrix for fsl_regfilt
fprintf('  Creating GLM design matrix...')
asl_noise(:,1) = wm_ts;
asl_noise(:,2) = csf_ts;
asl_noise(:,3:8) = asl_motion;
dlmwrite('asl_regfilt_design.txt',asl_noise,'delimiter',' ');
fprintf('complete'); disp(' ');

% regress out nuisance parameters
fprintf('  Running fsl_regfilt...')
 eval(['!',fsl,'fsl_regfilt -i ',[processdir,tmp_subj,'_ASL_MNI.nii.gz'],...
        ' -d asl_regfilt_design.txt -f "1,2,3,4,5,6,7,8" -o ',...
        [processdir,tmp_subj,'_ASL_denoised.nii.gz']])
fprintf('complete'); disp(' ');

%% Volume censoring (both bold and asl)
% remove any bold volumes with mean displacement >0.2mm - extract the same
% volumes from the asl data so the ts remain aligned.

disp('Running volume censoring (thresh=0.2mm)...')

% calculate mean displacement for bold data
bold_motion = load([processdir,tmp_subj,'_BOLD_mcf.par']);
fprintf('  Computing mean displacement...')
for aa = 2:length(bold_motion)
% calculate mean displacement
bold_mean_displacement(aa,:) = (((bold_motion(aa,1)-bold_motion(aa-1,1))^2)+...
((bold_motion(aa,2)-bold_motion(aa-1,2))^2)+((bold_motion(aa,3)-bold_motion(aa-1,3))^2)...
+((bold_motion(aa,4)-bold_motion(aa-1,4))^2)+((bold_motion(aa,5)-bold_motion(aa-1,5))^2)+...
((bold_motion(aa,6)-bold_motion(aa-1,6))^2))^(1/2);  
end
fprintf('complete'); disp(' ');

fprintf('  Volume censoring BOLD data...')
bold_data = load_image([processdir,tmp_subj,'_BOLD_denoised.nii.gz']);
censor_locs = find(bold_mean_displacement > 0.2);
counter = 1;
for tt = 1:size(bold_data,4)
    tmp_vol = bold_data(:,:,:,tt);
    censor = sum(tt == censor_locs);
    if censor == 0
       bold_censored_vols(:,:,:,counter) = tmp_vol;
       counter = counter + 1; 
    end
end
fprintf('complete'); disp(' ');

save_image([processdir,tmp_subj,'_BOLD_denoised.nii.gz'],bold_censored_vols,...
    [processdir,tmp_subj,'_BOLD_volcensored.nii.gz'])

fprintf('  Volume censoring ASL data...')
asl_data = load_image([processdir,tmp_subj,'_ASL_denoised.nii.gz']);
censor_locs = find(bold_mean_displacement > 0.2);
counter = 1;
for tt = 1:size(asl_data,4)
    tmp_vol = asl_data(:,:,:,tt);
    censor = sum(tt == censor_locs);
    if censor == 0
       asl_censored_vols(:,:,:,counter) = tmp_vol;
       counter = counter + 1; 
    end
end
fprintf('complete'); disp(' ');

save_image([processdir,tmp_subj,'_ASL_denoised.nii.gz'],asl_censored_vols,...
    [processdir,tmp_subj,'_ASL_volcensored.nii.gz'])

disp('Volume censoring complete.')

%% Smooth 4D data
% use 6mm full width at half-maximum (FWHM) Gaussian kernel.
disp('Smoothing 4D Data')

% BOLD data
fprintf('  Smoothing BOLD ts...')
smoothing_mm = 6;
constant_fsl = 2.3548; %FWHM (in mm) divide by 2.3548 to get sigma
sigma = smoothing_mm/constant_fsl;
eval(['!',fsl,'fslmaths ',[tmp_subj,'_BOLD_volcensored.nii.gz'],...
   ' -kernel gauss ',num2str(sigma),' -fmean ',[tmp_subj,'_BOLD_smoothed.nii.gz']])
fprintf('complete.'); disp(' ');

% ASL data
fprintf('  Smoothing ASL ts...')
smoothing_mm = 6;
constant_fsl = 2.3548; %FWHM (in mm) divide by 2.3548 to get sigma
sigma = smoothing_mm/constant_fsl;
eval(['!',fsl,'fslmaths ',[tmp_subj,'_ASL_volcensored.nii.gz'],...
   ' -kernel gauss ',num2str(sigma),' -fmean ',[tmp_subj,'_ASL_smoothed.nii.gz']])
fprintf('complete.'); disp(' ');

%% create mask in 2mm MNI
% required in next step for temporal filtering of 4D data
fprintf('Creating brain mask in 2mm MNI...');
file_path = [processdir,[tmp_subj,'_BOLD_MNI.nii.gz']];
save_dir = processdir;
save_name = [tmp_subj,'_MNI_2mm_mask.nii.gz'];
mask_data = create_mask(file_path,save_dir,save_name);
fprintf('complete.');  disp(' ');

%% Low pass filter bold data
% removes potential linear drift from scanner aquisition and perfusion
% weighting
% Neuronal signal , in general, will be below 0.15Hz - most of the higher frequency components are therefore noise.
% Filter at 1/2 nyquist frequency (1/4TR - Tak 2014) TR  = 4s

fprintf('Low-pass filtering BOLD signal (1/4*TR)...');
TR = 4; 
clear top bottom
bottom=0;
top = 1/(4*TR); % tak 2014

eval(['!',afni,'3dBandpass -mask ',[processdir,tmp_subj,'_MNI_2mm_mask.nii.gz'],...
     ' -overwrite -band ',num2str(bottom),' ',num2str(top),' -prefix ',...
     [processdir,tmp_subj,'_BOLD_lowpass'],' ',[processdir,tmp_subj,'_BOLD_smoothed.nii.gz']])
eval(['!',afni,'3dAFNItoNIFTI -overwrite -prefix ',[processdir,tmp_subj,'_BOLD_lowpass.nii.gz'],...
     ' ',[processdir,tmp_subj,'_BOLD_lowpass+*']])
eval(['!','rm ',[processdir,tmp_subj,'_BOLD_lowpass+*']])
fprintf('complete.'); disp(' ');

%% High pass filter ASL data
% removes T2* weighting from the bold contamination effect
% cuttoff frequency 1/4TR

fprintf('High-pass filtering ASL signal (1/4*TR)...');
TR = 4;
clear top bottom
bottom=1/(4*TR); % tak 2014
top = 999; % tak 2014    
eval(['!',afni,'3dBandpass -mask ',[processdir,tmp_subj,'_MNI_2mm_mask.nii.gz'],...
    ' -overwrite -band ',num2str(bottom),'  ',num2str(top),...
    ' -prefix ',[processdir,tmp_subj,'_ASL_highpass'],' ',[processdir,tmp_subj,'_ASL_smoothed.nii.gz']]) 
eval(['!',afni,'3dAFNItoNIFTI -overwrite -prefix ',[processdir,tmp_subj,'_ASL_highpass.nii.gz'],...
    ' ',[processdir,tmp_subj,'_ASL_highpass+*']])
eval(['!','rm ',[processdir,tmp_subj,'_ASL_highpass+*']])
fprintf('complete.'); disp(' ');

%% Demodulate high pass filtered asl data
% multiply each frame by cos[pi*n], where n denotes the volume number
% further attenuates BOLD contamination to accuratley estimate dynamic
% CBF changes

fprintf('Running demodulation for ASL data...');
tmp_img = load_image([processdir,tmp_subj,'_ASL_highpass.nii.gz']);
dyn = size(tmp_img,4);

 n = 0;
    for ii = 1:size(tmp_img,4)
        tmp_vol = tmp_img(:,:,:,ii);
        tmp_demod = arrayfun(@(x) x*cos(pi*ii),tmp_vol);
        demodulated_ts(:,:,:,ii) = tmp_demod;
        % progress update
        msg = sprintf('Processed %d/%d',ii,size(tmp_img,4));
        fprintf(repmat('\b',1,n));
        fprintf(msg);
        n=numel(msg);
    end

hdr_file = [processdir,tmp_subj,'_ASL_highpass.nii.gz'];
data  = demodulated_ts;
save_dir = [processdir,tmp_subj,'_ASL_demodulated.nii.gz'];
save_image(hdr_file,data,save_dir)
fprintf('...complete.'); disp(' ');

%% save pre-processed data

eval(['!','cp ',[processdir,tmp_subj,'_ASL_demodulated.nii.gz'],' ',[processdir,tmp_subj,'_ASL_preprocessed.nii.gz']])
eval(['!','cp ',[processdir,tmp_subj,'_BOLD_lowpass.nii.gz'],' ',[processdir,tmp_subj,'_BOLD_preprocessed.nii.gz']])

%% visualize preprocessed signals

close all
set(0,'DefaultFigureWindowStyle','docked');

plot1 = [processdir,tmp_subj,'_BOLD_preprocessed.nii.gz'];
label1 = char(extractAfter(plot1,'norm/'));
plot2 = [processdir,tmp_subj,'_ASL_preprocessed.nii.gz'];
label2 = char(extractAfter(plot2,'norm/'));

figure('name',label1,'numbertitle','off')
plot_image_ts(plot1)
title(label1,'interpreter','none')

figure('name',label2,'numbertitle','off')
plot_image_ts(plot2)
title(label2,'interpreter','none')

figure('name','combined','numbertitle','off')
plot_image_ts(plot1); hold on
plot_image_ts(plot2)
legend('BOLD','ASL');

choice = menu('Assess Preprocessed Signals. Continue?','Yes','No');
if choice==1 || choice==0
   close all
elseif choice==2
    close all
    disp('PRE-PROCESSING ERROR')
    return
end


end

disp(' '); disp('PREPROCESSING COMPLETE')
