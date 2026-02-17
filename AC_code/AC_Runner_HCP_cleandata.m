clear;
%% Specifying input directory
addpath('~/code/autocorrelation_scripts/AC_code');

homedir = '~/';
hcpdir= '/home/nichole/s3-drive/HCP_1200/'; % directly mounting the HCP project repository


sub_list = fopen(['~/code/autocorrelation_scripts/lists' filesep 'subj_list.txt']);
subfile = textscan(sub_list, '%s');
sub = subfile{1,1};

% run_list = fopen([home filesep 'lists' filesep 'run_list.txt']);
% runfile = textscan(run_list, '%s');
% run = runfile{1,1};

%% Reading pre-processed fMRI data

for i=1:length(sub)
    % create output folder
    
clear data;
clear func;

    output_dir = [homedir '/data/AC_maps' filesep sub{i}];

    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
    
func_name = ([hcpdir sub{i} filesep 'MNINonLinear/Results/rfMRI_REST1_7T_PA' filesep 'rfMRI_REST1_7T_PA_hp2000_clean.nii.gz']);   
data = MRIread(func_name);
func = data.vol;
 
%TR = .001*data.tr;
TR = 1;  %Specify TR if the header is wrong

%% Initial parameters
L = 5; %% Length of autocorrelation vector

%% Right Hippocampus mask
ROI_mask_name = ([homedir '~/code/autocorrelation_scripts/HippUnfold_subfields_from_jordan/func_bin_masks/thresholded_0.4' filesep 'sub-HCY' sub{i} '_R_wholeHipp_func.nii.gz' ]);
data = MRIread(ROI_mask_name);
ROI_mask = data.vol;

%% Running AC estimation and clustering function
[Map,Clusters] = AC_generator_HCP(func,ROI_mask,L,TR,1);

%% Writing output images
dataOut = data;
dataOut.vol = Clusters;
PTH = ([output_dir filesep sub{i} '_R_HPC_Clusters_regthresh0.4_cleandata.nii.gz' ]);
MRIwrite(dataOut,PTH); 

dataOut.vol = Map;
PTH = ([output_dir filesep sub{i} '_R_HPC_AC_Values_regthresh0.4_cleandata.nii.gz' ]);
MRIwrite(dataOut,PTH); 


%% Left Hippocampus mask
ROI_mask_name = ([homedir '~/code/autocorrelation_scripts/HippUnfold_subfields_from_jordan/func_bin_masks/thresholded_0.4' filesep 'sub-HCY' sub{i} '_L_wholeHipp_func.nii.gz' ]);
data = MRIread(ROI_mask_name);
ROI_mask = data.vol;

%% Running AC estimation and clustering function
[Map,Clusters] = AC_generator_HCP(func,ROI_mask,L,TR,1);

%% Writing output images
dataOut = data;
dataOut.vol = Clusters;
PTH = ([output_dir filesep sub{i} '_L_HPC_Clusters_regthresh0.4_cleandata.nii.gz' ]);
MRIwrite(dataOut,PTH); 

dataOut.vol = Map;
PTH = ([output_dir filesep sub{i} '_L_HPC_AC_Values_regthresh0.4_cleandata.nii.gz' ]);
MRIwrite(dataOut,PTH); 

fprintf('\n Done subject %s \n', sub{i});

clear func;
clear data;
end
