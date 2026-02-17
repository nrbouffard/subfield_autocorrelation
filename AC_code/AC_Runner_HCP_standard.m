clear;
%% Specifying input directory
addpath('~/code/autocorrelation_scripts/AC_code');

homedir = '~/';
hcpdir= '/home/nichole/s3-drive/HCP_1200/'; % mounting the HCP project repository directly


sub_list = fopen(['~/code/autocorrelation_scripts/lists' filesep 'subj_list.txt']);
subfile = textscan(sub_list, '%s');
sub = subfile{1,1};


%% Reading pre-processed fMRI data

for i=1:length(sub)
    % create output folder
    
    output_dir = [homedir 'AC_maps_standard' filesep sub{i}];

    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
func_name = ([homedir '/data/struct_to_func_registration/' sub{i} filesep 'registration.nii.gz']);
func = data.vol;
 
%TR = .001*data.tr;
TR = 1;  %Specify TR if the header is wrong

%% Initial parameters
L = 4; %% Length of autocorrelation vector

%% Right Hippocampus mask
ROI_mask_name = ([homedir '~/code/autocorrelation_scripts/full_hipp_mask/BNA_R_full_hippo_bin.nii.gz' ]);
data = MRIread(ROI_mask_name);
ROI_mask = data.vol;

%% Running AC estimation and clustering function
[Map,Clusters] = AC_generator_HCP(func,ROI_mask,L,TR,1);

%% Writing output images
dataOut = data;
dataOut.vol = Clusters;
PTH = ([output_dir filesep sub{i} '_R_HPC_Clusters.nii.gz' ]);
MRIwrite(dataOut,PTH); 

dataOut.vol = Map;
PTH = ([output_dir filesep sub{i} '_R_HPC_AC_Values.nii.gz' ]);
MRIwrite(dataOut,PTH); 


%% Left Hippocampus mask
ROI_mask_name = ([homedir '~/code/autocorrelation_scripts/BNA_L_full_hippo_bin.nii.gz']);
data = MRIread(ROI_mask_name);
ROI_mask = data.vol;

%% Running AC estimation and clustering function
[Map,Clusters] = AC_generator_HCP(func,ROI_mask,L,TR,1);

%% Writing output images
dataOut = data;
dataOut.vol = Clusters;
PTH = ([output_dir filesep sub{i} '_L_HPC_Clusters.nii.gz' ]);
MRIwrite(dataOut,PTH); 

dataOut.vol = Map;
PTH = ([output_dir filesep sub{i} '_L_HPC_AC_Values.nii.gz' ]);
MRIwrite(dataOut,PTH); 

fprintf('\n Done subject %s \n', sub{i});

clear func;
clear data;
end
