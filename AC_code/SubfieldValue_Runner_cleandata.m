clear;
%% Specifying input directory
addpath('~/code/autocorrelation_scripts/AC_code');

homedir = '~/';

sub_list = fopen(['~/code/autocorrelation_scripts/lists' filesep 'subj_list.txt']);
subfile = textscan(sub_list, '%s');
sub = subfile{1,1};

 hem = {'L','R'};
 roi = {'CA1','CA2CA3','CA4DG','subiculum', 'SRLM'};


%% Reading autocorrelation value maps

for i=1:length(sub)
    % create output folder
    
    output_dir = [homedir filesep 'subfield_AC_Values_cleandata' filesep sub{i}];
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
        
    for h=1:length(hem)
        
        for s=1:length(roi)

% Registration threshold 0.4   
func_name = ([homedir '/data/AC_maps' filesep sub{i} filesep sub{i} '_' hem{h} '_HPC_AC_Values_regthresh0.4_cleandata.nii.gz']);

    
    if ~exist(func_name, 'file')
        continues
    end

data = MRIread(func_name);
AC_Map = data.vol;
    
%% AC hippocampus cluster mask
% Registration threshold 0.4
ROI_mask_name = ([homedir '~/code/autocorrelation_scripts/HippUnfold_subfields_from_jordan/separate_bin_subfield_masks' filesep sub{i} filesep 'func_masks' filesep 'sub' sub{i} '_' hem{h} '_' roi{s} '_func_thresh0.4.nii.gz' ]);

data = MRIread(ROI_mask_name);
subfields = data.vol;

%% Running AC estimation and clustering function
Values = SubfieldValue_calculator(AC_Map,subfields);

%% Writing output 
xlswrite(strcat([output_dir filesep sub{i} '_' hem{h} '_' roi{s} '_regthresh0.4_AC_Values_cleandata.xlsx']),Values);


fprintf('\n Done subject %s run %s \n', sub{i})

         end % subfield
     end % hem
end % sub