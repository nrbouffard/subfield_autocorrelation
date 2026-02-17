clear 


sub_list = fopen(['~/code/autocorrelation_scripts/lists' filesep 'subj_list.txt']);
subfile = textscan(sub_list, '%s');
SUBJ = subfile{1,1};

%SUBJ = {'126426', '130114','130518'};
% 
ROIs = {'L','R'};

Masks = {'CA1','CA4DG', 'CA2CA3','subiculum','SRLM'};

for sbj = 1:length(SUBJ)

    for roi = 1:length(ROIs)

    for msk = 1:length(Masks)
        mm_name = strcat('~/code/autocorrelation_scripts/HippUnfold_subfields_from_jordan/separate_bin_subfield_masks/', SUBJ{sbj}, '/func_masks/sub', SUBJ{sbj}, '_', ROIs{roi}, '_', Masks{msk}, '_func_thresh0.4.nii.gz');
        data = MRIread(mm_name);
        mask = data.vol;
           
 
            mm_name = strcat('~/data/AC_maps/',SUBJ{sbj}, '/', SUBJ{sbj}, '_', ROIs{roi}, '_HPC_Clusters_regthresh0.4_cleandata.nii.gz');
            data = MRIread(mm_name);
            gCluster = data.vol;
            
            gCluster = gCluster.*mask;
            
            N_Voxel = sum(sum(sum(gCluster~=0)));
            XX = zeros(N_Voxel,1);
            YY = zeros(N_Voxel,1);
            ZZ = zeros(N_Voxel,1);

            X = data.height;    
            Y = data.width;
            Z = data.depth;
            Itr = 0;
            for x = 1:X
                Im(:,:) = gCluster(x,:,:);
                if sum(sum(Im))~=0
                    Itr = Itr + 1;
                    Ln(Itr) = sum(sum(Im));
                    Nm(Itr) = sum(sum(Im~=0));
                end
            end
            Ln_X_gCluster = Ln./Nm;
            LX = length(Ln_X_gCluster);

            clear Im Ln Nm;
            Itr = 0;
            for y = 1:Y
                Im(:,:) = gCluster(:,y,:);
                if sum(sum(Im))~=0
                    Itr = Itr + 1;
                    Ln(Itr) = sum(sum(Im));
                    Nm(Itr) = sum(sum(Im~=0));
                end
            end
            Ln_Y_gCluster = Ln./Nm;
            LY = length(Ln_Y_gCluster);

            clear Im Ln Nm;
            Itr = 0;
            for z = 1:Z
                Im(:,:) = gCluster(:,:,z);
                if sum(sum(Im))~=0
                    Itr = Itr + 1;
                    Ln(Itr) = sum(sum(Im));
                    Nm(Itr) = sum(sum(Im~=0));
                end
            end
            Ln_Z_gCluster = Ln./Nm;
            LZ = length(Ln_Z_gCluster);

            clear Im Ln Nm;
      
                
                file_name = strcat('~/data/AC_maps/', SUBJ{sbj}, '/', SUBJ{sbj}, '_', ROIs{roi}, '_HPC_AC_Values_regthresh0.4_cleandata.nii.gz');

                %% Specifying output directory
                data = MRIread(file_name);
                Values = data.vol;
                
                Values = Values.*mask;

                Itr = 0;
                for x = 1:X
                    Im(:,:) = Values(x,:,:,1);
                    if sum(sum(Im))~=0
                        Itr = Itr + 1;
                        Ln(Itr) = sum(sum(Im));
                        Nm(Itr) = sum(sum(Im~=0));
                    end
                end
                Ln(Itr+1:LX) = 0;
                Nm(Itr+1:LX) = 1;
                Ln_X(1,:) = Ln./Nm;

               % Ln_X(sbj,:) = (Ln_X2(sbj,:)-mean(Ln_X2(sbj,:)))/std(Ln_X2(sbj,:));

                clear Im Ln Nm;
                Itr = 0;
                for y = 1:Y
                    Im(:,:) = Values(:,y,:,1);
                    if sum(sum(Im))~=0
                        Itr = Itr + 1;
                        Ln(Itr) = sum(sum(Im));
                        Nm(Itr) = sum(sum(Im~=0));
                    end
                end
                Ln(Itr+1:LY) = 0;
                Nm(Itr+1:LY) = 1;
               % if sbj>2
               %     if length(Ln)<size(Ln_Y2,2)
              %          Ln = [Ln Ln(end)];
              %          Nm = [Nm Nm(end)];
              %      end
              %  end
                Ln_Y(1,:) = Ln./Nm;
             %   Ln_Y2(sbj,:) = (Ln_Y2(sbj,:)-mean(Ln_Y2(sbj,:)))/std(Ln_Y2(sbj,:));

                clear Im Ln Nm;
                Itr = 0;
                for z = 1:Z
                    Im(:,:) = Values(:,:,z,1);
                    if sum(sum(Im))~=0
                        Itr = Itr + 1;
                        Ln(Itr) = sum(sum(Im));
                        Nm(Itr) = sum(sum(Im~=0));
                    end
                end
                Ln(Itr+1:LZ) = 0;
                Nm(Itr+1:LZ) = 1;
                Ln_Z(1,:) = Ln./Nm;
                %Ln_Z2(sbj,:) = (Ln_Z2(sbj,:)-mean(Ln_Z2(sbj,:)))/std(Ln_Z2(sbj,:));

                clear Im Ln Nm;
         
             subnum = str2double(SUBJ{sbj});
             
            save(strcat('/output/', SUBJ{sbj},'_Lines_', ROIs{roi},'_HPC_', Masks{msk}),'subnum','Ln_X','Ln_Y','Ln_Z','Ln_X_gCluster','Ln_Y_gCluster','Ln_Z_gCluster')
            clear Ln_X Ln_Y Ln_Z;
    end
      
    end

end
