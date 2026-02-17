%% 
% This program read individual autocorrelation maps and creates 
% across-subject autocorrelation map and its clustering output.
% Written by ali Golestani, September 2018
% Edited by Nichole Bouffard, October 2018 -Loops through subject list, averages across subjects. 
% Need to manually change R/L hem and size of N


clear
%% Specifying output directory
addpath('~/code/autocorrelation_scripts/AC_code');

homedir = '~/';

sub_list = fopen(['~/code/autocorrelation_scripts/lists' filesep 'subj_list.txt']);
subfile = textscan(sub_list, '%s');
Subj = subfile{1,1};

output_dir = [homedir 'AC_group_standard/'];

%hemisphere = {'L', 'R'};
hemisphere = {'R'};


for hem =1:length(hemisphere)
    
    for sbj = 1:length(Subj)    
                  
    %% Reading data
    data_name = strcat(homedir,'~/data/AC_maps_standard/', Subj{sbj},'/',Subj{sbj}, '_', hemisphere{hem}, '_', 'HPC_AC_Values.nii.gz');
    
    if exist(data_name) == 0
    continue;
    end
    data = MRIread(data_name);
    mm = data.vol;
    X = data.height;    
    Y = data.width;
    Z = data.depth;
    L = data.nframes;
    N = sum(sum(sum(mm(:,:,:,1)~=0)));
    XX = zeros(N,1);
    YY = zeros(N,1);
    ZZ = zeros(N,1);
    mm_mtrx = zeros(length(Subj),N,size(mm,4));
    
    %% Averaging across subjects
    Itr = 0;
    for x = 1:X
        for y = 1:Y
            for z = 1:Z
                if mm(x,y,z,1) ~= 0
                    Itr = Itr + 1;
                    XX(Itr) = x;
                    YY(Itr) = y;
                    ZZ(Itr) = z;
                    mm_mtrx(sbj,Itr,:) = mm(x,y,z,:);
                end
            end
        end
    end

    end %sbj
end%hem
 
% Store average across subjects    
Ave_sub(:,:) = nanmean(mm_mtrx,1);


%% Clustering the average
SimMat = -sqrt(squareform(pdist(Ave_sub)));
SimMat = (SimMat - min(min(SimMat)))/(max(max(SimMat)) - min(min(SimMat)));
W = SimMat - eye(size(SimMat)); %% Initial similarity matrix

%% Modulairty Maximization
[COMTY, ~] = cluster_jl(W,1,1);
idx = COMTY.COM{size(COMTY.COM,2)};
N_cluster = max(idx);
AC = zeros(1,N_cluster);
Ac = zeros(1,N_cluster);
Sz = zeros(1,N_cluster);
Index = zeros(1,N_cluster);

%% Calculating size and average autocorrelation of each cluster
Itr = 0;
for i = 1:max(idx)
    Sz(i) = sum(idx==i);
    if Sz(i) > 0 
        Itr = Itr + 1;
        ttmp = zeros(1,L);
        Ac_tmp = 0;
        Sm = 0;
        for j = 1:N
            if idx(j) == i 
                ttmp = ttmp + Ave_sub(j,:);
                Ac_tmp = Ac_tmp + Ave_sub(j,1);
                Sm = Sm + 1;
            end
        end
        Ref(1,:) = ttmp./Sm;
        AC(Itr) = sqrt(sum(Ref.^2));
        Ac(Itr) = Ac_tmp/Sm;
        Index(Itr) = i;
    end
end

%% Sorting clusters based on their autocorrelation
[~,IDX] = sort(Ac);
idx_updated = zeros(size(idx));
for i = 1:N_cluster
    idx_updated(idx == Index(IDX(i))) = i;
end

%% Creating output data
Clusters = zeros(X,Y,Z);
Map = zeros(X,Y,Z,L);
for i = 1:N
    Clusters(XX(i),YY(i),ZZ(i)) = idx_updated(i);
    Map(XX(i),YY(i),ZZ(i),:) = (Ave_sub(i,:));%-min(ac(:,1)))/(max(ac(:,1))-min(ac(:,1)));
end

%% Writing values into a text file
fileID_output = fopen(strcat(output_dir,'Cluster_Values_', hemisphere{hem},'.txt'),'w');
for i = 1:N_cluster
    Cluster_masked = (Clusters==i);
    Cluster_masked(isnan(Cluster_masked)) = 0;
    Mp(:,:,:) = Map(:,:,:,1);
    Mp_masked = Mp.*Cluster_masked;
    Mp_masked(isnan(Mp_masked)) = 0;
    if sum(sum(sum(Cluster_masked))) ~= 0
        Value = sum(sum(sum(Mp_masked)))/sum(sum(sum(Cluster_masked)));
        fprintf(fileID_output,'Cluster %i, AC Value = %f\n',i,Value);
    end
end
fclose(fileID_output);
   
dataOut = data;
         
dataOut.vol = Clusters;
PTH = strcat(output_dir,'Group_', hemisphere{hem}, '_Clusters.nii.gz');
MRIwrite(dataOut,PTH); 

dataOut.vol = Map;
PTH = strcat(output_dir,'Group_', hemisphere{hem}, '_AC_Values.nii.gz');
MRIwrite(dataOut,PTH); 


