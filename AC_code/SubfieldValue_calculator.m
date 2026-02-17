function Values = ClusterValue_calculator(AC_Map,clusters)
%% This functuion returns the average AC values for each ROI
%% Inputs:
%%      AC_Map: 4D map of AC values at each lag
%%      Clusters: 3D mask of the Head/Body/Tail HPC ROIs
%% Outputs:
%%      Values: 2D matrix of the averaged AC values in each ROI at each lag 
%       Col 1 = Lag 1, Col 2 = Lag 2...
%       Row 1 = Cluster 3, Row 2 = Cluster 2, Row 3 = Cluster 1

X = size(AC_Map,1);
Y = size(AC_Map,2);
Z = size(AC_Map,3);
L = size(AC_Map,4);

%% Reading AC and calculating average AC
N_clusters = max(clusters(:));
Values = zeros(N_clusters,L);
N_Val = zeros(N_clusters,L);
for x = 1:X
    for y = 1:Y
        for z = 1:Z
            for n = 1:N_clusters
                for lag = 1:L
                    if clusters(x,y,z)== n
                        Values(n,lag) = Values(n,lag) + AC_Map(x,y,z,lag);
                        N_Val(n,lag) = N_Val(n,lag) + 1;
                    end
                end
            end
        end
    end
end
Values = Values./N_Val;

