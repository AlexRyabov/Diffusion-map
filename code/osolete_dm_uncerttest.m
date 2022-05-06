function  [RelError] = f_DM_SimilartyTest(UnstackedData, DistanceMetrics, k_max, M, GroupByID_idx)
%Alexey Ryabov 2021
%find mean relative error between the similarity matrix obtained with entire set of samples
%and similarity martix obtained for M bootstraped samples
%if M is a vector, then the errorr will be calculated for each Mi 
% UnstackedData has features in rows 
% DistanceMetrics what distance metrics to use
% k_max is the maximal number of nonzero components in the similarity
% matrix
% if GroupByID is defined, the samples will be grouped by this index, and
% bootstrabed will be done using unique values of this index
% size(SpUnStacked, 1) should equal length(GroupByID_idx)

if ~exist('GroupByID_idx', 'var')
   GroupByID_idx = 1:size(UnstackedData, 2);
end
GroupByIDUnique_idx = unique(GroupByID_idx);
%Get similarity matrixes for all data
aSFulls = {};
for iDM = 1:length(DistanceMetrics)
    DistanceMetric = DistanceMetrics{iDM};
    %[ev, aEV, aSFull, indPosFull, DiffDist] = f_DM_DMit(UnstackedData, DistanceMetric, false, k_max);
    [ev, aEV, aSFull, indPosFull, DiffDist] = f_DM_DMit(UnstackedData, DistanceMetric, false, k_max);
    aSFulls{iDM} = aSFull;
end

% calculate similarity matrixes for bootstraped samples
m_max = length(GroupByIDUnique_idx);
RelError = NaN(length(M), length(DistanceMetrics));
for im  = 1: length(M)
    for iDM = 1:length(DistanceMetrics)
        % get m randomly chosen samples from GroupByID
        RandInd = randi(m_max, M(im), 1);
        GroupByID_idx_sel = [];
        for ir = 1:length(RandInd)
            GroupByID_idx_sel = [GroupByID_idx_sel,  find(GroupByID_idx == GroupByIDUnique_idx(RandInd(ir)))];
        end
        % for these samples select from SpUnStacked usigng GroupByID_idx
        SpUnStacked_sel = UnstackedData(:, GroupByID_idx_sel);
        DistanceMetric = DistanceMetrics{iDM};
        [ev, aEV, aS, indPos, DiffDist] = f_DM_DMit(SpUnStacked_sel, DistanceMetric, false, k_max);
        aSFull = aSFulls{iDM};
        aS(isnan(aS))=0;
        %find relative deviation
        %Delta = (aSFull - aS)./(0.5*(aSFull + aS));
        %Delta = (aSFull - aS);
        
        %Delta(aSFull==0 & aS==0) = 0;
        %Delta(aSFull==0 & aS>0) = 1;
        %RelError(im, iDM) = (nanmean(abs(Delta(Delta~=0))));
        %%alterntive 1
        %Delta = (aSFull - aS)./max(aSFull,aS);
        %Delta(aSFull==0 & aS==0) = 0;
        %alternative 2
        %RelError(im, iDM) = (nanstd((Delta(Delta~=0))));
        RelError(im, iDM) = corr(aSFull(:), aS(:));
    end
end
end