function [eval, EVect, EVect2, Components, S, Ddiff] = dm_dmit(X, varargin)
%Alexey Ryabov, 2021
%Find diffusion map for the imput matrix X
% [eval, EVect, EVect2, Components, S, Ddiff] = dm_dmit(S) returns nonzero
% eigenvalues eval, corresponding eigenvectors as columns of matrix EVect,
% rescaled eigenvectors EVect2 = EVect/eval, which can be used to define
% traits, the number of Components in the diffusion map (the number of zero
% eigenvalues), similarity matrix S and the matrix of diffusion distances 
% in the space defined by EVect2 
%
%
% [...] = dm_dmit(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%     parameter name/value pairs to control the computation and handling of
%     special data types. Parameters are:
%
% 'Metric' -   Metric of the similarity matrix. Choices are:
%             'euclidean'; 'cityblock'; 'chebychev'; 'mahalanobis'; 'minkowski';
%             'cosine';  'spearman'; 'pearson'; 'kendall'; 'jaccard';'gaussian'
%             see pdist for their meaning. 
%             When metric is distance based, e.g. 'euclidean' 'cityblock' 
%             then similarity matrix S = 1/D, 
%             for correlation based metric  ('spearman'; 'pearson'; 'kendall') similarity S = (D + 1)/2,
%             for cosine, jaccard and gaussian based distance the similarity S = 1-D
% 'Norm' -    Normalization of input varaibles. Choices are:
%     'none' no normalization
%     'normalized'  Normilized along columns, so that sum(X.^2, 1) = 1
%     'zscore' each column standartized to have mean = 0 and std = 1, X = (X-mu)/sigma
% 'k_min' -   the minimal number of entries in the cleaned similarity
%             matrix, default k_min = 10
%     'Laplacian' -   Which laplacian matrix should be used. Choices are:
%        'rownorm' - row normalized Laplacian L = -Sij/sum(Sij) (i \ne j), Lii = 1;
%        'Lafon' - %from Lafon's presentation, L = Sij/sum(Sij)
%
%   Example:
%   [ev, aEV, aSClnd, indPos, D] = dm_dmit(X);
%   [ev, aEV, aSClnd, indPos, D] = dm_dmit(X, 'Laplacian', 'Lafon');
%
%  See also dm_dmsim, dm_plot, dm_simmat
% References:
%   [1]

%Get similarity matrix
S = dm_simmat(X, varargin{:});

if nargout >= 6 % we need to find the diffusion distance in EVect2 space
    [eval, EVect, EVect2, Components, Ddiff] = dm_dmsim(S,  varargin{:});
else
    [eval, EVect, EVect2, Components] = dm_dmsim(S, varargin{:});
end



end
