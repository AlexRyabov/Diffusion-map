function [S, Sfull] = dm_simmat(X, varargin)
%Alexey Ryabov, 2021
%Find similarity matrix beween columns (variables) of data matrix X and select 10 largest entries in each column. 
%Rows of X correspond to observations (e.g. dates or stations) 
%and columns to variables (e.g. species).
%
% [S] = dm_simty(X) returns the cleaned matrix S where
% each column contains k_max largest entries, by default k_max = 10
% for (n * p) matrix X. Rows of X correspond to observations and columns to
% variables.
%
% [S, Sfull] = dm_simmat(X) also returns the uncleaned similarity matrix Sfull
%
% [...] = dm_simmat(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%     parameter name/value pairs to control the computation and handling of
%     special data types. Parameters are:
%
% 'Metric' -   Metric of the similarity matrix. Choices are:
%             'euclidean'; 'cityblock'; 'chebychev'; 'mahalanobis'; 'minkowski';
%             'cosine';  'spearman'; 'pearson'; 'kendall'; 'jaccard';'gaussian'
%             see pdist for their meaning. 
%
%             When metric is distance based, e.g. 'euclidean' 'cityblock' 
%             then similarity matrix S = 1/D, 
%
%             for correlation based metric  ('spearman'; 'pearson'; 'kendall') similarity S = (D + 1)/2,
%             for correlation metrics NaN will be removed, and rho(i,j) 
%             will be computed using rows with no missing values in column
%             i or j (see corr(X, 'Row', 'pairwise') )
%             
%             for cosine, jaccard and gaussian based distance the similarity S = 1-D
% 'Norm' -    Normalization of input varaibles when distance metric is used. Choices are:
%     'none' no normalization
%     'normalized'  Normilized along columns, so that sum(X.^2, 1) = 1
%     'zscore' each column standartized to have mean = 0 and std = 1, 
%     X = (X-mu)/sigma ,
%     Note, normalization is not performed when correlation metrics are
%     used, to avoid problems with nomalization of NaNs, which are not used
%     for calculating correlations

% 'k_min' -   : the minimal number of entries in the cleaned similarity
% matrix, default k_min = 10
%
%   Example:
%  [Sclned, Sfull] = dm_simmat(X);
%   [Sclned] = dm_simmat(r, 'k_min', 10,  'Metric', 'cosine', 'Norm', 'normalized');
%
%  See also dm_dmit, dm_plot, pdist, zscore
% References:
%   [1]

%get paramters
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[n, p] = size(X);
if n == 0 || p == 0
    error('Input matrix X should not be empty');
end
internal.stats.checkSupportedNumeric('X',X); % not complex

% Parse arguments and check if parameter/value pairs are valid
%Norm  none/NORMALIZED/standardized
%Distance — Distance metric


paramNames = {'Metric',   'Norm', 'k_min'};
defaults   = {'spearman', 'normalized', 10};

[Metric, Normalization, k_min, sf, rest]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

% Validate String value for  Distance value
MetricNames = {'euclidean'; 'cityblock'; 'chebychev'; ...
    'mahalanobis'; 'minkowski'; 'cosine';  ...
    'spearman'; 'pearson'; 'kendall'; 'jaccard';'gaussian'};
Metric = internal.stats.getParamVal(Metric,MetricNames,...
    '''Similarity metric''');

% Validate String value for Normalization value
Normalizations = {'none'; 'normalized'; 'zscore'};
Normalization = internal.stats.getParamVal(Normalization,Normalizations,...
    '''Norm''');

% Validate integer value for 'k_min' option
if k_min > p
    warning('''k_min (%i)'' should be less or equal the number of columns in X', k_min)
end

% Validate the number of components option 'NumComponents'
if isempty(X)
    error(message('diffusion map:Input matrix X is empty'));
end


% indPos = find(nansum(X, 1)>0);
% XInit = X;
% 
% tPlAbndUns = (:, indPos); %take only positevely defined




%Find the similarity matrix
switch Metric
    case {'cosine', 'jaccard'}
        X = Normalize(X, Normalization);
        D = squareform(pdist(X',Metric));
        D(isnan(D)) = 1;  %Set distance to 1 when it is undefined
        Sfull = 1-D;
    case 'gaussian'  %Gaussian distance
        X = Normalize(X, Normalization);
        D = squareform(pdist(X','euclidean'));
        Sfull = exp(-D.^2);
    case {'pearson', 'spearman', 'kendall'}  %correlation distance.
        Sfull = (corr(X,'Type', Metric, 'Row', 'pairwise')+1)/2;
        Sfull(isnan(Sfull)) = 0;  %Set similarity to 0 when it is undefined
        D = 1-Sfull;
    otherwise %'euclidean', 'cityblock' ... 
        X = Normalize(X, Normalization);
        D = squareform(pdist(X',Metric));
        Sfull = 1./D;
end
%Cleaned similarity matrix
S = Sfull;

%set diagonal elements to zero
for iS = 1:size(S, 1)
    S(iS, iS) = 0;
end

%Remove infinite similatiry entiries
ind_inf = S == Inf;
if sum(ind_inf(:)) > 0
    [i_inf, j_inf] = find(ind_inf == 1);
    warning('The following pairs of variables have infinite similarity. We recommend either combining these variables or using a correlation, cosine, or gaussian similarity')
    display([i_inf, j_inf]);
    MaxSim = max(S(S(:) < Inf));
    S(S == Inf) = MaxSim;
    warning('Infinite similarity is set to max observed similarity of %f', MaxSim)
end

%clean the similarity matrix (leave k_min largest elements in each colum
B =(maxk(S,k_min));
S(S<repmat(B(end,:), size(S, 1), 1)) = 0;
S = max(S, S');

% %Check if there are variables with zero similarities to all other variables
% SimSum = sum(S, 1);
% ind = SimSum==0;
% if sum(ind) > 0
%     %get indexes of 0 related variables
%     Ind = 1:size(S, 2);
%     Ind = Ind(ind);
%     %show warning
%     if length(Ind) <= 10
%         warning('variables %s have zero similarity to other variables. \n These variables should be excluded from the diffusion map analysis.', num2str(Ind));
%     else
%         warning('variables %s ... (more than 10) have zero similarity to other variables. \n These variables should be excluded from the diffusion map analysis.', num2str(Ind(1:10)));
%     end
% end

end


function X = Normalize(X, Normalization)
%Normalize
switch Normalization
    case 'normalized'
        fNorm = sqrt(sum((X).^2, 1));
        fNorm(fNorm==0) = 1;
        X = X./ fNorm;
    case 'zscore'
        X = zscore(X);
    otherwise
        %nothing to do
end
end