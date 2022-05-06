function [Nspan, SUncert] = dm_learncrv(X, varargin)
%Alexey Ryabov, 2021
%Find learning curve for given data
% [Nspan, Scov, Lcov] = dm_learncrv(X) returns learning curve for input
% matrix X
%
%
% [...] = dm_learncrv(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%     parameter name/value pairs to control the computation and handling of
%     special data types. Parameters are:
%     contains the same paramters as dm_dmit,
% 'Metric' -   Metric of the similarity matrix or a cell array of metrics
% 'Nspan'  array with the number of smaples used for calculating the
% similarity, or integer, then this is the number of elements in this array
% Nspan (logdistributed)

%   Example:

%   [Nspan, Scov, Lcov] = dm_learncrv(S);
%   [Nspan, Scov, Lcov] = dm_learncrv(S, 'Nspan', [5, 10, 20, 40, 80]);  %find
%   correlations for 5, 10, 20, 40 and 80 samples
%   [Nspan, Scov, Lcov] = dm_learncrv(S, 'Bootstrap', true);  %use
%   bootstraping of samples (all rows can be used )
%   [Nspan, Scov, Lcov] = dm_learncrv(S, 'Metric', {'euclidean'; 'cosine';  'spearman'; 'gaussian'});
%
%  See also dm_dmit
% References:
%   [1]


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% Validate the matrix X
[XRows, XCols] = size(X);
if XRows == 0 || XCols == 0
    error('dm_learncrv:Input matrix X should not be empty');
end
internal.stats.checkSupportedNumeric('X',X); % not complex


% Parse arguments and check if parameter/value pairs are valid
%Norm  none/NORMALIZED/standardized
%Distance — Distance metric

paramNames = {'Metric',   'Nspan', 'Ndots', 'Bootstrap'};
defaults   = {{[]},  [], 20, false};

[Metrics, Nspan, Ndots, Bootstrap, ~, PassFurther]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});


% Validate Nspan and construct the numbers of samples, log distributed
% array

if isempty(Nspan) 
    if Bootstrap 
        Nspan = unique(round(logspace(log10(50), log10(XRows), Ndots)));
    else
        Nspan = unique(round(logspace(log10(50), log10(XRows/2), Ndots)));
    end
else
    %Nspan is vector
    if Bootstrap
        if  max(Nspan) > XRows
            error('dm_learncrv: The maximum value of in Nspan (%i) must be less than or equal to the number of rows in matrix X if Bootstrap = true', max(Nspan))
        end
    else
        if  max(Nspan) > XRows/2
            error('dm_learncrv: The maximum value of in Nspan (%i) must be less than or equal to half the number of rows in matrix X if Bootstrap = false', max(Nspan))
        end
    end
end





SUncert = NaN(length(Nspan), length(Metrics));
for im  = 1: length(Nspan)
    for iDM = 1:length(Metrics)
        Metric = Metrics{iDM};
        if Bootstrap
            % get m randomly chosen samples from GroupByID
            sampl1 =  randsample(XRows, Nspan(im), true);
            sampl2 =  randsample(XRows, Nspan(im), true);
        else
            %get two nonoverlapping samples
            sampl12 =  randsample(XRows, Nspan(im)*2);
            sampl1 = sampl12(1:Nspan(im));
            sampl2 = sampl12(Nspan(im)+1:end);
        end
        if isempty(Metric)
            [~, ~, ~, ~, S1] = dm_dmit(X(sampl1, :), PassFurther{:});
            [~, ~, ~, ~, S2] = dm_dmit(X(sampl2, :), PassFurther{:});
        else
            [~, ~, ~, ~, S1] = dm_dmit(X(sampl1, :), 'Metric', Metric, PassFurther{:});
            [~, ~, ~, ~, S2] = dm_dmit(X(sampl2, :), 'Metric', Metric, PassFurther{:});
        end
        S1(isnan(S1))=0;
        S2(isnan(S2))=0;
        SUncert(im, iDM) = corr(S1(:), S2(:));
    end
end
end




