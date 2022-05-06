function [eval, EVect, EVect2, Components, Ddiff] = dm_dmsim(S, varargin)
%Alexey Ryabov, 2021
%Diffusion map when similarity matrix is know
% S is the similarity matrix obtained with dm_simmat,
% [eval, EVect, EVect2, Components, Ddiff] = dm_dmsim(S) returns nonzero
% eigenvalues eval, corresponding eigenvectors as columns of matrix EVect,
% rescaled eigenvectors EVect2 = EVect/eval, which can be used to define
% traits, the number of Components in the diffusion map (the number of zero
% eigenvalues) and the matrix of diffusion distances in the space defined
% by EVect2
%
%
% [...] = dm_dmsim(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%     parameter name/value pairs to control the computation and handling of
%     special data types. Parameters are:
%
%     'Dist' -  Dclean cleaned distance matrix obtained with dm_simmat
%
%     'Laplacian' -   Which laplacian matrix should be used. Choices are:
%        'rownorm' - row normalized Laplacian L = -Sij/sum(Sij) (i \ne j), Lii = 1;
%        'Lafon'  %from Lafon's presentation, L = Sij/sum(Sij)
%
%   Example:
%   [ev, aEV, aSClnd, indPos, D] = dm_dmsim(S);
%   [ev, aEV, aSClnd, indPos, D] = dm_dmsim(S, 'Laplacian', 'Lafon');
%
%  See also dm_dmit, dm_plot, dm_simmat
% References:
%   [1]

internal.stats.checkSupportedNumeric('S',S); % not complex
[n, p] = size(S);
if n == 0 || p == 0
    error('Input matrix S should not be empty');
end
if n ~=p
    error('Input matrix S should squared');
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% Parse arguments and check if parameter/value pairs are valid
%Norm  none/NORMALIZED/standardized
%Distance — Distance metric


paramNames = {'Dist',   'Laplacian'};
defaults   = {[], 'rownorm'};

[D, Laplacian, sf, rest]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

% Validate String value for  Laplacian value
Laplacians = {'rownorm'; 'Lafon';};
Laplacian = internal.stats.getParamVal(Laplacian,Laplacians,...
    '''Laplacian''');

%Check Dist matrix
internal.stats.checkSupportedNumeric('Dist',D); % not complex
if ~isempty(D)
    [n1, p1] = size(D);
    if n ~=n1 || p1 ~= p
        error('Input matrix Dist should have the same size as matrix S');
    end
end

switch Laplacian
    case 'rownorm'
        %L = 1 or -Sij/sum(Sij)
        for i = 1:size(S, 2)
            S(i, i) = 0;
        end
        Lnorm = sum(S, 2);
        Lnorm(Lnorm==0) = 1;
        L = -S./Lnorm;
        
        for i = 1:size(S, 2)
            L(i, i) = 1;
        end
        %eigenvalues
        [EVect, EVals] = eig(L);
        eval = (diag(EVals));
        eval = eval(:)';  %row vector
        %ev = real(diag(EVals));
        %Sort eigevalues and eigenvectors in ascending order
        [eval, indEvals] = (sort(eval));
        EVect = EVect(:, indEvals);
        
        %remove zero eigenvalues
        ind = eval < 1e-10;
        Components = sum(ind);  %number of diffusion map components
        eval = eval(~ind);
        EVect = EVect(:, ~ind);
        
        %devide by eigenvalues
        EVect2 = EVect./eval;
        
    case 'normkern'  %normilized kernel matrix, diffusion map https://en.wikipedia.org/wiki/Diffusion_map
% Step 1. Given the similarity matrix L.
% Step 2. Normalize the matrix according to parameter {\displaystyle \alpha }\alpha : {\displaystyle L^{(\alpha )}=D^{-\alpha }LD^{-\alpha }}L^{{(\alpha )}}=D^{{-\alpha }}LD^{{-\alpha }}.
% Step 3. Form the normalized matrix {\displaystyle M=({D}^{(\alpha )})^{-1}L^{(\alpha )}}M=({D}^{{(\alpha )}})^{{-1}}L^{{(\alpha )}}.
% Step 4. Compute the k largest eigenvalues of {\displaystyle M^{t}}M^{t} and the corresponding eigenvectors.
% Step 5. Use diffusion map to get the embedding {\displaystyle \Psi _{t}}\Psi _{t}.
    case 'Lafon'  %from Lafon's presentation, L = Sij/sum(Sij)
        for i = 1:size(S, 2)
            S(i, i) = 0;
        end
        Lnorm = sum(S, 2);
        Lnorm(Lnorm==0) = 1;
        L = S./Lnorm;
        %eigenvalues
        
        %eigenvalues
        [EVect, EVals] = eig(L);
        eval = (diag(EVals));
        eval = eval(:)';  %row vector
        %ev = real(diag(EVals));
        %Sort eigevalues and eigenvectors in ascending order
        [eval, indEvals] = sort(eval, 'descend');
        EVect = EVect(:, indEvals);
        
        %remove zero eigenvalues
        ind = abs(1-eval) < 1e-10;
        Components = sum(ind);  %number of diffusion map components
        eval = eval(~ind);
        EVect = EVect(:, ~ind);
        
        %multiply by eigenvalues
        EVect2 = EVect.*eval;
end

if nargout >= 5 % we need to find the diffusion distance in EVect2 space
    Ddiff =  squareform(pdist(EVect2(:, :)));
else
    Ddiff = [];
end

end
