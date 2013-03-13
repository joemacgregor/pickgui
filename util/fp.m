function varargout          = fp(varargin)
% FP Make a quick, clean figure.
% 
%   FP(X) plots vector X cleanly in a new figure ("first plot" or
%   "figure/plot"). If X is monotonically increasing and uniformly spaced,
%   then it is plotted as a black line. Otherwise it is plotted as black
%   dots. If resolvable by MATLAB, the variable name of X will be given on
%   the y-axis.
%   
%   FP(X,Y) plots X versus Y. If resolvable by MATLAB, the variable name of
%   X will be given on the x-axis and that of Y will be given on the y-axis.
%   
%   H = FP(...) returns the figure handle.
%   
%   [H, P] = FP(...) returns both the figure and plot handles.
%   
%   See also HP.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/13/13

if ((nargin < 1) || (nargin > 2))
    error('fp:nargin', 'Incorrect number of inputs to FP (must be 1 or 2).')
end
if (nargout > 2)
    error('fp:nargout', 'Too many outputs requested from FP (NARGOUT > 2).')
end
if ((nargin == 2) && (length(varargin{1}) ~= length(varargin{2})))
    error('fp:lengths', 'Input vectors not the same length.')
end

h                           = figure;
set(gca, 'fontsize', 20)
switch nargin
    case 1
        p                   = plot(varargin{1}, 'k', 'linewidth', 2);
        xlabel('Index')
        ylabel(inputname(1), 'interpreter', 'none')
    case 2
        if (all(diff(varargin{1}) > 0) && ~all(diff(diff(varargin{1})))) % check if x is monotonically increasing and uniformly spaced
            p               = plot(varargin{1}, varargin{2}, 'k', 'linewidth', 2);
        else
            p               = plot(varargin{1}, varargin{2}, 'k.', 'markersize', 14);
        end
        xlabel(inputname(1), 'interpreter', 'none')
        ylabel(inputname(2), 'interpreter', 'none')
end
grid on
box on
hold on

switch nargout
    case 1
        varargout{1}        = h;
    case 2
        [varargout{1}, varargout{2}] ...
                            = deal(h, p);
end