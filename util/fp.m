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
%   [H,P] = FP(...) returns both the figure and plot handles.
%   
%   See also HP.
% 
% Joe MacGregor (NASA)
% Last updated: 12 August 2024

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
addToolbarExplorationButtons(h)
switch nargin
    case 1
        if ~isempty(find(isnan(varargin{1}), 1))
            p               = line(1:length(varargin{1}), varargin{1}, 'Color', 'k', 'Marker', '.', 'MarkerSize', 14);
        else
            p               = line(1:length(varargin{1}), varargin{1}, 'Color', 'k', 'LineWidth', 2);
        end
    case 2
        if (all(diff(varargin{1}) > 0) && ~all(diff(diff(varargin{1})))) % check if x is monotonically increasing and uniformly spaced
            p               = line(varargin{1}, varargin{2}, 'Color', 'k', 'LineWidth', 2);
        else
            p               = line(varargin{1}, varargin{2}, 'Color', 'k', 'Marker', '.', 'MarkerSize', 14);
        end
end
set(gca, 'FontSize', 20)
switch nargin
    case 1
        xlabel('Index')
        ylabel(inputname(1), 'Interpreter', 'none')
    case 2
        xlabel(inputname(1), 'Interpreter', 'none')
        ylabel(inputname(2), 'Interpreter', 'none')
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