function varargout          = hp(varargin)
% HP Add another plot to an existing figure.
% 
%   HP(X) plots vector X onto of the most recent figure (or a new figure if
%   none presently exist). If X is monotonically increasing and uniformly
%   spaced, then it is plotted as a red line. Otherwise it is plotted as
%   red dots.
%   
%   HP(X,Y) plots vector X vs Y.
%   
%   HP(...,C) plots using the color C instead of red, where C must be a
%   single-character string that is one of the standard MATLAB colors (see
%   HELP PLOT).
%   
%   P = HP(...) returns the plot handle.
%   
%   See also FP.
%   
% Joe MacGregor (UTIG)
% Last updated: 03/13/13

if ((nargin < 1) || (nargin > 3))
    error('hp:nargin', 'Incorrect number of inputs to HP (must be between 1 and 3).')
end
if (nargout > 1)
    error('hp:nargout', 'Too many outputs requested from HP (NARGOUT > 1).')
end

if ischar(varargin{end})
    if ~any(strcmp(varargin{end}, {'b' 'g' 'r' 'c' 'm' 'y' 'k' 'w'}))
        error('hp:colors', 'Plot color must be one of the MATLAB standards (see HELP PLOT).')
    end
    num_var = nargin - 1;
else
    num_var = nargin;
end

if ((num_var == 2) && (length(varargin{1}) ~= length(varargin{2})))
    error('hp:lengths', 'Input vectors not the same length.')
end

hold on
switch num_var
    case 1
        p                   = plot(varargin{1}, 'linewidth', 2);
        if (num_var < nargin)
            set(p, 'color', varargin{end})
        else
            set(p, 'color', 'r')
        end
    case 2
        if (all(diff(varargin{1}) > 0) && ~all(diff(diff(varargin{1})))); % check if x is monotonically increasing and uniformly spaced
            p               = plot(varargin{1}, varargin{2}, 'linewidth', 2);
        else
            p               = plot(varargin{1}, varargin{2}, '.', 'markersize', 14);
        end;
        if (num_var < nargin)
            set(p, 'color', varargin{end})
        else
            set(p, 'color', 'r')
        end
end

if nargout
    varargout{1}            = p;
end