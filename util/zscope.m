% ZSCOPE Plot a quick z-scope echogram
% 
% Joe MacGregor
% Last updated: 12 August 2024

function varargout          = zscope(varargin)

switch nargin
    case 1
        data                = varargin{1};
        disp_scale          = [min(data(:), [], 'omitnan') max(data(:), [], 'omitnan')];
    case 2
        data                = varargin{1};
        disp_scale          = varargin{2};
        if isscalar(disp_scale)
            disp_scale      = [-disp_scale disp_scale];
        end
    case 3
        x                   = varargin{1};
        y                   = varargin{2};
        data                = varargin{3};
        disp_scale          = [min(data(:), [], 'omitnan') max(data(:), [], 'omitnan')];
        if isempty(x)
            x               = 1:size(data, 2);
        end
        if isempty(y)
            y               = 1:size(data, 1);
        end
    case 4
        x                   = varargin{1};
        y                   = varargin{2};
        data                = varargin{3};
        disp_scale          = varargin{4};
        if isempty(x)
            x               = 1:size(data, 2);
        end
        if isempty(y)
            y               = 1:size(data, 1);
        end
        if isscalar(disp_scale)
            disp_scale      = [-disp_scale disp_scale];
        end
end

figure
if (nargin <= 2)
    imagesc(data, disp_scale)
else
    imagesc(x, y, data, disp_scale)
end
colormap(bone)
set(gca, 'FontSize', 20)
if (nargin <= 2)
    title(inputname(1), 'Interpreter', 'none')
else
    xlabel(inputname(1), 'Interpreter', 'none')
    ylabel(inputname(2), 'Interpreter', 'none')
    title(inputname(3), 'Interpreter', 'none')
end
colorbar('FontSize', 20)
grid on
box on

switch nargout
    case 1
        varargout{1}        = gcf;
    case 2        
        [varargout{1}, varargout{2}] ...
                            = deal(gcf, gca);
end