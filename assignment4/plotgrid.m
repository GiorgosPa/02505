function plotgrid(X,Y,varargin)

% Simple parser for optional variables
for narg=1:2:numel(varargin)
    eval([varargin{narg} '=varargin{' num2str(narg+1) '};']);
end
if ~exist('prune','var'), prune=1; end
if ~exist('imagesc','var'), imagesc=0; end

% Transpose X and Y and inverse Y so that the grid is displayed like
% imagesc
if imagesc
    X=X'; Y=Y';
    Y=-Y;
end

hold on;

% plot horizontal lines
for nx=1:prune:size(X,1)
    line([X(nx,1:end-1); X(nx,2:end)],[Y(nx,1:end-1); Y(nx,2:end)],'Color','r')
end

% plot vertical lines
for ny=1:prune:size(X,2)
    line([X(1:end-1,ny) X(2:end,ny)]',[Y(1:end-1,ny) Y(2:end,ny)]','Color','r')
end

axis image;



