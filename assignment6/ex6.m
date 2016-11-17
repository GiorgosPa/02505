clear all;
close all;
load abdomen_biascorrected

%% select a rectangular area of interest
x = find(max(abdomencorrected,[],2));
y = find(max(abdomencorrected));
aoi = abdomencorrected(x,y);
imagesc(aoi)

% polar representation
% Allocate arrays for pixel values and coordinates of lines extending
% from image center to the border as the spokes of a wheel
[nrows, ncols] = size(aoi);
maxlen = ceil(0.5*norm([nrows ncols]));
nang = 360;
lines = zeros(nang,maxlen); % Pixel values along wheel spokes
lcoords = zeros(nang,maxlen,2);


for I=1:nang,
  % Compute unit vector in direction of spoke
  t = [cos(2*pi*I/nang+pi/2) sin(2*pi*I/nang+pi/2)];
  % Compute coordinates of pixels along spoke
  linecoords = ones(maxlen,1)*[ncols nrows]/2 + [0:maxlen-1]'*t;
  % Store coordinates as pixel index in original image
  lcoords(I,:,:) = linecoords;
  % Extract pixel values using bilinear interpolation
  lines(I,:) = interp2(aoi,linecoords(:,1),linecoords(:,2),'linear')';
end

% Make sure we don?t have NaN values anywhere
lines(find(isnan(lines))) = 0;
% Re-order the angles
lines = fliplr(lines);
lcoords = flipdim(lcoords,2);

% extract the radial derivatives
edges = diff(lines')';
edges = abs(edges);
edges = 1 - edges;

path = findpath(edges);

% plot polar
figure;
imagesc(lines)
hold on
plot(path, 1:360, '.r','markersize', 12)

% plot regular image
pixels = zeros(360,2);
for i=1:360
    pixels(i,:) = lcoords(i,path(i),:);
end
figure;
imagesc(aoi);
hold on
plot(pixels(:,1), pixels(:,2), '.r', 'markersize', 12)

%% modify polar representation
mpolar = .5*ones(size(lines));
for i=1:size(lines,1)
    mpolar(i,1:(128 - path(i))) = lines(i,path(i)+2:end);
    lcoords(i,1:(128 - path(i)),:) = lcoords(i,path(i)+2:end,:);
end


edges = diff(mpolar')';
edges = abs(edges);
edges = 1 - edges;
%edges = edges(:,2:end);

figure;
imagesc(mpolar);

path = findpath(edges);

% plot polar
figure;
imagesc(mpolar)
hold on
plot(path, 1:360, '.r','markersize', 12)

% plot regular image
pixels = zeros(360,2);
for i=1:360
    pixels(i,:) = lcoords(i,path(i),:);
end
figure;
imagesc(aoi);
hold on
plot(pixels(:,1), pixels(:,2), '.r', 'markersize', 12)


