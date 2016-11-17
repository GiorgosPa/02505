clear; close all;
load 'correctedData'
load 'segmentData'
maxiter = 100;
loglik = [];

%% plot the image and its histogram
figure;imagesc(data)
figure;imagesc(correctedData)

% initialize em parameters
minimum = min(correctedData(mask));
maximum = max(correctedData(mask));
p = [1/3, 1/3, 1/3];
width = (maximum - minimum)/3;
mu = [(minimum + width)/2,(2*minimum + 3*width)/2,(2*minimum + 5*width)/2];
sigma = [width^2,width^2,width^2];

% Compute some reasonable mean and variance for the Gaussian distribution
intensities = correctedData( find( mask ) );
mean = sum( intensities ) / length( intensities );
variance = sum( ( intensities - mean ).^2 ) / length( intensities );
% Compute the histogram and its properties
[ histogram binCenters ] = hist( intensities, 64 );
pdf = histogram / sum( histogram );
binSize = binCenters(2) - binCenters(1);
% Display
figure
bar( binCenters, pdf );
grid
hold on
gauss = 1/sqrt( 2 * pi * variance ) * ...
        exp( -( binCenters - mean ).^2 / 2 / variance );
plot( binCenters, gauss * binSize, 'g' )

figure
for iter=1:maxiter
    subplot(2,3,1)
    bar( binCenters, pdf );
    grid
    colors = ['r','g','y'];

    for i=1:3
        hold on
        gauss = 1/sqrt( 2 * pi * sigma(i) ) * ...
                exp( -( binCenters - mu(i) ).^2 / 2 / sigma(i) );
        plot( binCenters, gauss * binSize*p(i), colors(i) )
    end
    hold off

    % compute wik
    pr = zeros(numel(correctedData(mask)),3);
    for i=1:3
        pr(:,i) = 1/sqrt( 2 * pi * sigma(i) ) * ...
            exp( -( correctedData(mask) - mu(i) ).^2 / 2 / sigma(i) ) * p(i);
    end
    loglik(end+1) = sum(log(sum(pr,2)));
    w = pr./repmat(sum(pr,2),1,3);
    % display wik
    grshow=@(x) imshow(mat2gray(x));
    for i=1:3
        wim = zeros(240,160);
        wim(mask) = w(:,i);
        subplot(2,3,3+i); grshow(wim)
    end
    
    mu = sum(w.*repmat(correctedData(mask),1,3))./ sum(w);
    sigma = sum(w.*bsxfun(@minus,repmat(correctedData(mask),1,3),mu).^2) ./ sum(w);
    p = sum(w) / numel(correctedData(mask));
    subplot(2,3,2); plot(loglik)
    drawnow
end

%% vary only m2 from m1 to m3
m2 = linspace(mu(1),mu(3),100);

loglik2 = [];
lower_bound = [];

mu(2) = mu(1);
for i=1:3
    pr(:,i) = 1/sqrt( 2 * pi * sigma(i) ) * ...
        exp( -( correctedData(mask) - mu(i) ).^2 / 2 / sigma(i) ) * p(i);
end
w = pr./repmat(sum(pr,2),1,3);
for m = m2
    mu(2)=m;
    pr = zeros(numel(correctedData(mask)),3);
    for i=1:3
        pr(:,i) = 1/sqrt( 2 * pi * sigma(i) ) * ...
            exp( -( correctedData(mask) - mu(i) ).^2 / 2 / sigma(i) ) * p(i);
    end
    loglik2(end+1) = sum(log(sum(pr,2)));
    lower_bound(end+1) = sum(sum(w.*log(pr./w)));
end
munew = sum(w.*repmat(correctedData(mask),1,3))./ sum(w);
[~,ind] = min(abs(m2 - munew(2)));
figure; plot(m2,loglik2);
hold on
plot(m2,lower_bound,'-.');
plot(munew(2),lower_bound(ind),'xb','markersize',10)


