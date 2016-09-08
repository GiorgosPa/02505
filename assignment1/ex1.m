%% Simulation

close all;

% Set up parameters
y=[3;1;1;2];
x=[1 1; 1 -1; -1 1; -1 -1];
N=length(y);
p=size(x,2);
K=computeK(x,x);
computeP=@(x) [ones(1,size(x,1)); x'];
P=computeP(x);
    
[X,Y]=meshgrid(-1.5:0.05:1.5,-1.5:0.05:1.5);
    
for lambda=[0 1 10 100 1e7]
    
    % Compute the interpolating and smooth TPS
    S=[K+lambda*eye(N) P'; P zeros(p+1,p+1)];
    coeff=S\[y;zeros(p+1,1)];
    Slam=[K P']*pinv(S);
    
    % Plot f on a densely sampled grid
    Z=reshape([computeK([X(:) Y(:)],x) computeP([X(:) Y(:)])']*coeff,size(X));
    figure;
    surf(X,Y,Z);
    title(['lambda=' num2str(lambda) ', Deg Freedom=' num2str(trace(Slam(:,1:N)))]);
end

%% Bias-field correction

% Pick 100 coordinates in fat voxels
load('abdomen.mat'); abdomen=double(abdomen); roi=logical(roi);
figure;
imagesc(abdomen);
[px,py]=ginput(1);
save('coordinates.mat','px','py');

% Show images with selected points
load('coordinates.mat');
figure; imagesc(abdomen); hold on; plot(px,py,'x');
x=round([px py]);
y=double(abdomen(sub2ind(size(abdomen),x(:,2),x(:,1)))); % Inverse x,y to get correct values

N=length(y);
p=size(x,2);
K=computeK(x,x);
computeP=@(x) [ones(1,size(x,1)); x'];
P=computeP(x);

% Find lambda for 5, 10 and 20 deg freedom
lambda=700;
opt_lambda = [];
for target_df=[20 10 5]
    df=inf;
    while df>target_df
        % Compute smooth TPS
        S=[K+lambda*eye(N) P'; P zeros(p+1,p+1)];
        coeff=S\[y;zeros(p+1,1)];
        Slam=[K P']*pinv(S);        
        df=trace(Slam(:,1:N));
        lambda=lambda+1;
    end
    opt_lambda(end+1) = lambda; %#ok<SAGROW>
    disp(['Lambda=' num2str(lambda) ', df=' num2str(df)]);
end

%% Plot solution for 5, 10 and 20 deg freedom
[X,Y]=meshgrid(1:size(abdomen,1),1:size(abdomen,2));
basis=[computeK([X(:) Y(:)],x) computeP([X(:) Y(:)])'];

grshow=@(x) imshow(mat2gray(x));
figure; grshow(abdomen);

for lambda=opt_lambda
    S=[K+lambda*eye(N) P'; P zeros(p+1,p+1)];
    coeff=S\[y;zeros(p+1,1)];
    bias=reshape(basis*coeff,size(X));
    figure; grshow(bias);
    abdomen_corr=abdomen; abdomen_corr=abdomen_corr.*roi; 
    abdomen_corr(roi)=abdomen(roi)./bias(roi);
    figure; grshow(abdomen_corr);
end