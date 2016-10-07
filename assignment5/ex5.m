%% Compute squared distance matrix
load meta
c=meta(:,1);
d=meta(:,2);
[C,T]=compute_correspondance(c,d);

figure; hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C) d].','r');
axis image;

%% Subsample lung2
load lung
c = lung1;
d = lung2(1:4:200);
[C,T]=compute_correspondance(c,d);

% Plot actual matches (not dummy)
figure; hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C(1:numel(d))) d].','r');
axis image;

%% Compute shape context (rfrac=0.25)
rfrac=0.25;
sc1=shapecontext(c,rfrac);
sc2=shapecontext(d,rfrac);

% Plot the context of points having similar context

% Compute cost matrix
dummycost=1;
A=computecostmatrix(sc1,sc2,150,dummycost);
figure;
subplot(1,2,1); imagesc(A); % Here we can see that we do not match dummy points
subplot(1,2,2); imagesc(A(1:200,1:50));

% Compute correspondance using Hungarian algorithm
[C,T] = hungarian(A);

figure; hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C(1:numel(d))) d].','r');
axis image;

%% Now, we have some extra information we can use, namely that the lung dataset
% consist of a left and a right part and the points are sampled equidistantly along
% the outlines.
% Implement these constraints by setting the cost in the cost matrix
% such that points from the left lung cannot match to the right lung and vice
% versa, and such that points from one lung dataset cannot shift more than ±20%
% along the outline.

% Define indices for left and right lungs
c0 = c-mean(c); % Centered shape
cs0 = c0/sqrt(c0'*c0); % Centered and unit size shape
c_dist=abs(cs0(1:end-1)-cs0(2:end));
[~,ind]=max(c_dist);
c_left_lung=1:(ind-1);
c_right_lung=ind:size(c,1);

d0 = d-mean(d); % Centered shape
ds0 = d0/sqrt(d0'*d0); % Centered and unit size shape
d_dist=abs(ds0(1:end-1)-ds0(2:end));
[~,ind]=max(d_dist);
d_left_lung=1:(ind-1);
d_right_lung=ind:size(d,1);

ind=ismember(c_right_lung,C(d_left_lung));
disp([num2str(sum(ind)) ' points of left lung were matched to right lung']);
ind=ismember(c_left_lung,C(d_right_lung));
disp([num2str(sum(ind)) ' points of right lung were matched to left lung']);

for nc=c_left_lung, A(nc,d_right_lung)=100; end
for nc=c_right_lung, A(nc,d_left_lung)=100; end

imagesc(A);

% Compute correspondance using Hungarian algorithm
[C,T] = hungarian(A);

figure; hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C(1:numel(d))) d].','r');
axis image;

% Turn of points that move more than 20% along the outline

% Plot new cost matrix
% imagesc(Anew);

%% Optimize registration

computeP=@(x) [ones(1,size(x,1)); x'];
dold=zeros(size(d)); dnew=d; 
while sum(abs(dnew-dole))>0.001

% Warp second shape to first one using TPS warp
K=computeK(c(C(1:50)),d);
S=[K+lambda*eye(N) P'; P zeros(p+1,p+1)];
coeff=S\[y;zeros(p+1,1)];
dold=dnew;
dnew=

% Re-compute correspondences between the first shape and the warped
% second shape using the Hungarian algorithm.
[C,T] = hungarian(A);

% Repeat until convergence
end