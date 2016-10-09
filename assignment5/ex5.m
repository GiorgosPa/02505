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

% TO BE ADDED

% Compute cost matrix
dummycost=1;
A=computecostmatrix(sc1,sc2,0,dummycost);
figure;
subplot(1,3,1); imagesc(A); axis image; % Here we can see that we do not match dummy points
subplot(1,3,2); imagesc(A(1:200,1:50)); axis image;

% Compute correspondance using Hungarian algorithm
[C,T] = hungarian(A);
% ind=sub2ind(size(A),C(1:numel(d)),1:numel(d));
% A(ind)

subplot(1,3,3); hold on;
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
c_left_lung=1:ind;
c_right_lung=(ind+1):size(c,1);

d0 = d-mean(d); % Centered shape
ds0 = d0/sqrt(d0'*d0); % Centered and unit size shape
d_dist=abs(ds0(1:end-1)-ds0(2:end));
[~,ind]=max(d_dist);
d_left_lung=1:ind;
d_right_lung=(ind+1):size(d,1);

figure;
subplot(1,2,1); plot(c_dist); xlabel('Points difference'); ylabel('Distance');
subplot(1,2,2); plot(d_dist); xlabel('Points difference'); ylabel('Distance');

ind=ismember(c_right_lung,C(d_left_lung));
disp([num2str(sum(ind)) ' points of left lung were matched to right lung']);
ind=ismember(c_left_lung,C(d_right_lung));
disp([num2str(sum(ind)) ' points of right lung were matched to left lung']);

Alung=zeros(size(A));
for nc=c_left_lung, Alung(nc,d_right_lung)=1; end
for nc=c_right_lung, Alung(nc,d_left_lung)=1; end
Amasked=A;
Amasked(Alung==1)=1;

figure; subplot(1,2,1);
imagesc(Amasked); axis image;

% Compute correspondance using Hungarian algorithm
[C,T] = hungarian(Amasked);
% ind=sub2ind(size(A),C(1:numel(d)),1:numel(d));
% A(ind)

subplot(1,2,2); hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C(1:numel(d))) d].','r');
axis image;

%% Turn off points that move more than 20% along the outline

contour_pts_diff=@(x) abs([x(1:end-1)-x(2:end); x(end)-x(1)]);
a2polar=@(x) cos(x)+1i*sin(x);

thr=0.2;
l=a2polar(thr*2*pi);

lung={'left','right'};
total_d=0;
Acontour=zeros(size(A));
for nl=1:numel(lung)    
    % Map lungs points to circle
    d_lung_ind=eval(['d_' lung{nl} '_lung']);
    d_lung=ds0(d_lung_ind);
    d_diff=contour_pts_diff(d_lung);
    d_contour=sum(d_diff);
    d_circ=cumsum([0; d_diff(1:end-1)])/d_contour*2*pi;
    d_img=a2polar(d_circ);
    
    c_match=C(d_lung_ind(1));
    c_lung_ind=eval(['c_' lung{nl} '_lung']);
    ind=find(c_lung_ind==c_match);
    c_lung_new_ind=c_lung_ind([ind:end 1:(ind-1)]);
    c_lung=cs0(c_lung_new_ind);
    c_diff=contour_pts_diff(c_lung);
    c_contour=sum(c_diff);
    c_circ=cumsum([0; c_diff(1:end-1)])/c_contour*2*pi;
    c_img=a2polar(c_circ);
    
    for nd=1:numel(d_lung)
        ind=imag(c_img.*conj(d_img(nd))*l)>=0 & imag(c_img.*conj(d_img(nd))*conj(l))<=0;
        Acontour(c_lung_new_ind(~ind),total_d+nd)=1;       
    end
    total_d=numel(d_lung);
end

figure; subplot(1,2,1);
Amasked=A;
Amasked(Alung==1|Acontour==1)=1;
imagesc(Amasked); axis image;

% Compute correspondance using Hungarian algorithm
[C,T] = hungarian(Amasked);

subplot(1,2,2); hold on;
plot(c,'bx');
plot(d,'g+');
plot([c(C(1:numel(d))) d].','r');
axis image;

%% Optimize registration

computeP=@(x) [ones(size(x,1),1) real(x) imag(x)];
lambda=10;
N=numel(d);

src=[real(d) imag(d)];
dest=[real(c(C(1:50))) imag(c(C(1:50)))];
Cold=nan(50,1);
h=figure;
while ~isequal(sort(Cold),sort(C(1:50)))
    
    % Warp second shape to first one using TPS warp 
    P=[ones(N,1) src]';    
    K=computeK(src,dest);
    Y=[dest;zeros(3,2)];
    S=[K+lambda*eye(N) P'; P zeros(3,3)];
    coeff=S\Y;
    
    src=[K+lambda*eye(N) P']*coeff;
    src_img=src(:,1)+1i*src(:,2);

%     figure; hold on;
%     plot(c,'bx');
%     plot(src(:,1),src(:,2),'g+');
%     axis image;
    
    % Re-compute correspondences between the first shape and the warped
    % second shape using the Hungarian algorithm.
    dummycost=1;
    sc1=shapecontext(c,rfrac);
    sc2=shapecontext(src_img,rfrac);
    A=computecostmatrix(sc1,sc2,0,dummycost);
    Amasked=A;
    Amasked(Alung==1|Acontour==1)=1;
    Cold=C(1:50);
    [C,T] = hungarian(Amasked);
    dest=[real(c(C(1:50))) imag(c(C(1:50)))];
    
    figure(h);
    subplot(1,2,1); cla;
    imagesc(Amasked);
    
    subplot(1,2,2); cla; hold on;
    plot(d,'g+');
    plot(c,'bx');
    plot(src_img,'ms');    
    plot([d src_img].','r');
    axis image;
    legend({'source','target','transformed'})
    
    % Repeat until convergence
end