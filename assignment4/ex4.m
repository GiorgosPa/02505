%% Construct a set of 2D tensor B-spline basis functions

%image dimensions
m = [20 20];
% number of knots
p1 = 3; p2 = 3;

% knot sequence 1
k1 = linspace(1,m(1),p1); k1 = augknt(k1,3);

% knot sequence 2
k2 = linspace(1,m(2),p2); k2 = augknt(k2,3);

% Form B-spline matrix Q which contains basic functions
B1 = spmak(k1,eye(p1)); Q1 = fnval(B1,1:m(1))';
%
B2 = spmak(k2,eye(p2)); Q2 = fnval(B2,1:m(2))';

Q = kron(speye(2),kron(Q2,Q1));

% Set all elements of w1 to a values and elements of w2 to 0
w = [10*ones(9,1); zeros(9,1)];
[x1,x2] = meshgrid(1:m(1),1:m(2));
y = [x1(:); x2(:)] + Q*w;
y1 = reshape(y(1:end/2),size(x1,1),size(x1,2));
y2 = reshape(y(end/2+1:end),size(x1,1),size(x1,2));
% figure;
% subplot(1,4,1); plotgrid(x1,x2);
% subplot(1,4,2); plotgrid(y1,y2);

%plot 1D splines
figure; plot(Q1);
%plot 2D splines

for i=1:3
    for j=1:3
        figure; surf(Q1(:,i)*Q1(:,j)');
    end
end

figure; 
plotgrid(x1,x2);axis square;
figure;
plotgrid(y1,y2);axis square;

% Set all elements of w1 to a values and elements of w2 to 0
w = [zeros(9,1); 10*ones(9,1)];
[x1,x2] = meshgrid(1:m(1),1:m(2));
y = [x1(:); x2(:)] + Q*w;
y1 = reshape(y(1:end/2),size(x1,1),size(x1,2));
y2 = reshape(y(end/2+1:end),size(x1,1),size(x1,2));
% subplot(1,4,3); plotgrid(y1,y2);

figure;
plotgrid(y1,y2);axis square;

% Set all elements of w1 and w2 to 0 except for the 9th basis
w = [zeros(8,1); 10; ones(8,1); 10];
[x1,x2] = meshgrid(1:m(1),1:m(2));
y = [x1(:); x2(:)] + Q*w;
y1 = reshape(y(1:end/2),size(x1,1),size(x1,2));
y2 = reshape(y(end/2+1:end),size(x1,1),size(x1,2));
% subplot(1,4,4); plotgrid(y1,y2);

figure;
plotgrid(y1,y2);axis square;

%% Perform non-linear intensity-based registration with Gauss-Newton optimization
close all
load mr.mat
m=size(mr1);
mr_diff=abs(mr1-mr2);

% Specify voxels coordinates
[x1,x2]=meshgrid(1:m(2),1:m(1));

% Create basis functions
m = size(mr1);
alpha=10;
p1 = 5; p2 = 5; % Set the number of basis functions
k1 = linspace(1,m(1),p1); k1 = augknt(k1,3);
k2 = linspace(1,m(2),p2); k2 = augknt(k2,3);
B1 = spmak(k1,eye(p1)); Q1 = fnval(B1,1:m(1))';
B2 = spmak(k2,eye(p2)); Q2 = fnval(B2,1:m(2))';
Q = kron(speye(2),kron(Q2,Q1));

% Set regularizer matrix
G=eye(size(Q,2));

w=zeros(size(Q,2),1);


thr=0.0001;
maxiter=500;
newTy=[];
I=nan(1,maxiter);
iter=1;

doPlot=false; % Turn on or off diagnostic plotting

if doPlot
    nupdate=250; %#ok<*UNRCH>
    h=figure;
    subplot(2,3,1); imagesc(mr1); axis image; axis off;
end

% figure; imagesc(mr1); axis image; axis off;
% figure; imagesc(mr2); axis image; axis off;

while true
    
    % If Ty was computed from previous iteration, reuse it
    if isempty(newTy) 
        Ty=computeTy(mr2,x1,x2,Q,w);
    else
        Ty=newTy;
    end
    [dxTy,dyTy]=gradient(Ty);
    dTy=[dxTy(:); dyTy(:)];
    
    % Solve for update dw using Eq. (2.34)
    A=repmat(dTy,1,size(Q,2)).*Q;    
    dw=pinv(A'*A+alpha*(G'*G))*(A'*repmat(mr1(:)-Ty(:),2,1)-alpha*(G'*G)*w);
    
    % Update w := w + dw
    
    % Compute new value of the objective function
    [newTy,y1,y2]=computeTy(mr2,x1,x2,Q,w+dw);
    I(iter)=0.5*(norm(mr1(:)-newTy(:))^2+alpha*norm(G*(w+dw))^2);
        
    % If objective function increases, do half-step
    if iter>1 && I(iter-1)<I(iter); 
        w=w+dw/2;
        [newTy,y1,y2]=computeTy(mr2,x1,x2,Q,w);
        I(iter)=0.5*(norm(mr1(:)-newTy(:))^2+alpha*norm(G*w)^2); % Update the objective function
    else
        w=w+dw;
    end
    if iter>1
        disp(['Iter ' num2str(iter) ', I=' num2str(I(iter)) ', Idiff=' num2str(abs(I(iter-1)-I(iter)))]);
    end

    % Plot things for diagnostic purposes
    if doPlot && mod(iter,nupdate)==1 && iter~=1
        figure(h);
        subplot(2,3,2); imagesc(Ty); axis image; axis off;
        img_diff=abs(mr1-newTy);
        diff_scale=max([mr_diff(:); img_diff(:)]);
        subplot(2,3,4); imshow(mr_diff/diff_scale); axis image; axis off;
        subplot(2,3,5); imshow(img_diff/diff_scale); axis image; axis off; 
        ha=subplot(2,3,3); cla; plotgrid(y1,y2,'prune',7,'imagesc',1); axis off;
        subplot(2,3,6); plot(I); axis tight;
        drawnow;
    end
    
    % End-loop conditions
    if iter>1 && (abs(I(iter-1)-I(iter))/abs(I(iter))<thr || iter==maxiter)
        break;
    end
    iter=iter+1;
end
if iter==maxiter&&abs(I(iter-1)-I(iter))>thr
    warning('Did not converge');
end

figure; imagesc(newTy); axis image; axis off;
img_diff=abs(mr1-newTy);
diff_scale=max([mr_diff(:); img_diff(:)]);
figure; imshow(mr_diff/diff_scale); axis image; axis off;
figure; imshow(img_diff/diff_scale); axis image; axis off;
figure; cla; plotgrid(y1,y2,'prune',7,'imagesc',1); axis off;
figure; plot(I); axis tight;
fusionRGB=cat(3,mr1/256,newTy/256,zeros(size(newTy)));
figure; imshow(fusionRGB);

fusionRGB=cat(3,mr1/256,mr2/256,zeros(size(newTy)));
figure; imshow(fusionRGB);


disp(['Initial squared norm ' num2str(norm(mr1(:)-mr2(:))^2)]);
disp(['Final squared norm ' num2str(norm(mr1(:)-newTy(:))^2)]);

