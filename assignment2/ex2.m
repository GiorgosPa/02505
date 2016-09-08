% select 10 points
load('visiblehuman');

% Each member selects points select
cpselect(head,head_frozen);
cpselect(head,head_fresh);
cpselect(head,head_mri);

save('selected_points_head1_frozen','ctrl_head1','ctrl_frozen', ...
    'ctrl_head2','ctrl_fresh','ctrl_head3','ctrl_mri');

% Compute FLE
N=10; J=3; D=2;
X=cat(3,ctrl_head1,ctrl_head2,ctrl_head3);
u=mean(X,3);
var_FLE=sum(reshape((X-repmat(u,[1,1,J])).^2,1,N*J*D))/(2*N*J);

% Estimate parameters from MRI to frozen CT and reverse
x=ctrl_mri;
xc=bsxfun(@minus,x,mean(x));
y=ctrl_frozen;
yc=bsxfun(@minus,y,mean(y));
H=xc'*yc;
[U,S,V]=svd(H);
R=V*diag([1,det(V*U)])*U';
s=sum(diag(xc*R'*yc'))/sum(diag(x*x'));
t=mean(y)'-s*R*mean(x)';
yhat=(R*x'+repmat(t,1,size(x,1)))';
FRE_mri=sum(reshape((yhat-y).^2,1,numel(y)));

figure;
imagesc(head_frozen);
hold on;
plot(y(:,1),y(:,2),'mx');
hold on;
plot(yhat(:,1),yhat(:,2),'gx');
