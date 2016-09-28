C=chol([3 2; 2 3]);
x=randn(2000,2)*C;
x=bsxfun(@minus,x,mean(x));
[V,D]=eigs(x'*x/size(X,1));
D=sqrt(diag(D));

figure; hold on;
plot(x(:,1),x(:,2),'.b');
plot([0 D(1)*V(1,1)],[0 D(1)*V(2,1)],'g','LineWidth',2);
plot([0 D(2)*V(1,2)],[0 D(2)*V(2,2)],'r','LineWidth',2);
axis square;
xlim([-6 6]); ylim([-6 6]);

x2=x*100; x2=round(x2); x2=bsxfun(@minus,x2,min(x2))+1;
img=zeros(max(x2(:,2)),max(x2(:,1)));
for n=1:size(x2,1)
    img(x2(n,2),x2(n,1))=1;
end

com_pc(img);