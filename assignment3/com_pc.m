function [u,V,D]=com_pc(img)

[x,y]=meshgrid(0:(size(img,2)-1),0:(size(img,1)-1));
coord=[x(:) y(:)];
u=sum(repmat(img(:),1,2).*coord)/sum(img(:));
S=(repmat(img(:),1,2).*bsxfun(@minus,coord,u))'*bsxfun(@minus,coord,u)/sum(img(:));
[V,D]=eigs(S);
D=sqrt(diag(D));

figure;
imagesc(img); hold on;
plot(u(1),u(2),'xr');
plot([u(1) u(1)+D(1)*V(1,1)],[u(2) u(2)+D(1)*V(2,1)],'g');
plot([u(1) u(1)+D(2)*V(1,2)],[u(2) u(2)+D(2)*V(2,2)],'r');