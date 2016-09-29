function [Ty,y1,y2]=computeTy(img,x1,x2,Q,w)

min_x1=min(x1(:)); max_x1=max(x1(:));
min_x2=min(x2(:)); max_x2=max(x2(:));

m=size(img);
y=[x1(:); x2(:)] + Q*w;
y1 = reshape(y(1:end/2),m);
y2 = reshape(y(end/2+1:end),m);
y1(y1<min_x1)=min_x1; y1(y1>max_x1)=max_x1;
y2(y2<min_x2)=min_x2; y2(y2>max_x2)=max_x2;
ind = sub2ind(m,round(y2(:)),round(y1(:)));
Ty=img(ind);
Ty=reshape(Ty,m);