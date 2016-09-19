function [Um,Sm]=principal_axis(X,m)

Xn=(X-min(X(:)))/(max(X(:))-min(X(:)));
coord=(0:(size(X,2)-1))-m(1);
Wx=repmat(coord,size(X,1),1).*Xn;
coord=(0:(size(X,1)-1))-m(2);
Wy=repmat(coord',1,size(X,2)).*Xn;
W=[Wx(:) Wy(:)];
[Um,Sm,~]=svd(W','econ');
Sm=diag(Sm);
figure; plot(W(:,1),W(:,2),'.'); hold on;
plot([0 sqrt(Sm(1))*Um(1,1)],[0 sqrt(Sm(1))*Um(2,1)],'g');
plot([0 sqrt(Sm(2))*Um(1,2)],[0 sqrt(Sm(2))*Um(2,2)],'m');
axis square