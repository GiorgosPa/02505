function xx=com(X)

W=repmat(0:(size(X,2)-1),size(X,1),1).*X;
sumX=sum(X(:));
xx(1)=sum(W(:))/sumX;
W=repmat((0:(size(X,1)-1))',1,size(X,2)).*X;
xx(2)=sum(W(:))/sumX;
