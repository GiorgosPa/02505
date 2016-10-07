function K=computeK(x,y)

p=size(x,2);
l2=@(x) sqrt(sum(x.^2,2)); % assumes rows are samples and columns dimensions
switch p
    case 2
        eta=@(r) r.^2.*log(r);
    case 3
        eta=@(r) r.^3;
    otherwise
        error('Invalid dimension p');
end

K=nan(size(x,1),size(y,1));
for j=1:size(y,1)
    K(:,j)=eta(l2(bsxfun(@minus,x,y(j,:))));
end
K(isnan(K))=0;