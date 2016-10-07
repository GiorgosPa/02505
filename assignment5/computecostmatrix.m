function cost = computecostmatrix(ctx1,ctx2,ndummy,dummycost)
%
% this function computes the cost matrix for shape context matching
% 'ctx1' and 'ctx2' are shape contexts for two points sets. 'ndummy' 
% is the number of dummy points, and 'dummycost' is teh cost for matching a
% point in one set to a dummy in the other
%
% Author: rl@imm.dtu.dk
%

k1= size(ctx1,1);
k2= size(ctx2,1);
maxk = max(k1,k2);
k = maxk+ndummy;
cost = [zeros(k1,k2) dummycost*ones(k1,k-k2); dummycost*ones(k-k1,k2) zeros(k-k1,k-k2)];

for i=1:k1,
    for j=1:k2,
       di = ctx1(i,:) - ctx2(j,:);
       su = ctx1(i,:) + ctx2(j,:);
       ndx = find(di ~= 0);
       cost(i,j) = 0.5*sum((di(ndx).^2)./su(ndx));
    end
end

