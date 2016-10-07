function ctx = shapecontext(shape,rfrac)
% 
% this function computes the shape context for one shape (set of points)
% the set of points 'shape' is a complex vector containing x-coordinates in
% the real part and the y-coordinates in the imaginary part.
% 'rfrac' is the size of the context as a
% proportion of the maximum range of the coordinates. 
%
% Author: rl@imm.dtu.dk
%

%n= size(shape,1);
%k=n/2;

%x = shape(1:k,1);
%y = shape(k+1:n,1);

k= size(shape,1);

x = real(shape);
y = imag(shape);


    Rmax = rfrac*max(range(x),range(y));
    Rlevels = 3;
    Rlimits = Rmax*(exp((0:Rlevels)/Rlevels)-1)/(exp(1)-1);
    Rlimits = Rlimits.*Rlimits;
    Rup = Rlimits(2:end);
    Rlow = Rlimits(1:end-1);
    thlevels = 8;
    thlimits = (-4:4)*pi/4;
    thup = thlimits(2:end);
    thlow =thlimits(1:end-1);
   
    c = x(:) + sqrt(-1)*y(:);
    a = repmat(c,1,k);
    b = repmat(c.',k,1);
    d = a-b;
    R = d.*conj(d);
    theta = atan2(imag(d),real(d));
    
ctx = zeros(k,Rlevels*thlevels);    
for j=1:k,    
    i = 1;
    for s=1:Rlevels,
        ndx = find(Rlow(s)<R(j,:) & R(j,:)<=Rup(s));
        for t=1:thlevels,
           ctx(j,i) = length(find(thlow(t)<theta(j,ndx) & theta(j,ndx)<=thup(t)));
           i = i+1;
        end
    end
end
ctx = ctx./(repmat(sum(ctx,2),1,size(ctx,2)));