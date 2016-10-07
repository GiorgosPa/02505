function [C,T]=compute_correspondance(c,d)

c0 = c-mean(c); % Centered shape
cs0 = c0/sqrt(c0'*c0); % Centered and unit size shape

d0 = d-mean(d); % Centered shape
ds0 = d0/sqrt(d0'*d0); % Centered and unit size shape

% Compute squared distance matrix
Nc = size(c,1);
Nd = size(d,1);
a = repmat(cs0,1,Nd);
b = repmat(ds0.',Nc,1);
e = a-b;
A = e.*conj(e);

% If A is non-square, add dummy points
if Nc>Nd
    A = [A 0.5*ones(Nc,Nc-Nd)];
elseif Nd>Nc
    A = [A; 0.5*ones(Nd-Nc,Nd)];
end

% Compute correspondances using Hungarian algorithm
[C,T] = hungarian(A);

