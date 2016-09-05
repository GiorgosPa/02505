y = [3;1;1;2];
X = [1 1;1 -1;-1 1;-1 -1];

P = [ones(size(X,1),1) X]';
lambda = 0;

N = size(X, 1);
squaredDistance = zeros(N, N);
for coordinateNumber = 1 : 2
    squaredDistance = squaredDistance + ...
        (repmat(X(:, coordinateNumber), [1 N]) - ...
         repmat(X(:, coordinateNumber)', [N 1])).^2;
end
K = log( squaredDistance + eps ) .* squaredDistance / 2;

ab = [K P'; P zeros(size(P,1))] \ [y; zeros(size(P,1),1)];
ab2 = inv([K P'; P zeros(size(P,1))])*[y; zeros(size(P,1),1)];
ab3 = pinv([K P'; P zeros(size(P,1))])*[y; zeros(size(P,1),1)];

a = ab(1:size(X,1));
b = ab(size(X,1) + 1:size(ab,1));

%y2 = [K P']*ab;

points = 100;
Xnew = linspace(-1.5,1.5,points);
[X1, X2] = meshgrid(Xnew,Xnew);
Xnew = [X1(:) X2(:)];

N_a = size(Xnew, 1);
N_b = size(X,1);
squaredDistance = zeros(N_a, N_b);
for coordinateNumber = 1 : 2
    squaredDistance = squaredDistance + ...
        (repmat(Xnew(:, coordinateNumber), [1 N_b]) - ...
         repmat(X(:, coordinateNumber)', [N_a 1])).^2;
end

Kx = log( squaredDistance + eps ) .* squaredDistance / 2;
Px = [ones(size(Xnew,1),1) Xnew]';

y2 = [Kx Px'] * ab;

surf(X1,X2,reshape(y2,size(X1)));
hold on
plot3(X(1,1),X(1,2),y(1),'rx');
hold on
plot3(X(2,1),X(2,2),y(2),'rx');
hold on
plot3(X(3,1),X(3,2),y(3),'rx');
hold on
plot3(X(4,1),X(4,2),y(4),'rx');

al = zeros(size(X,1),4);
bl = zeros(size(X,2) + 1,4);
i = 1;
for lambda=[0.1,1,10,100]
    KL = K + lambda*eye(size(K,1));
    abl = [KL P'; P zeros(size(P,1))] \ [y; zeros(size(P,1),1)];
    al(:,i) = abl(1:size(X,1));
    bl(:,i) = abl(size(X,1) + 1:size(abl,1));
    i = i + 1;
end