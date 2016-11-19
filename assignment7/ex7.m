%% 1-Select 20 points inside tumor (label 1) and compute mean and standard deviation

clear all; close all;

K=2;

load('tumorData');
data=double(data);

load('points');
imgray(data);
% [x,y]=ginput(20);
hold on;
plot(x,y,'m.','LineWidth',4);
% save('points.mat','x','y');

vals=data(sub2ind(size(data),ceil(y),ceil(x)));
m(1)=mean(vals);
s(1)=std(vals);

%% 2-Compute mean and standard deviation of whole image (label 2)

ind=data(:)>10; % Ignore obvious background voxels
m(2)=mean(data(ind));
s(2)=std(data(ind));

%% Show PDFs
figure; hold on;
[N,edges]=hist(data(:),50);
bar(edges,N/trapz(edges,N));
x=min(data(:)):0.1:max(data(:));
prior=[0.2 0.8];
y=nan(numel(x),2);
for n=1:K
    y(:,n)=normpdf(x,m(n),s(n))*prior(n);
end
l=plot(x,y,'LineWidth',2);
l(3)=plot(x,sum(y,2),'k','LineWidth',2);
xlim([0, max(data(:))]);
legend(l,{'Tumor PDF','Image PDF','Total PDF'});

%% 3-Fit GMM (K=2) to the data

% Compute posteriors according to mean and standard deviation estimated in
% 1-2
prior=[0.2 0.8];
post=nan(numel(data),K);
for n=1:K
    norm_const=arrayfun(@(x,y) normpdf(data(:),m(y),s(y))*prior(x),n*ones(1,K),1:K,'UniformOutput',0);
    post(:,n)=(normpdf(data(:),m(n),s(n))*prior(n))./sum(cat(2,norm_const{:}),2);
end

% Make final label assignment according to MAP
[~,labels]=max(post,[],2);
labels=reshape(labels,size(data));

figure;
for n=1:K
    subplot(1,3,n);
    imagesc(reshape(post(:,n),size(data))); axis image; axis off;
end
subplot(1,3,3);
tmp = repmat(data,[ 1 1 3 ]);
red = zeros([size(data) 3]);
red(:,:,1) = 255;
tmp = tmp.*repmat(labels==2,[1 1 3]) + ...
    red .* repmat(labels==1, [1 1 3]);
imshow(uint8(tmp)); axis image; axis off;

%% 4-5-Compute mean-field approximation of the segmentatio posterior using the iterative minimization algorithm

% Define neighbor indices
tmp1 = ones(size(data,1),1);
tmp1(1:2:end) = 0;
tmp2 = ones(1,size(data,2));
tmp2(1:2:end) = 0;
whiteIndices = logical(kron(tmp1,tmp2));
blackIndices = logical(kron(~tmp1,~tmp2));
redIndices = logical(kron(~tmp1,tmp2));
greenIndices = logical(kron(tmp1,~tmp2));
indicesToUpdateCell{1}=find(whiteIndices);
indicesToUpdateCell{2}=find(blackIndices);
indicesToUpdateCell{3}=find(redIndices);
indicesToUpdateCell{4}=find(greenIndices);

neighborMap=zeros([numel(data),3]);
neighborMap(whiteIndices,:)=1;
neighborMap(redIndices,1)=1;
neighborMap(greenIndices,2)=1; 
neighborMap=reshape(neighborMap,[size(data),3]);
figure; imshow(neighborMap);

% for beta=[0 0.3 0.6 2.5];
for beta=[5];
        
    % Use GMM posterior as initialization for mean-field approximation q. This
    % assumes that the responsabilities are initially equal to the label prior
    q=post;

    % Compute responsabilities and update q until convergence (label assignment for is voxel is unchanged)
    lpost=labels(:); lpre=zeros(numel(data),1); diffNeigh=zeros(numel(data),K); resp=zeros(numel(data),1);
    niter=1;
    while sum(lpre~=lpost)>500 || niter<100
        % Compute neighbor differences
        if niter==1
            for n=1:K
                nonPadded=1-q(:,n);
                padded=zeros(size(nonPadded)+2);
                padded(2:end-1,2:end-1)=nonPadded;
                diffNeigh(:,n)=padded([2:end-1]-1,[2:end-1]) + ...
                    padded([2:end-1]+1,[2:end-1])+ ...
                    padded([2:end-1],[2:end-1]-1)+ ...
                    padded([2:end-1],[2:end-1]+1)+ ...
                    padded([2:end-1]-1,[2:end-1]-1)+ ...
                    padded([2:end-1]-1,[2:end-1]+1)+ ...
                    padded([2:end-1]+1,[2:end-1]-1)+ ...
                    padded([2:end-1]+1,[2:end-1]+1); %#ok<NBRAK>
            end
        end
        
        % Update responsabilities and posterior of non neighboring voxels
        for n=1:numel(indicesToUpdateCell)
            ind=indicesToUpdateCell{n};
            
            % Update responsabilities
            neigh_prior=arrayfun(@(x) prior(x)*exp(-beta*diffNeigh(ind,x)),1:K,'UniformOutput',0);
            norm_const=sum(cat(2,neigh_prior{:}),2);
            resp=cat(2,neigh_prior{:})./repmat(norm_const,1,K);
            
            % Update posterior
            like_prior=arrayfun(@(x) resp(:,x).*normpdf(data(ind),m(x),s(x)),1:K,'UniformOutput',0);
            norm_const=sum(cat(2,like_prior{:}),2);
            q(ind,:)=cat(2,like_prior{:})./repmat(norm_const,1,K);
            
            % Update neighborhood differences with new posterior
            for k=1:K
                nonPadded=1-q(:,k);
                padded=zeros(size(nonPadded)+2);
                padded(2:end-1,2:end-1)=nonPadded;
                diffNeigh(:,k)=padded([2:end-1]-1,[2:end-1]) + ...
                    padded([2:end-1]+1,[2:end-1])+ ...
                    padded([2:end-1],[2:end-1]-1)+ ...
                    padded([2:end-1],[2:end-1]+1)+ ...
                    padded([2:end-1]-1,[2:end-1]-1)+ ...
                    padded([2:end-1]-1,[2:end-1]+1)+ ...
                    padded([2:end-1]+1,[2:end-1]-1)+ ...
                    padded([2:end-1]+1,[2:end-1]+1); %#ok<NBRAK>
            end
        end
        
        % Update label to MAP
        lpre=lpost;
        [~,lpost]=max(q,[],2);
        
        niter=niter+1;
    end
    
    fprintf('Niter=%d\n',niter);
    
    % Plot posteriors and overlay
    figure;
    for n=1:K
        subplot(1,3,n);
        imagesc(reshape(q(:,n),size(data))); axis image; axis off;
    end
    subplot(1,3,3);
    tmp = repmat(data,[ 1 1 3 ]);
    red = zeros([size(data) 3]);
    red(:,:,1) = 255;
    rlpost=reshape(lpost,size(data));
    tmp = tmp.*repmat(rlpost==2,[1 1 3]) + ...
        red .* repmat(rlpost==1, [1 1 3]);
    imshow(uint8(tmp)); axis image; axis off;
       
    % 6-Expected area of tumor
    V=sum(q(:,1));
    fprintf('Expected volume of tumor (B=%0.1f)=%0.4f\n',beta,V);
end
