load 'visiblehuman'
head_mri=double(head_mri);
head_frozen=double(head_frozen);

% Get center of mass, principal axis and lengths
[uf,Vf,Df]=com_pc(head_frozen);
y=[uf; uf+Df(1)*Vf(:,1)'; uf+Df(2)*Vf(:,2)'];
hold on;
for i=1:3
    plot(y(i,1),y(i,2),'m+')
end
[um,Vm,Dm]=com_pc(head_mri);
x=[um; um+Dm(1)*Vm(:,1)'; um+Dm(2)*Vm(:,2)'];
hold on;
for i=1:3
    plot(x(i,1),x(i,2),'m+')
end

% Estimate parameters for alignment between MRI to frozen CT using
% principal axis and center of mass
xm=mean(x);
xc=bsxfun(@minus,x,xm);
ym=mean(y);
yc=bsxfun(@minus,y,mean(y));
H=xc'*yc;
[U,D,V]=svd(H);
R=V*diag([1,det(V*U)])*U'
S=sum(diag(xc*R'*yc'))/sum(diag(xc*xc'))
T=ym'-S*R*xm'
yhat=(S*R*x'+repmat(T,1,size(x,1)))';
S_mri_ct=S; R_mri_ct=R; T_mri_ct=T;

% Transform MR to CT using bilinear interpolation
V=head_mri;
[X_mri,Y_mri]=meshgrid(0:(size(head_mri,2)-1),0:(size(head_mri,2)-1));
coords=S_mri_ct*R_mri_ct*[X_mri(:) Y_mri(:)]'+repmat(T_mri_ct,1,numel(X_mri));
X=reshape(coords(1,:),size(X_mri)); Y=reshape(coords(2,:),size(Y_mri));
[Xq,Yq]=meshgrid(0:(size(head_frozen,2)-1),0:(size(head_frozen,1)-1));
Vq=griddata(X,Y,V,Xq,Yq,'linear'); % interp2 doesn't work with non-uniform grid like X and Y here
fusionRGB=cat(3,head_frozen/256,Vq/256,zeros(size(Vq)));
figure; imshow(fusionRGB);

%% Find optimal parameters using Powell's method

V=head_mri;
[X_mri,Y_mri]=meshgrid(0:(size(head_mri,2)-1),0:(size(head_mri,2)-1));
[Xq,Yq]=meshgrid(0:(size(head_frozen,2)-1),0:(size(head_frozen,1)-1));

newT=[0;0]; newR=eye(2); newS=0;
Tscales=10.^(1:-1:-2);
Sscales=10.^(-2:-1:-4);
Rscales=(2*pi/180)*[2 1 0.1];
modes={'tx','ty','s','r'};
maxiter=25;
diff_tx=inf; diff_ty=inf; diff_s=inf; diff_r=inf;
while diff_tx>0.01 || diff_ty >0.01 || diff_s > 0.01 || diff_r > 0.01
    
    init_newT=newT;
    init_newS=newS;
    init_newR=newR;
    
    for nm=1:numel(modes)
        switch modes{nm}
            case {'tx','ty'}
                scales=Tscales;
            case 's'
                scales=Sscales;
            case 'r'
                scales=Rscales;
        end
        
        for ns=1:numel(scales)
            psearch=[-scales(ns) 0 scales(ns)];
            MI=nan(numel(psearch),1);
            niter=0;
            while true                
                for np=find(isnan(MI))'
                    T=T_mri_ct+newT;
                    S=S_mri_ct+newS;
                    R=R_mri_ct*newR;
                    switch modes{nm}
                        case 'tx'
                            T=T+[psearch(np); 0];
                        case 'ty'
                            T=T+[0; psearch(np)];
                        case 's'
                            S=S+psearch(np);
                        case 'r'
                            R=R*[cos(psearch(np)) -sin(psearch(np)); sin(psearch(np)) cos(psearch(np))];
                    end
                    coords=S*R*[X_mri(:) Y_mri(:)]'+repmat(T,1,numel(X_mri));
                    X=reshape(coords(1,:),size(X_mri)); Y=reshape(coords(2,:),size(Y_mri));
                    Vq=griddata(X,Y,V,Xq,Yq,'linear'); % interp2 doesn't work with non-uniform grid like X and Y here
                    jointHistogram = histogram2(double(head_frozen(:)'), double(Vq(:)'), ...
                        [0 256 64; 0 256 64]);
                    jointPdf = jointHistogram / sum(jointHistogram(:));
                    firstMarginalPdf = sum(jointPdf, 1);
                    secondMarginalPdf = sum( jointPdf, 2);
                    MI(np) = log(jointPdf(:) + eps)' * jointPdf(:) - ...
                        log(firstMarginalPdf(:) + eps)' * firstMarginalPdf(:) - ...
                        log(secondMarginalPdf(:) + eps)' * secondMarginalPdf(:);
                end
                
                [~,ind]=max(MI);
                if ind==1
                    psearch=[psearch(1)-scales(ns) psearch];
                    MI=[nan; MI];
                elseif ind==numel(MI)
                    psearch(end+1)=psearch(end)+scales(ns);
                    MI(end+1)=nan;
                else
                    switch modes{nm}
                        case 'tx'
                            newT=[psearch(ind); 0]+newT;
                        case  'ty'
                            newT=[0; psearch(ind)]+newT;
                        case 's'
                            newS=newS+psearch(ind);
                        case 'r'
                            newR=newR*[cos(psearch(ind)) -sin(psearch(ind)); sin(psearch(ind)) cos(psearch(ind))];
                    end
                    break;
                end
                
                niter=niter+1;
                if niter==maxiter, error('Maxiter reached'); end                
            end
        end
    end
    
    diff_tx=abs(init_newT(1)-newT(1));
    diff_ty=abs(init_newT(2)-newT(2));
    diff_s=abs(init_newS-newS);
    diff_r=abs(acos(init_newR(1,1))*180/2*pi-acos(newR(1,1))*180/2*pi);
end

% Transform MR to CT using bilinear interpolation
V=head_mri;
T=T_mri_ct+newT;
S=S_mri_ct+newS;
R=R_mri_ct*newR;
coords=S*R*[X_mri(:) Y_mri(:)]'+repmat(T,1,numel(X_mri));
X=reshape(coords(1,:),size(X_mri)); Y=reshape(coords(2,:),size(Y_mri));
[Xq,Yq]=meshgrid(0:(size(head_frozen,2)-1),0:(size(head_frozen,1)-1));
Vq=griddata(X,Y,V,Xq,Yq,'linear'); % interp2 doesn't work with non-uniform grid like X and Y here
fusionRGB=cat(3,head_frozen/256,Vq/256,zeros(size(Vq)));
figure; imshow(fusionRGB);





