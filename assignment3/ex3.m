load 'visiblehuman'
head_mri=double(head_mri);
head_frozen=double(head_frozen);

% Find center of mass for MRI and CT
mm=com(head_mri);
mf=com(head_frozen);
[Um,Sm]=principal_axis(head_mri,mm);
[Uf,Sf]=principal_axis(head_frozen,mf);

% Plot principal axis
figure;
subplot(1,2,1);
imagesc(head_mri);
hold on;
plot(mm(1),mm(2),'mx');
hold on;
plot([mm(1) sqrt(Sm(1))*Um(1,1)+mm(1)],[mm(2) sqrt(Sm(1))*Um(2,1)+mm(2)],'g-');
hold on;
plot([mm(1) sqrt(Sm(2))*Um(1,2)+mm(1)],[mm(2) sqrt(Sm(2))*Um(2,2)+mm(2)],'m-');
axis('square')

subplot(1,2,2);
imagesc(head_frozen);
hold on;
plot(mf(1),mf(2),'mx');
hold on;
plot([mf(1) sqrt(Sf(1))*Uf(1,1)+mf(1)],[mf(2) sqrt(Sf(1))*Uf(2,1)+mf(2)],'g-');
hold on;
plot([mf(1) sqrt(Sf(2))*Uf(1,2)+mf(1)],[mf(2) sqrt(Sf(2))*Uf(2,2)+mf(2)],'m-');
axis('square')

