function path = findpath(image)
    % find the acumulator image
    aim = zeros(size(image));
    aim(1,:) = image(1,:);
    bim = zeros(size(image));
    for i=2:size(image,1)
        for j=1:size(image,2)
            cols=[j-1 j j+1];
            cols=cols(cols>0&cols<=size(image,2));
            [val,ind]=min(aim(i-1,cols));
            aim(i,j)=image(i,j)+val;
            bim(i,j)=cols(ind);
        end
    end

    % find the path in polar coordinates
    [~,ind] = min(aim(end,:));
    path = zeros(size(aim,1),1);
    path(end) = ind;
    for i=size(bim,1):-1:2
        path(i-1) = bim(i,ind);
        ind = bim(i,ind);
    end
end