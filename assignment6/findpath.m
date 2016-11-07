function path = findpath(image)
    % find the acumulator image
    aim = zeros(size(image));
    aim(1,:) = image(1,:);
    bim = zeros(size(image));
    for i=2:size(image,1)
        for j=1:size(image,2)
            if j == 1
                [val, ind] = min([image(i-1,j),image(i-1,j+1)]);
                cols = [j, j+1];
                aim(i,j) = image(i,j) + val;
                bim(i,j) = cols(ind);
            elseif j == size(image,2)
                [val, ind] = min([image(i-1,j-1),image(i-1,j)]);
                cols = [j-1,j];
                aim(i,j) = image(i,j) + val;
                bim(i,j) = cols(ind);
            else
                [val, ind] = min([image(i-1,j-1),image(i-1,j),image(i-1,j+1)]);
                cols = [j-1,j,j+1];
                aim(i,j) = image(i,j) + val;
                bim(i,j) = cols(ind);
            end
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