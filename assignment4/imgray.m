function imgray(img)

img=(img-min(img(:)))/(max(img(:))-min(img(:)));
imshow(img);
