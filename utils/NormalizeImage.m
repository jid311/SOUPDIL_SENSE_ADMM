function DATA = NormalizeImage(DATA)  %图像归一化
im = ifft2c(DATA);
scale = max(abs(im(:)));
DATA = DATA./scale;
end