function p = psnr(I, Ir)

if max(abs(Ir(:)))>100
    MaxPower = 255.^2;
else
    MaxPower = 1;
end

I = double(I);
Ir = double(Ir);
mse = mean((I(:) - Ir(:)).^2);
p = 10*log10(MaxPower/mse);

end
