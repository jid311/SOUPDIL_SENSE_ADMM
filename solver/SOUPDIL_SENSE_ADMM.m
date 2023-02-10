function [res, out] = SOUPDIL_SENSE_ADMM(kData, opts)
% Copyright ? 2022 by Yan Su and Jizhong Duan, Kunming University of
% Science and Technology. All Rights Reserved.
% This is implementation for SOUPDIL_SENSE_ADMM from undersampled k-space.

DATA = opts.DATA;
[sx,sy,~] = size(DATA);
IterOG = opts.IterOG;
mask = opts.mask;
maskroi = opts.maskroi;
maps = opts.maps;
x = opts.x0; 
mu1 = opts.mu1; 
mu2 = opts.mu2; 
num2 = opts.num2;
wd = opts.wd;
a = opts.a;
a2 = opts.a2;
n = opts.n;
K = 4*n;
et = [logspace(log10(a), log10(a2),IterOG-10) a2*ones(1,10)];

uZ = zeros(size(maps));
SHS  = SHSmap(maps);
imm = opts.im.*maskroi;

randn('state',0); rand('state',0);

if(K<=n)
    D0 = (kron(dctmtx(sqrt(n)),dctmtx(sqrt(n))))';
    D0 = D0(:,1:K);
else
    D0 = (kron(dctmtx(sqrt(n)),dctmtx(sqrt(n))))';
    for jj=n+1:K
        D0(:,jj) = randn(n,1);
        D0(:,jj) = D0(:,jj)/norm(D0(:,jj),2);
    end
end

D = D0; %initial dictionary

%% initalization
out.snr = []; 
out.time = [];
t0 = tic;
t00 = cputime;

%% Iteration
for iter = 1:IterOG
    xp = x;
    %create image patches including wrap around patches
    Ib = [x x(:,1:(sqrt(n)-1));x(1:(sqrt(n)-1),:) x(1:(sqrt(n)-1),1:(sqrt(n)-1))];
    [TE,idx] = my_im2col(Ib,[sqrt(n),sqrt(n)],wd);
    N2 = size(TE,2); %total number of image patches
    
    if(iter==1)
        W = sparse(zeros(K,N2));
        [rows,cols] = ind2sub(size(Ib)-sqrt(n)+1,idx);
    end
    
    [D,W] = SOUPDIL(TE,D,W,et(iter),num2,0,0,0);
        
    %image update
    IMoutR = zeros(size(Ib));
    IMoutI = zeros(size(Ib));
    bbb = sqrt(n);
    ZZ = (D*W);
    Lt = 10000;
    for jj = 1:Lt:N2
        jumpSize = min(jj+Lt-1,N2);
        block = reshape(real(ZZ(:,jj:jumpSize)),bbb,bbb,jumpSize-jj+1);
		blockc = reshape(imag(ZZ(:,jj:jumpSize)),bbb,bbb,jumpSize-jj+1);
        for ii  = jj:jumpSize
            col = cols(ii); 
			row = rows(ii);
            IMoutR(row:row+bbb-1,col:col+bbb-1) = IMoutR(row:row+bbb-1,col:col+bbb-1)+block(:,:,ii-jj+1);
            IMoutI(row:row+bbb-1,col:col+bbb-1) = IMoutI(row:row+bbb-1,col:col+bbb-1)+blockc(:,:,ii-jj+1);
        end
    end
    
    IMout = IMoutR + (0+1i)*IMoutI;
    
    IMout2 = zeros(sx,sy);
    IMout2(1:sx,1:sy)=IMout(1:sx,1:sy);
    IMout2(1:(sqrt(n)-1),:) = IMout2(1:(sqrt(n)-1),:)+ IMout(sx+1:size(IMout,1),1:sy);
    IMout2(:, 1:(sqrt(n)-1)) = IMout2(:, 1:(sqrt(n)-1)) + IMout(1:sx,sy+1:size(IMout,2));
    IMout2(1:(sqrt(n)-1),1:(sqrt(n)-1)) = IMout2(1:(sqrt(n)-1),1:(sqrt(n)-1))+ IMout(sx+1:size(IMout,1),sy+1:size(IMout,2));
    
    Lb = n*ones(sx,sy);
    
    Z  = ( mu1*(fft2c(Smap(maps, x)) + uZ) + kData) ./ (mask + mu1 + eps);
    Xr = SmapH(maps,ifft2c(Z - uZ));
    I2 = IMout2;
    
    x = ( mu1*Xr + mu2*I2 ) ./ ( mu1*SHS + mu2* Lb + eps );
    uZ = uZ + 1*(fft2c(Smap(maps, x)) - Z);
    
    res = x;
    xm = abs(x).*maskroi;
    out.psnr(iter) = psnr(xm,imm);
    out.snr(iter) = snr(xm,imm);
    out.hfen(iter) = hfen(xm,imm);
    out.nrmse(iter) = nrmse(xm,imm);
    out.ssim(iter) = ssim(xm,imm);
    out.re(iter) = re(abs(x).*maskroi,abs(xp).*maskroi);
       
    out.timet(iter) = toc(t0);
    out.time(iter) = cputime - t00;
end % iteration end

out.opts = opts;
out.Iout = abs(x);
out.I1 = opts.im;

end
