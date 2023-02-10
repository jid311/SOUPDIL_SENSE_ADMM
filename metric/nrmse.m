function NRmse = nrmse(x, ref)
% 归一化均方根误差(Normalized Root Mean Square Error, NRMSE)
% sqrt( sum(||x-ref||^2) / N ) / (max(ref)-min(ref)).
% 文献来源: p.11 of [1] M. Lustig and J. M. Pauly, "SPIRiT: Iterative
% self-consistent parallel imaging reconstruction from arbitrary k-space,"
% Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 457-471, 2010.
%
% 代码来源: Analytical MRI Phantom\MRIPhantomv0-8\exp\isbi_exp2.m
% Error = m-m_num;
% MSE = Error(:)'*Error(:)/numel(Error);
% NRMSE{b}(l) = sqrt(MSE)/(max(abs(m(:))-min(abs(m(:)))));
%
% 代码来源: Analytical MRI Phantom\MRIPhantomv0-8\SensFitting.m
% residue = data-M*s;
% mse = residue'*residue/numel(residue);
% nrmse = sqrt(mse)/(max(abs(data)-min(abs(data))));
%
% fessler\irt\nufft\greengard 中提供了另一种实现.
% nn = norm(xhat(:) - xtrue(:)) / norm(xtrue(:));

IsReal = isreal(x) && isreal(ref);
mse = mean((abs(ref(:)-x(:))).^2);
NRmse = sqrt(mse) / (max(ref(:))-min(ref(:)));
end