function y = re(x,xp)
y = norm(x-xp,'fro')/norm(xp,'fro');
end