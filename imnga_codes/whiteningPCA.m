function z = whiteningPCA(data, dnum)
cc = (data'*data)./size(data,1);
[u, s] = svd(cc);
% s = diag(s);

[J,M] = size(data);
eig_val = diag(s);
[eig_val,idx] = sort(eig_val);
s = s(idx,idx);
u = u(:,idx);
s = s(M-dnum:M-1,M-dnum:M-1);
u = u(:,M-dnum:M-1);
s = diag(s);

w = u*diag(sqrt(1./(s+1e-6)));
z = data*w;
end
