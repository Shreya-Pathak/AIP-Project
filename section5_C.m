clc;
clear;
nv=zeros(100,1);
for i=1:100
    if rem(i,10)==0
        i
    end
N=200;
m=round(N/4);
p = round(log(m));
k=round(m/5);
Phi1=randn(m+p,N);
Phi=Phi1(1:m,:);
Psi = Phi1(m+1:end, :);
% n = sqrt(sum(Phi.^2,1)); 
% Phi= bsxfun(@rdivide,Phi,n);
x0=randn(N,1);
x0(randperm(N,N-k))=0;
an = randn(m+p, 1);
sig = 0.01;
n = sig*an;
y1 = Phi1*x0+n;
y = y1(1:m);
ycv = y1(m+1:end);
[lamfin,xfincv]=get_lambda(y,Phi,ycv,Psi,x0);
x=homotopy_cont(y1,Phi1,lamfin);
nv(i)=norm(x-x0)/norm(x0);
end