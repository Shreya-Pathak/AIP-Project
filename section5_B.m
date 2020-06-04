clc;
clear;
N=1000;
m=400;
k=50;
d=150;
sig=0.01;
A=randn(m,N)/sqrt(m);
n = sqrt(sum(A.^2,1)); 
A= bsxfun(@rdivide,A,n);
%Acv=randn(mcv,N)/sqrt(m);
%n = sqrt(sum(Acv.^2,1)); 
%Acv= bsxfun(@rdivide,Acv,n);
%Acv=Acv*sqrt(mcv/m);
x=randn(N,1);
x(randperm(N,N-k))=0;
an = randn(m, 1)*sqrt(1/m);
%ancv = randn(mcv,1)*sqrt(1/m);
n = sig*an;
%ncv = sig*ancv;
y = A*x+n;
%ycv = Acv*x+ncv;
mlist=[40:40:320];
%x2=OMP_ward(A,y,0.01,k+50);
[x_out, j] = num_meas(y, A, sig, mlist, 0.03,size(mlist,2),k);
figure(1)
plot(x)
title('Original signal')
xlabel('time')
ylabel('value')
figure(2)
plot(x_out)
title('Reconstructed signal')
xlabel('time')
ylabel('value')