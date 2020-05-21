clc;
clear;
N=1000;
k=50;
sig = sqrt(0.1);
M=400;
cv = [100 90 80 70 60 50 40 30 20 10];
mv = zeros(10,1);
j = 1;
rmsecv = zeros(10,1);
rmse = zeros(10,1);
for mcv=cv
    i=1;
    avg = zeros(300, 1);
    avgcv = zeros(300, 1);
    m = M-mcv;
    mv(j) = m;
    for numexpt = 1:300
        if rem(numexpt, 100)==0
            numexpt
        end
     A=randn(m,N)/sqrt(m);
    n = sqrt(sum(A.^2,1)); 
    A= bsxfun(@rdivide,A,n);
    Acv=randn(mcv,N)/sqrt(m);
    n = sqrt(sum(Acv.^2,1)); 
    Acv= bsxfun(@rdivide,Acv,n);
    Acv=Acv*sqrt(mcv/m);
    x=randn(N,1);
    x(randperm(N,N-k))=0;
    an = randn(m, 1)*sqrt(1/m);
    ancv = randn(mcv,1)*sqrt(1/m);
    n = sig*an;
    ncv = sig*ancv;
    y = A*x+n;
    ycv = Acv*x+ncv;
    x_est = OMPCV(A, Acv, y, ycv, 100);
    x1 = OMP(A, y, sig^2); 
    avgcv(i) = sum((x_est-x).^2);
    avg(i) = sum((x1-x).^2);
    i = i+1;
    end
    rmsecv(j) = mean(avgcv);
    rmse(j) = mean(avg);
    j=j+1;
end
 Ar=randn(M,N)/sqrt(M);
    n = sqrt(sum(Ar.^2,1)); 
    Ar= bsxfun(@rdivide,Ar,n);
anr = randn(M,1)*sqrt(1/M);
nr = sig*anr;
yr = Ar*x+nr;
x1 = OMP(Ar, yr, sig^2);
er = sum((x1-x).^2);
rmsecv'
rmse'
er
plot(mv,10*log10(rmsecv), 'Marker',"*")
hold on;
plot(mv, 10*log10(rmse), 'Marker','+')
hold on;
plot(mv, 10*log10(er), 'Marker','>')
legend('OMP-CV', 'OMP-residual', 'm=400')
hold off;