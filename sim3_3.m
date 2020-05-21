clc;
clear;
N=1000;
k=50;
sigv = [sqrt(0.02) sqrt(0.04) sqrt(0.06) sqrt(0.08) sqrt(0.1) sqrt(0.12) sqrt(0.14) sqrt(0.16) sqrt(0.18) sqrt(0.2)];
M=400;
m=352;
mcv = 48;
%mv = zeros(10,1);
j = 1;
rmsecv = zeros(10,1);
rmsew = zeros(10,1);
rmse = zeros(10,1);
for sig=sigv
    i = 1;
    avgcv = zeros(100, 1);
    avg = zeros(100, 1);
    avgw = zeros(100, 1);
    for numexp=1:100
        if rem(numexp, 100)==0
            numexp
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
     Ar=randn(M,N)/sqrt(M);
    n = sqrt(sum(Ar.^2,1)); 
    Ar= bsxfun(@rdivide,Ar,n);
anr = randn(M,1)*sqrt(1/M);
nr = sig*anr;
yr = Ar*x+nr;
    x_est = OMPCV(A, Acv, y, ycv, 100);
    x1 = OMP(Ar, yr, sig^2); 
    x2 = OMPw(Ar, yr, M);
    avgcv(i) = sum((x_est-x).^2);
    avg(i) = sum((x1-x).^2);
    avgw(i) = sum((x2-x).^2);
    i = i+1;
    end
    rmsecv(j) = mean(avgcv);
    rmse(j) = mean(avg);
    rmsew(j) = mean(avgw);
    j = j+1;
end
rmsecv'
rmse'
rmsew'
mv = sigv.^2;
plot(mv,10*log10(rmsecv), 'Marker',"*");
hold on;
plot(mv, 10*log10(rmse), 'Marker','+');
hold on;
plot(mv, 10*log10(rmsew), 'Marker','>');
hold off;