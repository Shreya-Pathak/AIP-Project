clc;
clear;
N=1000;
k=50;
sig = sqrt(0.1);
M=400;
cv = [100 90 80 70 60 50 40 30 20 10];
mv = zeros(10,1);
i = 1;
rmsecv = zeros(10,1);
rmse = zeros(10,1);
for mcv=cv
    m = M-mcv;
    mv(i) = m;
    x = zeros(N,1);
    indices = randperm(N); x(indices(1:k)) = 10*randn(k,1);
    A = randn(m,N)*sqrt(1/m);
    Au = columnNormalise(A);
    Acv = randn(mcv, N)*sqrt(1/m);
    Aucv = columnNormalise(Acv).*sqrt(mcv/m);
    an = randn(m, 1)*sqrt(1/m);
    ancv = randn(mcv,1)*sqrt(1/m);
    n = sig*an;
    ncv = sig*ancv;
    y = Au*x+n;
    ycv = Aucv*x+ncv;
    x_est = OMPCV(Au, Aucv, y, ycv, 100);
    x1 = OMP(Au, y, sig^2); 
    rmsecv(i) = sqrt(sum((x_est-x).^2))/sqrt(sum(x.^2));
    rmse(i) = sqrt(sum((x1-x).^2))/sqrt(sum(x.^2));
    i = i+1;
end
Ar = randn(M,N)*sqrt(1/M);
Aur = columnNormalise(Ar);
anr = randn(M,1)*sqrt(1/M);
nr = sig*anr;
yr = Aur*x+nr;
x1 = OMP(Aur, yr, sig^2);
er = sqrt(sum((x1-x).^2))/sqrt(sum(x.^2));
rmsecv'
rmse'
er
plot(mv,10*log(rmsecv), 'Marker',"*");
hold on;
plot(mv, 10*log(rmse), 'Marker','+');
hold on;
plot(mv, 10*log(er), 'Marker','>');
hold off;