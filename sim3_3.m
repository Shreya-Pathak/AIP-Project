clc;
clear;
N=1000;
k=50;
sigv = [sqrt(0.02) sqrt(0.04) sqrt(0.06) sqrt(0.08) sqrt(0.1) sqrt(0.12) sqrt(0.14) sqrt(0.16) sqrt(0.18) sqrt(0.2)];
M=400;
m=352;
mcv = 48;
%mv = zeros(10,1);
i = 1;
rmsecv = zeros(10,1);
rmsew = zeros(10,1);
rmse = zeros(10,1);
for sig=sigv
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
    x2 = OMPw(Au, y, M);
    rmsecv(i) = sqrt(sum((x_est-x).^2))/sqrt(sum(x.^2));
    rmse(i) = sqrt(sum((x1-x).^2))/sqrt(sum(x.^2));
    rmsew(i) = sqrt(sum((x2-x).^2))/sqrt(sum(x.^2));
    i = i+1;
end
rmsecv'
rmse'
rmsew'
mv = sigv.^2;
plot(mv,10*log(rmsecv), 'Marker',"*");
hold on;
plot(mv, 10*log(rmse), 'Marker','+');
hold on;
plot(mv, 10*log(rmsew), 'Marker','>');
hold off;