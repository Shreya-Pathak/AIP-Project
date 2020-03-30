N=1000;
k=50;
sig = sqrt(0.1);
m=360;
cv = [10 20 30 40 50 60 70 80];
i = 1;
rmse = zeros(8,1);
for mcv=cv
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
    rmse(i) = sqrt(sum((x_est-x).^2))/sqrt(sum(x.^2));
    i = i+1;
end
rmse
plot(cv,10*log(rmse));