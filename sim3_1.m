N=1000;
k=50;
sig = sqrt(0.1);
m=360;
%cv = 5:85;
cvplot = [10 20 30 40 50 60 70 80]
i = 1;
rmse = zeros(8,1);
for mcv=cvplot
    avg = zeros(300, 1);
    j=1;
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
    avg(j) = sum((x_est-x).^2);
    j=j+1;
    end
    rmse(i) = mean(avg);
    i=i+1;
end
rmse
plot(cvplot,10*log10(rmse))
legend('OMP-CV')