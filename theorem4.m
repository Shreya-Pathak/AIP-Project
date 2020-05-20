clear;
N=1000;
m=400;
k=50;
d=150;
sig=0.01;
lamnot=4;
beta5=0.0376;
for mcv=[48:16:48]
    C0=beta5*(lamnot^2/(mcv-2*lamnot^2));
    C1=2*C0+1+2*sqrt(C0^2+C0);
    gc=0; % greater count ecvp/ecvo>c1
    totc3=0;
    cvec = zeros(1000, 1);
    for numexpt=1:1000
        if rem(numexpt,100)==0
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
    [xcap,cas,xo,o,ecvp,ecvo,egp,ego]=OMPCV_fort4(A, Acv, y, ycv, d,x,sig);
    if cas==3 && ~any(x~=0 & xo==0)
        totc3=totc3+1;
        cvec(totc3) = egp/ego;      
    end
    end
    c = cvec(1:totc3, :);
    size(find(c<=C1))
end