x=randn(1000,1);
x(randperm(1000,950))=0;
data=[];
sigman=0.1;
mcv=50;
m=150;
N=1000;
A=randn(m,N)/sqrt(m);
n = sqrt(sum(A.^2,1)); 
A= bsxfun(@rdivide,A,n);
Acv=randn(mcv,N)/sqrt(m);
n = sqrt(sum(Acv.^2,1)); 
Acv= bsxfun(@rdivide,Acv,n);
Acv=Acv*sqrt(mcv/m);
xcap=OMPCV(A,Acv,A*x+sigman*randn(m,1)/sqrt(m),Acv*x+sigman*randn(mcv,1)/sqrt(m),70);
for i =[0:1:1000]
    i
    Acv=randn(mcv,N)/sqrt(m);
    n = sqrt(sum(Acv.^2,1)); 
    Acv= bsxfun(@rdivide,Acv,n);
    Acv=Acv*sqrt(mcv/m);
    data=[data; (norm(Acv*x-Acv*xcap+sigman*randn(mcv,1)/sqrt(m)))^2];
end
pde=fitdist(data,'Normal');
ex=(norm(x-xcap))^2;
mu=(ex+sigman^2)*mcv/m;
sigma=(ex+sigman^2)*sqrt(2*mcv)/m;
gaus=@(x) exp(-(x-mu)^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
xvals=[mu-3*sigma:6*sigma/100:mu+3*sigma]';
pd=pdf(pde,xvals);
theor=arrayfun(gaus,xvals);
plot(xvals,pd,xvals,theor);
legend({'experimental','theoritical'});
