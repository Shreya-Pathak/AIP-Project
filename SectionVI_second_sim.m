clear;
N=3600;
d=100;
k=200;
C=1;
%r=70;
m=800;
sig=sqrt(0.05);
rmin=5;
rmax=90;
total_expt=10;
avgor=zeros(rmax-rmin+1,1);
avgcv=zeros(rmax-rmin+1,1);
avgomp=zeros(rmax-rmin+1,1);
for r=[rmin:rmax]
    r
    tor=0;
    tcv=0;
    tomp=0;
    for ne=1:total_expt
        n=m-r;
        eps=3/sqrt(r);
        x0=zeros(N,1);
        x0(randperm(N,d))=1;
        xa=x0;
        %xa=xa/sqrt(sum(xa.^2));
        A=randn(n,N);
        n1 = sqrt(sum(A.^2,1));
        A= bsxfun(@rdivide,A,n1);
        Acv=randn(r,N)/sqrt(r);
        y=A*xa + randn(n,1)*sig;
        ycv=Acv*xa;
        [xcap,netaomp,netaor,netacv] = OMPCV_VI(A, Acv, y, ycv, k,x0);
        tor=tor+netaor;
        tcv=tcv+netacv;
        tomp=tomp+netaomp;
    end
    avgor(r-rmin+1)=tor/total_expt;
    avgcv(r-rmin+1)=tcv/total_expt;
    avgomp(r-rmin+1)=tomp/total_expt;
end
