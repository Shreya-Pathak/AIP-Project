clc;
clear;
I = double(imread('barbara.png'));
I=imresize(I,[256,256]);
% n = (2*randn(size(I)));
% nImage = double(I)+n;
% % imshow(nImage, []);
phi = randn(30,64);
n23 = sqrt(sum(phi.^2,1)); 
phi= bsxfun(@rdivide,phi,n23);
bas = kron(dctmtx(8), dctmtx(8));
H = size(I,1);
W = size(I,2);
fimg=zeros(H,W);
counter=zeros(H,W);
A = phi*bas;
mlist = [17,19,21,23,25];
m=30;
sig=0.01;
k=15;
% [a, c] = max(eigs(A'*A));
% alpha = a;
% % X=bas*reshape(nImage(1:8, 1:8), 64, 1);
% Y=ISTA(A, phi*reshape(nImage(1:8, 1:8), 64, 1), 1, alpha);
for i=[1:1:H-7]
    i
    for j=[1:1:W-7]
        %j
%         dmat=[];
%         dp=c(i:i+7,t*W+j:t*W+j+7);
%         dmat=[dmat diag(reshape(dp',64,1))*bas'];
        an = randn(m, 1)*sqrt(1/m);
        n = sig*an;
        %y = phi*reshape(I(i:i+7, j:j+7),64 , 1) + n;
        y=phi*reshape(dct2(I(i:i+7, j:j+7)),64,1)+n;
%         bas*reshape(nImage(i:i+7, j:j+7), 64, 1)
%         dct2(reshape(nImage(i:i+7, j:j+7), 64, 1))
        %Theta = zeros(64,1);
        [Theta, ~] = num_meas(y, phi , sig, mlist, 0.02,size(mlist,2),55);
       % [Theta, ~] = OMP(phi, y, 0.01);
%         Z = bas*reshape(nImage(i:i+7, j:j+7), 64, 1);
%         [m1,m2]=OMP(dmat,reshape(coded(i:i+7,j:j+7)',64,1),eps);
        counter(i:i+7,j:j+7)=counter(i:i+7,j:j+7)+1;
        imt=idct2(reshape(Theta,8,8));
        %fimg(i:i+7,j:j+7)=fimg(i:i+7,j:j+7)+reshape(imt,8,8);
        fimg(i:i+7,j:j+7)=fimg(i:i+7,j:j+7)+imt;
    end
end
imshow(fimg./counter,[]);
error1 = sqrt(sum((I-fimg./counter).^2)/sum(I.^2));