function [xcap,netaomp,netaor,netacv] = OMPCV_VI(A, Acv, y, ycv, d,xtrue)
    epscv = sum(ycv.^2);
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = zeros(n,1);
%     B =  sqrt(sum(A.^2,1));
%     Au =  bsxfun(@rdivide,A,B);
    minepscv = epscv;
    minx = Theta1;
    netaor=norm(xtrue);
    for i = 1:d
        a = A'*r;
        [M,I] = max(abs(a));
        T1 = union(T1, I(1));
        Theta1(T1) = pinv(A(:,T1))*y;
        r = y - A*Theta1;
        epscv = sum((Acv*Theta1 - ycv).^2);
        netaor=min(netaor,norm(xtrue-Theta1));
        if epscv < minepscv
            minepscv = epscv;
            minx = Theta1;
        end
    end
    xcap = minx;
    netaomp=norm(xtrue-Theta1);
    netacv=sqrt(minepscv);
end