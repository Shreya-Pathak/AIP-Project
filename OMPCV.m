function [xcap] = OMPCV(A, Acv, y, ycv, d)
    epscv = sum(ycv.^2);
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = zeros(n,1);
%     B =  sqrt(sum(A.^2,1));
%     Au =  bsxfun(@rdivide,A,B);
    minepscv = epscv;
    minx = Theta1;
    for i = 1:d
        a = A'*r;
        [M,I] = max(abs(a));
        T1 = union(T1, I(1));
        Theta1(T1) = pinv(A(:,T1))*y;
        r = y - A(:,T1)*Theta1(T1);
        epscv = sum((Acv*Theta1 - ycv).^2);
        if epscv < minepscv
            minepscv = epscv;
            minx = Theta1;
        end
    end
    xcap = minx;
end