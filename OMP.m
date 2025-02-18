function [Theta, T] = OMP(A, y, eps)
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = zeros(n,1);
    i = 0;
%     B =  sqrt(sum(A.^2,1));
%     Au =  bsxfun(@rdivide,A,B);
    while sum(r.^2) > eps
        a = A'*r;
        [M,I] = max(abs(a));
        T1 = union(T1, I(1));
        %T1;
        i = i+1;
        Theta1(T1) = pinv(A(:,T1))*y;
        r = y - A*Theta1;
        sum(r.^2);
    end
    Theta = Theta1;
    T = T1;
end
