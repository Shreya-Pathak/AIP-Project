function [Theta, T] = OMP_ward(A, y, eps, prev,k)
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = prev;
    i = 0;
%     B =  sqrt(sum(A.^2,1));
%     Au =  bsxfun(@rdivide,A,B);
    while (sum(r.^2) > eps ) && i<k
        sum(r.^2)
        a = A'*r;
        [~,I] = max(abs(a));
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