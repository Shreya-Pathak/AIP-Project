function [xcap_seq,netaor] = OMP_for_secVI(A,y,d,xtrue)
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = zeros(n,1);  
    netaor=norm(xtrue);
    xcap_seq=zeros(d+1,size(xtrue,1));
    for i = 1:d
        a = A'*r;
        [~,I] = max(abs(a));
        T1 = union(T1, I(1));
        Theta1(T1) = pinv(A(:,T1))*y;
        netaor=min(netaor,norm(xtrue-Theta1));
        r = y - A*Theta1;
        xcap_seq(i+1,:)=Theta1;
    end
end