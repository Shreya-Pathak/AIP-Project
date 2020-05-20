function [xcap,cas,xo,o,ecvp,ecvo,egp,ego] = OMPCV_fort4(A, Acv, y, ycv, d,xtrue,sig)
    % ego = epsilon generalised o (oracle)
    % egp = epsilon generalised p
    % cas is 1 if o==p (dicard for theorem 4) , and see below for other
    % cases
    epscv = sum(ycv.^2);
    r = y;
    T1 = [];
    n = size(A,2);
    Theta1 = zeros(n,1);
%     B =  sqrt(sum(A.^2,1));
%     Au =  bsxfun(@rdivide,A,B);
    minepscv = epscv;
    minx = Theta1;
    ego=sum(xtrue.^2)+sig^2;
    xo=zeros(size(xtrue));
    oidx=0;
    p=0;
    for i = 1:d        
        a = A'*r;
        [~,I] = max(abs(a));
        T1 = union(T1, I(1));
        Theta1(T1) = pinv(A(:,T1))*y;
        r = y - A*Theta1;
        epscv = sum((Acv*Theta1 - ycv).^2);
        if epscv < minepscv
            minepscv = epscv;
            minx = Theta1;
            p = i;
            if any(xtrue~=0 & Theta1==0) % T\Tp isnt phi 
                cas=2;
            else
                cas=3;
            end
            %ecvp=sum((Acv*Theta1 - ycv).^2);
            egp=sum((Theta1 - xtrue).^2)+sig^2;
        end
%         if i==p
%             xp=Theta1;
%             if any(xtrue~=0 & Theta1==0) % T\Tp isnt phi 
%                 cas=2;
%             else
%                 cas=3;
%             end
%             ecvp=sum((Acv*Theta1 - ycv).^2);
%             egp=sum((Theta1 - xtrue).^2)+sig^2;
%         end
        %sum((Theta1 - xtrue).^2)+sig^2;
        if sum((Theta1 - xtrue).^2)+sig^2<ego
            ego=sum((Theta1 - xtrue).^2)+sig^2;
            ecvo=sum((Acv*Theta1 - ycv).^2);
            oidx=i;
            xo=Theta1;
        end
    end
    xcap = minx;
    if oidx==p
        cas=1;
    end
    o=oidx;
    ecvp = minepscv;
    %size(setdiff(find(xtrue~=0),find( minx==0)))
end