function [x_out, j] = num_meas(y, Phi, gamma, mlist, tau,p,k)
    m = size(Phi,1);
    N = size(Phi, 2);
    m1 = mlist(1);
    Phi1 = Phi(1:m1, :);
    y1 = y(1:m1);
    j=p+1;
    %x_prev = OMP_ward(Phi1, y1, gamma,k+50);
    %x=x_prev;
    for i=1:p
        i
        %x_prev = x;
        m1 = mlist(i);
        m_rem = m-m1;
        Phi1 = Phi(1:m1, :);
        y1 = y(1:m1);
        Psi = Phi(m1+1:end, :);
        x = OMP_ward(Phi1, y1, gamma,k+50);
        y2 = y(m1+1:end);
        res = y2-Psi*x;
        rat = sqrt(m_rem)*(norm(res))/(norm(y)*(sqrt(m_rem)-3*log(p)));
        rat
        if(rat<=tau) 
            j=i;
            x_out = x;
            break;
        end
    end
    if(j==p+1)
        x_out = OMP_ward(Phi,y,gamma,k+50);
    end
end