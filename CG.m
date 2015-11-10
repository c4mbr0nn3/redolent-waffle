% Conjugate Gradient (CG)
function CG(pvec, vec, sysmat, ia, ja, nequ)
x=pvec;
b=vec;
for i=1:nequ
    Ax(i)=0;
end
for i=1:nequ
    k=ia(i);
    Ax(i)=Ax(i)+sysmat(k)*x(i);
    for k=ia(i)+1:ia(i+1)-1
        j=ja(k);
        Ax(i)=Ax(i)+sysmat(k)*x(j);
        Ax(j)=Ax(j)+sysmat(k)*x(i);
    end
end
r=b-Ax;
p=r;
tau=norm(r);
tol=10e-10;
iter=0;
while tau>tol
    % t=A*p
    for i=1:nequ
        t(i)=0;
    end
    for i=1:nequ
        k=ia(i);
        t(i)=t(i)+sysmat(k)*p(i);
        for k=ia(i)+1:ia(i+1)-1
            j=ja(k);
            t(i)=t(i)+sysmat(k)*p(j);
            t(j)=t(j)+sysmat(k)*p(i);
        end
    end
    a=(r*p.')/(p*t.');
    x=x+(a*p);
    r=r-(a*t);
    b=(-r*t.')/(p*t.');
    p=r+(b*p);
    tau=norm(r);
    iter=iter+1;
end

%display iteration and result
disp(x)
disp(iter)