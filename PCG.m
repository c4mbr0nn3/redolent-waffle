% Preconditioned Conjugate Gradient (PCG)
function PCG(x, b, sysmat, ia, ja, nequ, nterm, lmat)
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
vec=r;
pvec=lsolve(nequ,nterm,ia,ja,lmat,vec);
p=pvec;
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
    vec=r;
    pvec=lsolve(nequ,nterm,ia,ja,lmat,vec);
    v=pvec;
    b=(-v*t.')/(p*t.');
    p=v+(b*p);
    tau=norm(r);
    iter=iter+1;
end

%display iteration and result
disp(x)
disp(iter)