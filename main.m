% Homework Numerical Methods (Ferronato)

% full matrix load
prompt='File matrice: ';
str=input(prompt,'s');
A=csvread(str);

% ask if is CRS or full matrix
choice=menu('La matrice � salvata in formato CRS?','Si','No');
if choice==1
    % request SYSMAT and JA columns number
    A_crs=A;
    prompt1='Numero di colonne vettori SYSMAT e JA: [a,b] ';
    B=input(prompt1);
    
    % creating: n, nz, SYSMAT, JA, IA
    % matrix rank
    n=A_crs(1,1);
    
    % non-zero terms number
    nz=A_crs(1,2);
    
    % SYSMAT vector sort, eliding zero terms at the end
    SYSMAT=A_crs(2:ceil(nz/B(1))+1,1:B(1)).';
    SYSMAT=SYSMAT(:);
    SYSMAT=SYSMAT(1:nz,1);
    
    % JA vector sort, eliding zero terms at the end
    JA=A_crs(ceil(nz/B(1))+2:ceil(nz/B(1))+1+ceil(nz/B(2)),1:B(2)).';
    JA=JA(:);
    JA=JA(1:nz,1);
    
    % IA vector sort, eliding zero terms at the end
    IA=A_crs(ceil(nz/B(1))+ceil(nz/B(2))+2:end,1:B(2)).';
    IA=IA(:);
    IA=IA(1:n+1,1);
    
    sym=1;
else
    % check if matrix is symmetric 
    % if true retain only upper part of matrix
    sym=isequal(A,A.');
    if (sym==1)
    A=triu(A);
    end;

    % convert matrix to CRS
    % credit to: Roger B. Sidje (rbs@maths.uq.edu.au)
    prompt2='Nome file CRS: ';
    filename=input(prompt2,'s');
    mat2crs(A, filename)
    
    % CRS matrix to vector:  n, nz, SYSMAT, JA, IA
    A_crs=csvread(filename);
    n=A_crs(1,1);
    nz=A_crs(1,2);
    IA=A_crs(2:2+n,1);
    JA=A_crs(3+n:2+n+nz,1);
    SYSMAT=A_crs(3+n:2+n+nz,2);
end

% vector v load (v unitary)
prompt3='File vettore v: ';
str2=input(prompt3,'s');
v=csvread(str2);

% calculate vector vec=b (Ax=b)
% matrix by vector product
% check if matrix is symmetric or not
if (sym==1)
    for i=1:n
        vec(i)=0;
    end
    for i=1:n
        k=IA(i);
        vec(i)=vec(i)+SYSMAT(k)*v(i);
        for k=IA(i)+1:IA(i+1)-1
            j=JA(k);
            vec(i)=vec(i)+SYSMAT(k)*v(j);
            vec(j)=vec(j)+SYSMAT(k)*v(i);
        end
    end
else
    for i=1:n
        vec(i)=0;
        for k=IA(i):IA(i+1)-1
            j=JA(k);
            vec(i)=vec(i)+SYSMAT(k).'*v(j);
        end
    end
end

% clear useless variables
clear A A_crs B choice i j k prompt prompt1 prompt2 prompt3 str str2 filename sym

% incomplete Cholesky factorization (Kershaw method)
nequ=n;
nterm=nz;
ia=IA.';
ja=JA.';
sysmat=SYSMAT.';
prec=kersh(nequ,nterm,ia,ja,sysmat);

% calculate initial point pvec=x0
lmat=prec;
pvec=lsolve(nequ,nterm,ia,ja,lmat,vec);
clear IA JA n nz SYSMAT prec

% residual correction method
prompt4='Numero iterazioni CRM: ';
maxiter=input(prompt4);
iter1=0;
b=vec;
x=pvec;
while iter1<maxiter
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
    vec=b-Ax;
    pvec=lsolve(nequ,nterm,ia,ja,lmat,vec);
    x=x+pvec;
    iter1=iter1+1;
end

% ask if use CG or PCG
choice=menu('Quale algoritmo vuoi usare?','CG','PCG');
if choice==1
    CG(x, b, sysmat, ia, ja, nequ)
else
    PCG(x, b, sysmat, ia, ja, nequ, nterm, lmat)
end

clear choice