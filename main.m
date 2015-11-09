% Homework Numerical Methods (Ferronato)

% ask if is CRS of full matrix

% full matrix load
prompt='File matrice: ';
str=input(prompt,'s');
A=csvread(str);

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

% CRS matrix to vector
A_crs=csvread(filename);
n=A_crs(1,1);
nz=A_crs(1,2);
IA=A_crs(2:2+n,1);
JA=A_crs(3+n:2+n+nz,1);
SYSMAT=A_crs(3+n:2+n+nz,2);

% vector v load
prompt3='File vettore v: ';
str2=input(prompt3,'s');
v=csvread(str2);

% clear workspace
clear A prompt prompt2 prompt3 str str2 filename A_crs

% matrix by vector product
% check if matrix is symmetric or not
if (sym==1)
    for i=1:n
        w(i)=0;
    end
    for i=1:n
        k=IA(i);
        w(i)=w(i)+SYSMAT(k).'*v(i);
        for k=IA(i)+1:IA(i+1)-1
            j=JA(k);
            w(i)=w(i)+SYSMAT(k).'*v(j);
            w(j)=w(j)+SYSMAT(k).'*v(i);
        end
    end
else
    for i=1:n
        w(i)=0;
        for k=IA(i):IA(i+1)-1
            j=JA(k);
            w(i)=w(i)+SYSMAT(k).'*v(j);
        end
    end
end