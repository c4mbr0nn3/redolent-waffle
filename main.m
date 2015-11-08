% Homework Numerical Methods (Ferronato)
%
% full matrix load
prompt='File matrice: ';
str=input(prompt,'s');
A=csvread(str);
% convert matrix to CRS
prompt2='Nome file formato CRS: ';
filename=input(prompt2,'s');
mat2crs(A, filename) %credit to: Roger B. Sidje (rbs@maths.uq.edu.au)
% CRS matrix to vector
A_crs=csvread(filename);
n=A_crs(1,1);
nz=A_crs(1,2);
IA=A_crs([2:2+n],1);
JA=A_crs([3+n:2+n+nz],1);
SYSMAT=A_crs([3+n:2+n+nz],2);
% clear workspace
clear A prompt prompt2 str filename A_crs
% matrix by vector product
