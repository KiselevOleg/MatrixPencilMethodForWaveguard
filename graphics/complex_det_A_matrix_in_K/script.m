clear all;

S=load('sizes.data');
ReSize=S(1); ImSize=S(2);

A=load('AY.data');
X=reshape(A(:,1),ImSize,ReSize);
Y=reshape(A(:,2),ImSize,ReSize);
Z=reshape(A(:,3),ImSize,ReSize);

surf(X,Y,Z);

