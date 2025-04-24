clear all;

S=load('wavelet_transform_sizes.data');
M=load('wavelet_transform.data');

D=reshape(M(:,3),S(2),S(1));
A=reshape(M(:,1),S(2),S(1));
B=reshape(M(:,2),S(2),S(1));

surf(A,B,D,'EdgeColor','none');
%surf(A,B,D);

view(2);

