clear all;

A=load('dispersion_curve.data');
plot(A(:,1),A(:,2),'.b');
hold('on');
plot(A(:,1),A(:,3),'.r');

