clear all;

%A=load('xi=1_signal.data');
%A=load('xi=Nx_div_2_signal.data');
%A=load('xi=Nx_signal.data');
%A=load('xi=1_spectrum.data');
A=load('xi=Nx_div_2_spectrum.data');
%A=load('xi=Nx_spectrum.data');

plot(A(:,1),A(:,2),'.-');

