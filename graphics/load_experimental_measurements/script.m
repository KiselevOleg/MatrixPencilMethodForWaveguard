clear all;

%A=load('xi=1_signal.data');
A=load('xi=Nx_div_2_signal.data');
%A=load('xi=Nx_signal.data');
%A=load('xi=1_spectrum.data');
%A=load('xi=Nx_div_2_spectrum.data');
%A=load('xi=Nx_spectrum.data');

plot(A(:,1)/2.8,A(:,2),'.-');

B=load('u_smoothing.data');
T=load('t_smoothing.data');
hold('on');
%plot(T(2:length(T))*1000000,B(:)/40,'.-');
%plot(A(1:length(B),1)*2,B(:,1)/20,'.-');

