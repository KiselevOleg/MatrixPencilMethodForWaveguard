clear all;

%A=load('xi=1_signal.data');
%A=load('xi=Nx_div_2_signal.data');
%A=load('xi=Nx_signal.data');
%A=load('xi=1_spectrum.data');
A=load('xi=Nx_div_2_spectrum.data');
%A=load('xi=Nx_spectrum.data');

plot(A(:,1)/2.8,A(:,2),'.-');

B=load('u_smoothing.data');
T=load('t_smoothing.data');
%hold('on');
%plot(T(2:length(T))*1000000,B(:)/40,'.-');
%plot(A(1:length(B),1)*2,B(:,1)/20,'.-');

##A=load('xi=Nx_div_2_spectrum_1.data');
##plot(A(:,1)/2.8,A(:,2),'.-');
##hold('on');
##A=load('xi=Nx_div_2_spectrum_10.data');
##plot(A(:,1)/2.8,A(:,2),'.-');
##A=load('xi=Nx_div_2_spectrum_2000.data');
##plot(A(:,1)/2.8,A(:,2),'.-');


##A=load('xi=Nx_div_2_signal_1.data');
##plot(A(:,1)/2.8,A(:,2),'.-');
##hold('on');
##A=load('xi=Nx_div_2_signal_10.data');
##plot(A(:,1)/2.8,A(:,2),'.-');
##A=load('xi=Nx_div_2_signal_2000.data');
##plot(A(:,1)/2.8,A(:,2),'.-');

##B=load('xi=Nx_div_2_signal_1.data');
##B=B((1:1:1026)*2-1,1);
##A=load('xi=Nx_div_2_signal_smoothing_1.data');
##plot(B(:,1)/2.8,A(:,1),'.-');
##hold('on');
##A=load('xi=Nx_div_2_signal_smoothing_10.data');
##plot(B(:,1)/2.8,A(:,1),'.-');
##A=load('xi=Nx_div_2_signal_smoothing_2000.data');
##plot(B(:,1)/2.8,A(:,1),'.-');

set(gca,'FontSize',24,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

