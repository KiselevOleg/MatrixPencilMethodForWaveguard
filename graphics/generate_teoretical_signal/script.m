clear all;

A=load('xi=Nx_div_2_spectrum.data');
time=(1:1:length(A))*16e-9*1e6;
plot(A(:,1),A(:,2).*A(:,1)*10,'.-c');
hold('on');

##A=load('u_smoothing           2           2 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##time=load('t_smoothing.data');
##time=time(2:length(A)+1)*1e6;
##plot(time-7,A,'.-m');
##hold('on');



A=load('signal.data');
plot(A(:,1)/2/3.14159,sqrt(A(:,2).^2+A(:,3).^2)*0.6,'.-b');

%A=load('signal3.data');
%plot(A(:,1),A(:,2),'x');

xlim([0 2]);
ylim([-0.01 2]);

