clear all;

A=load('px1,25py0.data');
time=(1:1:length(A))*16e-9*1e6-53;
plot(time,A/5,'.-c');
hold('on');
##A=load('px1py0.data');
##time=(1:1:length(A))*16e-9*1e6-53;
##plot(time,A/5,'.-c');
##hold('on');
##A=load('px1,5py0.data');
##time=(1:1:length(A))*16e-9*1e6-53;
##plot(time,A/5,'.-c');
##hold('on');
##A=load('px1,75py0.data');
##time=(1:1:length(A))*16e-9*1e6-53;
##plot(time,A/5,'.-c');
##hold('on');
##A=load('px2py0.data');
##time=(1:1:length(A))*16e-9*1e6-53;
##plot(time,A/5,'.-c');
##hold('on');




##A=load('u_smoothing           2           2 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##plot(time,A,'.-m');
##hold('on');

##A=load('u_smoothing           1           1 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##plot(time,A,'.-m');
##hold('on');
##A=load('u_smoothing           3           3 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##plot(time,A,'.-m');
##hold('on');
##A=load('u_smoothing           4           4 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##plot(time,A,'.-m');
##hold('on');
##A=load('u_smoothing           5           5 _600.data');
##time=(1:1:length(A))*16e-9*1e6*8-7;
##plot(time,A,'.-m');
##hold('on');

A=load('u_smoothing          10          10 _600.data');
A=load('u_smoothing           2           2 _600.data');
time=(1:1:length(A))*16e-9*1e6*8-7;
time=load('t_smoothing.data');
time=time(2:length(A)+1)*1e6;
plot(time-7,A,'.-m');
hold('on');



A=load('signal.data');
plot(A(:,1)-2,A(:,2)*0.89,'.-b');

%A=load('signal3.data');
%plot(A(:,1),A(:,2),'x');

xlim([-20 50]);
ylim([-4 6]);

