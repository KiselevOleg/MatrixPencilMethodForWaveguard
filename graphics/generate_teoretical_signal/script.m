clear all;

A=load('signal.data');
plot(A(:,1)/2/pi,sqrt(A(:,2).^2+A(:,3).^2)/500,'.-');

hold('on');

A=load('specrum_sin-1mks-load.data');
A=load('specrum_sin-1mks-20mm-3000-12-44.data');
plot(A(:,1),A(:,2),'.-');

hold('on');

A=load('specrum_sin-1mks-load.data');
plot(A(:,1),A(:,2)/100,'.-');

xlim([0 3]);

