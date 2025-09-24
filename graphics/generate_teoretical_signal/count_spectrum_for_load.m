clear all;

A=load('sin-1mks-load.data');
time=(1:1:length(A))*16e-9*1e6-50;
plot(time,A(:,1),'.-');

%xlim([0 2]);
%ylim([-0.01 2]);



f=0:0.01:4;
omega=f*2*pi;
FA=omega;
for omega_i=1:1:length(omega)
  omega_=omega(omega_i);
  FA(omega_i)=0;

  for time_i=1:1:length(time)
    time_=time(time_i);

    FA(omega_i)=FA(omega_i)+A(time_i)*exp(1i*omega_*time_)*(time(2)-time(1));
  end

  omega_
end

plot(f,abs(FA),'.-');


r=[f; abs(FA); real(FA); imag(FA)];
r=r';
dlmwrite('specrum_sin-1mks-load.data', r, 'delimiter', ' ');

