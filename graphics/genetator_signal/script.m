clear all;

%A=load("generator_signal_with_load-2025-03-10.data");
f=load("generator_signal_without_load-2025-03-10.data");
time=(1:1:length(f))*16e-9*1e6;
%plot(time,f,'.-');

omega=0.01:0.025:30;
F=omega;
for omegai=1:1:length(omega)
  F(omegai)=0;
  for ti=1:1:length(time)
    F(omegai)=F(omegai)+f(ti)*exp(1i*omega(omegai)*time(ti));
  end
end

plot(omega,abs(F),'.-');
hold('on');

Ft=2*sin(0.5*omega)./omega*87;
plot(omega,abs(Ft));

