clear all;

x0 = 0;%mm
x1 = 75;%mm
dx = 0.25;%mm
Nx = (x1-x0)/dx+1;

dt = 16e-9;%s
Nt = 15000;
Nt_before_signal=3000;



data = zeros(Nx,Nt);
data_normalized = zeros(Nx,Nt);

x_index=1;
for x = x0:dx:x1+1e-10
     signal = load(strcat(
         'px',strrep(num2str(x),'.',','),
         'py0.txt'
     ));

    data(x_index,:) = signal;
    data_normalized(x_index,:) = signal/norm(signal);
    x_index = x_index+1;
end

% surf(data,'lineStyle','none');
% view([0,90]);
surf(data_normalized,'lineStyle','none');
view([0,90]);
ylim([1 Nx]); xlim([1 Nt]);

coord_ = linspace(x0,x1,Nx)/1000;
time_ = linspace(-Nt_before_signal*dt,(Nt-Nt_before_signal)*dt,Nt);
save scan.mat coord_ time_ data;

