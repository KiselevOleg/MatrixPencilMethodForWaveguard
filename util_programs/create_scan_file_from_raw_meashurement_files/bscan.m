clear all;

x0 = 0;%mm
x1 = 50;%mm
dx = 0.5;%mm
Nx = floor((x1-x0)/dx)+1;

dt = 16e-9;%s
Nt = 15000;
Nt_before_signal=3000;



data = zeros(Nx,Nt);

x_index=1;
for x = x0:dx:x1+1e-10
    data(x_index,:) = load(strcat(
        'px',strrep(num2str(x),'.',','),
        'py0.txt'
    ));

    x_index = x_index+1;
end

surf(data,'lineStyle','none');
view([0,90]);
xlim([1 Nt]); ylim([1 Nx]);

coord_ = linspace(x0,x1,Nx)/1000;
time_ = linspace(-Nt_before_signal*dt,(Nt-Nt_before_signal)*dt,Nt);
save scan.mat coord_ time_ data;

