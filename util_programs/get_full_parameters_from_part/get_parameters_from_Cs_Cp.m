clear all;

rho=1;
Cp=1;
Cs=0.3;

lambda=rho*Cp*Cp-2*rho*Cs*Cs;
mu=rho*Cs*Cs;

E=Cs*Cs*rho*2*(1+(2*Cs*Cs-Cp*Cp)/(2*Cs*Cs-2*Cp*Cp));
nu=(2*Cs*Cs-Cp*Cp)/(2*Cs*Cs-2*Cp*Cp);

disp("rho");
disp(rho);

disp("lambda");
disp(lambda);
disp("mu");
disp(mu);

disp("E");
disp(E);
disp("nu");
disp(nu);

disp("Cp");
disp(Cp);
disp("Cs");
disp(Cs);

