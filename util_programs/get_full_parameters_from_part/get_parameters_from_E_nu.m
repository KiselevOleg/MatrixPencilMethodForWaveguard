clear all;

rho=1;
E=1;
nu=0.3;

lambda=nu*E/(1+nu)/(1-2*nu);
mu=E*0.5/(1+nu);

Cp=sqrt(E/rho)*sqrt((1-nu)/(1+nu)/(1-2*nu));
Cs=sqrt(E/rho)*sqrt(0.5/(1+nu));

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

