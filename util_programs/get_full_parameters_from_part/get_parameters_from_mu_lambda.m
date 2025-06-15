clear all;

rho=1;
mu=1;
lambda=2;

E=lambda*0.5/(lambda+mu);
nu=sqrt(lambda*0.5/(lambda+mu));

Cp=sqrt((lambda+mu+mu)/rho);
Cs=sqrt(mu/rho);

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

