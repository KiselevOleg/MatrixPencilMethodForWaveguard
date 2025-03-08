clear all;

A=load('dispersion_curves.data');
plot(A(:,1),A(:,2),'.b');
hold('on')
plot(A(:,1),A(:,3),'.r');

ylim([-5 40]);
xlim([0 2]);

xlabel("f MGh");

set(gca,'FontSize',24,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

