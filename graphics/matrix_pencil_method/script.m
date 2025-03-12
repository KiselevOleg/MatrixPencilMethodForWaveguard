clear all;

A=load('dispersion_curve.data');
%A=load('dispersion_curve_2dx_2dL_filter_glass_L=40_full_f_smoothing.data');
%A=load('dispersion_curve_2dx_2dL_filter_glass_L=40_full_f.data');
%A=load('dispersion_curve_2dx_2dL_filter_glass_L=40.data');
plot(A(:,1),A(:,2),'xb');
hold('on');
plot(A(:,1),A(:,3),'.r');

A=load('../distersion_curve_for_K/dispersion_curves.data');
plot(A(:,1),A(:,2),'.k');

xlim([0 3]);
%ylim([0 140/4]);
%xlim([0 1.5]);
ylim([0 35]);
set(gca,'FontSize',24,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

