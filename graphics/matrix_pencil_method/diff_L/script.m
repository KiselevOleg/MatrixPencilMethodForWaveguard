clear all;

L3=load("dispersion_curve_2dx_2L_filter_L=           3 .data");
L5=load("dispersion_curve_2dx_2L_filter_L=           5 .data");
L7=load("dispersion_curve_2dx_2L_filter_L=           7 .data");
L10=load("dispersion_curve_2dx_2L_filter_L=          10 .data");
L15=load("dispersion_curve_2dx_2L_filter_L=          15 .data");
L20=load("dispersion_curve_2dx_2L_filter_L=          20 .data");
L25=load("dispersion_curve_2dx_2L_filter_L=          25 .data");
L30=load("dispersion_curve_2dx_2L_filter_L=          30 .data");
L35=load("dispersion_curve_2dx_2L_filter_L=          35 .data");
L40=load("dispersion_curve_2dx_2L_filter_L=          40 .data");
L50=load("dispersion_curve_2dx_2L_filter_L=          50 .data");
L60=load("dispersion_curve_2dx_2L_filter_L=          60 .data");
L70=load("dispersion_curve_2dx_2L_filter_L=          70 .data");
L80=load("dispersion_curve_2dx_2L_filter_L=          80 .data");
L90=load("dispersion_curve_2dx_2L_filter_L=          90 .data");
L100=load("dispersion_curve_2dx_2L_filter_L=         100 .data");
L110=load("dispersion_curve_2dx_2L_filter_L=         110 .data");
L120=load("dispersion_curve_2dx_2L_filter_L=         120 .data");
L130=load("dispersion_curve_2dx_2L_filter_L=         130 .data");
L140=load("dispersion_curve_2dx_2L_filter_L=         140 .data");
L150=load("dispersion_curve_2dx_2L_filter_L=         150 .data");
L160=load("dispersion_curve_2dx_2L_filter_L=         160 .data");
L180=load("dispersion_curve_2dx_2L_filter_L=         180 .data");
L200=load("dispersion_curve_2dx_2L_filter_L=         200 .data");

width=640*2;
height=480*2;
figure(1, 'position',[200,200,width,height]);

for animation_type=1:1:2

if(animation_type==1) filename = 'dispersion_curves.gif'; end
if(animation_type==2) filename = 'dispersion_curves_with_K_curves.gif'; end
iLast=24;
for i=1:1:iLast
  if(i==1) Ld=L3; L=3; end
  if(i==2) Ld=L5; L=5; end
  if(i==3) Ld=L7; L=7; end
  if(i==4) Ld=L10; L=10; end
  if(i==5) Ld=L15; L=15; end
  if(i==6) Ld=L20; L=20; end
  if(i==7) Ld=L25; L=25; end
  if(i==8) Ld=L30; L=30; end
  if(i==9) Ld=L35; L=35; end
  if(i==10) Ld=L40; L=40; end
  if(i==11) Ld=L50; L=50; end
  if(i==12) Ld=L60; L=60; end
  if(i==13) Ld=L70; L=70; end
  if(i==14) Ld=L80; L=80; end
  if(i==15) Ld=L90; L=90; end
  if(i==16) Ld=L100; L=100; end
  if(i==17) Ld=L110; L=110; end
  if(i==18) Ld=L120; L=120; end
  if(i==19) Ld=L130; L=130; end
  if(i==20) Ld=L140; L=140; end
  if(i==21) Ld=L150; L=150; end
  if(i==22) Ld=L160; L=160; end
  if(i==23) Ld=L180; L=180; end
  if(i==24) Ld=L200; L=200; end

  if(animation_type==1) plot(Ld(:,1),Ld(:,2),'.b'); end
  if(animation_type==2) plot(Ld(:,1),Ld(:,2),'xb'); end
  hold('on');
  plot(Ld(:,1),Ld(:,3),'.r');

  A=load('../../distersion_curve_for_K/dispersion_curves.data');
  if(animation_type==2) plot(A(:,1),A(:,2),'.k'); end
  hold('off');

  title(["L=" int2str(L)]);

  xlim([0 2]);
  ylim([0 140/4]);
  xlabel('\omega MGh');
  %xlim([0 1.5]);
  %ylim([0 35]);
  set(gca,'FontSize',24,'fontWeight','bold');
  set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

  pause(0.1);

  drawnow
  frame = getframe(1);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im);
  if i == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf, "DelayTime", .9);
  elseif i == iLast
    imwrite(imind,cm,filename,'gif','WriteMode','append', "DelayTime", .9);
  else
    imwrite(imind,cm,filename,'gif','WriteMode','append', "DelayTime", .9);
  end
end

end

