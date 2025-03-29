clear all;

load("scan.mat");

width=640*2;
height=480*2;
figure(1, "position",[200,200,width,height]);

coord_=coord_*1e3;
time_=time_*1e6;



N=4;

data_=data(5:length(coord_),2900:length(time_)/1.35);
coord_=coord_(5:length(coord_));
time_=time_(2900:length(time_)/1.35);

coord=coord_;
time=[];
data=[];

disp("start");
if(N==1)
  time=time_;
  data=data_;
else
  i_next=1; di=1000;
  for i=1:1:length(time_)/N
    i_=1+(i-1)*N;
    time(i)=time_(i_);

    if(i_>i_next)
      disp(i_);
      i_next=i_next+di;
    end

    for j=1:1:length(coord_)
      v=0;
      for k=1:1:N
        v=v+data_(j,i_+k-1);
      end
      v=v/N;
      data(j,i)=v;
    end
  end
end



coord_=coord;
time_=time;



filename_signal = "signal.gif";
filename_signal_viewed_by_x = "signal_viewed_by_x.gif";

tind=1:25:length(time_);
for tj=tind
  plot(coord_,data(:,tj),'.-');
  title(['t=' num2str(time_(tj))]);
  xlabel('x');
  ylim([-100 100]);
  xlim([min(coord_) max(coord_)]);
  set(gca,'FontSize',24,'fontWeight','bold');
  set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

  pause(0.1);

  drawnow
  frame = getframe(1);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im);
  if tj == 1
    imwrite(imind,cm,filename_signal,'gif', 'Loopcount',inf, "DelayTime", .2);
  elseif tj == tind(length(tind))
    imwrite(imind,cm,filename_signal,'gif','WriteMode','append', "DelayTime", .2);
  else
    imwrite(imind,cm,filename_signal,'gif','WriteMode','append', "DelayTime", .2);
  end
end

xind=1:1:length(coord_);
for xi=xind
  plot(time_(tind),data(xi,tind),'.-');
  title(['x=' num2str(coord_(xi))]);
  xlabel('t');
  ylim([-100 100]);
  xlim([min(time_) max(time_)]);
  set(gca,'FontSize',24,'fontWeight','bold');
  set(findall(gcf,'type','text'),'FontSize',24,'fontWeight','bold');

  pause(0.1);

  drawnow
  frame = getframe(1);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im);
  if tj == 1
    imwrite(imind,cm,filename_signal_viewed_by_x,'gif', 'Loopcount',inf, "DelayTime", .5);
  elseif tj == tind(length(tind))
    imwrite(imind,cm,filename_signal_viewed_by_x,'gif','WriteMode','append', "DelayTime", .5);
  else
    imwrite(imind,cm,filename_signal_viewed_by_x,'gif','WriteMode','append', "DelayTime", .5);
  end
end

