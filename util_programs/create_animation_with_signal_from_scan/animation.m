clear all;

scan_file_name="scan.mat";

min_animation_image_oxy_value=-30;
max_animation_image_oxy_value=30;

animation_window_width=640*2;
animation_window_height=480*2;

get_point_for_t_from=3000;
get_point_for_t_to=8000;
get_only_Nth_point_for_t_where_N=25;



figure(1, "position",[200,200,animation_window_width,animation_window_height]);

load(scan_file_name);

coord_=coord_*1e3;
time_=time_*1e6;

N=get_only_Nth_point_for_t_where_N;

data_=data(1+1:length(coord_),get_point_for_t_from:get_point_for_t_to);
coord_=coord_(1+1:length(coord_));
time_=time_(get_point_for_t_from:get_point_for_t_to);

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

tind=1:1:length(time_);
for tj=tind
  plot(coord_,data(:,tj),'.-');
  title(["t=" num2str(time_(tj))]);
  xlabel("x");
  ylim([min_animation_image_oxy_value max_animation_image_oxy_value]);
  xlim([min(coord_) max(coord_)]);
  set(gca,"FontSize",24,"fontWeight","bold");
  set(findall(gcf,"type","text"),"FontSize",24,"fontWeight","bold");

  pause(0.1);

  drawnow
  frame = getframe(1);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im);
  if tj == 1
    imwrite(imind,cm,filename_signal,"gif", "Loopcount",inf, "DelayTime", .2);
  elseif tj == tind(length(tind))
    imwrite(imind,cm,filename_signal,"gif","WriteMode","append", "DelayTime", .2);
  else
    imwrite(imind,cm,filename_signal,"gif","WriteMode","append", "DelayTime", .2);
  end
end

xind=1:1:length(coord_);
for xi=xind
  plot(time_,data(xi,:),'.-');
  title(["x=" num2str(coord_(xi))]);
  xlabel("t");
  ylim([min_animation_image_oxy_value max_animation_image_oxy_value]);
  xlim([min(time_) max(time_)]);
  set(gca,"FontSize",24,"fontWeight","bold");
  set(findall(gcf,"type","text"),"FontSize",24,"fontWeight","bold");

  pause(0.1);

  drawnow
  frame = getframe(1);
  im = frame2im(frame);
  [imind,cm] = rgb2ind(im);
  if tj == 1
    imwrite(imind,cm,filename_signal_viewed_by_x,"gif", "Loopcount",inf, "DelayTime", .5);
  elseif tj == tind(length(tind))
    imwrite(imind,cm,filename_signal_viewed_by_x,"gif","WriteMode","append", "DelayTime", .5);
  else
    imwrite(imind,cm,filename_signal_viewed_by_x,"gif","WriteMode","append", "DelayTime", .5);
  end
end

