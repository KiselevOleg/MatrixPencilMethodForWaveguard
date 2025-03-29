clear all;

N=4;
%load('Glass_rect_1mks-2025-03-19.mat');
load('Glass_rect_1mks-2025-03-13_small.mat');

data_=data(20:length(coord_),50:length(time_)/2.5);
data2_=data2(20:length(coord_),50:length(time_)/2.5);
coord_=coord_(20:length(coord_));
time_=time_(50:length(time_)/2.5)/1000;

coord=coord_;
time=[];
data=[];
data2=[];

disp("start")
if(N==1)
  time=time_;
  data=data_;
  data2=data2_;
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
      v2=0;
      for k=1:1:N
        v=v+data_(j,i_+k-1);
        v2=v2+data2_(j,i_+k-1);
      end
      v=v/N;
      v2=v2/N;
      data(j,i)=v;
      data2(j,i)=v2;
    end
  end
end

f=fopen('_parameters.data','w');
fprintf(f,'h=%d\n',0.274);
fprintf(f,'rho=%d\n',2.419);
fclose(f);
disp('parameters');

f=fopen('_x.data','w');
fprintf(f,'%d\n',length(coord));
for i=1:1:length(coord)
    fprintf(f,'%d\n',coord(i)/1000);
end
fclose(f);
disp('x');

f=fopen('_t.data','w');
fprintf(f,'%d\n',length(time));
for i=1:1:length(time)
    fprintf(f,'%d\n',time(i)-time(1));
end
fclose(f);
disp('t');

f=fopen('_u.data','w');
for i=1:1:length(coord)
    for j=1:1:length(time)
        fprintf(f,'%d\n',data(i,j));
    end
    disp(i);
    disp(length(coord));
end
fclose(f);
disp('u');

f=fopen('_u2.data','w');
for i=1:1:length(coord)
    for j=1:1:length(time)
        fprintf(f,'%d\n',data2(i,j));
    end
    disp(i);
    disp(length(coord));
end
fclose(f);
disp('u2');

