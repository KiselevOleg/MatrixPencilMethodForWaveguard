clear all;

scan_file_name="scan.mat";

get_point_for_t_from=3000;
get_point_for_t_to=8000;
get_only_Nth_point_for_t_where_N=1;

skip_first_point_for_x=0;

known_parameters=[
  "material=Al",
  "h=0.329cm",
  "rho=null"
];



N=get_only_Nth_point_for_t_where_N;
load(scan_file_name);

data_=data(skip_first_point_for_x+1:length(coord_),get_point_for_t_from:get_point_for_t_to);
coord_=coord_(skip_first_point_for_x+1:length(coord_));
time_=time_(get_point_for_t_from:get_point_for_t_to);

coord=coord_;
time=[];
data=[];

disp("start")
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

f=fopen("parameters.data","w");
for i=1:1:size(known_parameters,1)
  fprintf(f,strcat(known_parameters(i,:),"\n"));
end
%fprintf(f,"h=%d\n",0.274);
%fprintf(f,"rho=%d\n",2.419);
fclose(f);
disp("parameters");

f=fopen("x.data",'w');
fprintf(f,"%d\n",length(coord));
for i=1:1:length(coord)
    fprintf(f,"%d\n",coord(i));
end
fclose(f);
disp("x");

f=fopen("t.data","w");
fprintf(f,"%d\n",length(time));
for i=1:1:length(time)
    fprintf(f,"%d\n",time(i)-time(1));
end
fclose(f);
disp("t");

f=fopen("u.data","w");
for i=1:1:length(coord)
    for j=1:1:length(time)
        fprintf(f,"%d\n",data(i,j));
    end
    disp(i);
    disp(length(coord));
end
fclose(f);
disp("u");

