clear all;

d = load("px0,5py0.txt");
d = load("px49py0.txt");
%plot(d,'.-');
%hold('on');

fh=fopen("rewrite.txt","w");
for x=0:0.5:50
    intx=floor(x);
    floatx=x-floor(x);
    while(abs(floor(floatx)-floatx)>0.0001)
        floatx=floatx*10;
    end
    if(floatx<0.00001)
        file = "px"+string(intx)+"py0.txt";
    else
        file = "px"+string(intx)+","+string(floatx)+"py0.txt";
    end
    d=load(file);
    
    if(max(abs(d(1:2500)))>2.3)
        disp(x);
        fwrite(fh,string(x));
        fwrite(fh," ");
    end
end
fclose(fh);
