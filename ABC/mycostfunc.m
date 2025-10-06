function err= mycostfunc(ind,hd,wp,ws,degree)
[h,w]=freqz(ind,1,1024);
w=w/pi;
h=abs(h);
hd=hd';


err=sum((h-hd).^2);
end
