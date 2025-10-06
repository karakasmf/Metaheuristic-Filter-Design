function err= mycostfunc(ind,hd,wp,ws,degree)



    [h,w]=freqz(ind,1,1024);
    w=w/pi;
    h=abs(h);
    [pr1 locsp1]=max(findpeaks(h(find(w<wp))));
    if locsp1>=1
        pr1=abs(1-max(pr1));
    else
        pr1=0.8;
    end
    [pr2 locsp2]=min(findpeaks(-h(find(w<wp))));
    if locsp2>=1
        pr2=abs(1-abs(max(pr2)));
    else
        pr2=0.8;
    end
    
    [sr locss1]=max(findpeaks(h(find(w>ws))));
    if locss1>=1
        sr=abs(max(sr));
    else
        sr=0.8;
    end
    
    Ep=sum((1-h(find(w<wp))).^2);
    Es=sum(h(find(w>ws)).^2);
    Et=(median(h(ws>w>wp))-0.707).^2;
    cc=corrcoef(hd,h);
    cm2=(1-cc(1,2));
    C=sum([pr1 pr2 sr Ep Es Et cm2]);
    err= C;
end
