function err= cost(Pop,Hd,wp,ws,degree)

Cost = zeros(size(Pop,1),1);


for ii = 1:size(Pop,1)
    p = Pop(ii,:);
[h,w]=freqz(p,1,1024);
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
% % % %  pr1=abs(1-max(h(find(w<wp))));
% % % %   pr2=abs(1-min(h(find(w<wp))));
% % % %    sr=max(h(find(w>ws)));


    Ep=sum((1-h(find(w<wp))).^2);
    Es=sum(h(find(w>ws)).^2);
    Et=(median(h(ws>w>wp))-0.707).^2;
    cc=corrcoef(Hd,h);
    cm=mse(abs(Hd-h));
    cm2=(1-cc(1,2));
%     C=pr1+pr2+sr+Ep+Es+Et+cm2;
%     C=max([pr1 pr2 sr Ep Es Et cm2]);
C=sum([cm2 cm]);
    err(ii,1) = C;
end
end
