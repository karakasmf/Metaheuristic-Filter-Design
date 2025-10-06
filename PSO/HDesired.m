function [hd,wd] = HDesired(wp,ws,N)

wd=1/N:1/N:1;
H0=zeros(1,length(wd));
[~,n]=find(wd<=wp);
na=max(n)+1;
[~,n2]=find(wd>ws);
nu=min(n2)-1;
uz=(nu-1)-(na+1);
hx=linspace(1.00001,0.00001,uz);
H0(1,1:na)=1;
for i=1:uz
    H0(1,na+i)=hx(1,i);
end
H0(1,nu:end)=0;
hd=H0;

end

