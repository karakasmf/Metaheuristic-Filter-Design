clear all
close all
clc

for dimension=20:5:20
    tic
rng('default')

wp=0.0624;
ws=0.0626;
% dimension =30;
degree=dimension;
x=designfilt('lowpassfir', 'FilterOrder', degree, 'PassbandFrequency', wp, 'StopbandFrequency', ws);
[h1,w1]=freqz(x.Coefficients,1,1024);
h1=abs(h1);
w1=w1/pi;

%% Desired Filter


N=1024;
w0=1/N:1/N:1;
%%%%%%%%%%%%%%%%%
H0=zeros(1,length(w0));
[~,n]=find(w0<=wp);
na=max(n)+1;
[~,n2]=find(w0>ws);
nu=min(n2)-1;
uz=(nu-1)-(na+1);
hx=linspace(1.00001,0.00001,uz);
H0(1,1:na)=1;
for i=1:uz
    H0(1,na+i)=hx(1,i);
end
H0(1,nu:end)=0;

% H0(w0<wc)=1;
Hd=H0;

%% SSA & Optimization Parameters

s_count=100;
iter =2500;
predator=0.1;
gliding_constant=1.9;
dg=0.5+(1.11-0.5)*rand();
acorntree=ceil(s_count*6/100);
hickorytree=1;
normaltree=s_count-acorntree-hickorytree;
lowerbound=-1;
upperbound=1;
current_iter=1;
sc=0;
smin=0;


wintercount=0;
summercount=0;

%% START

position=zeros(s_count,dimension);
best=zeros(iter,dimension);

for i=1:s_count
    for j=1:dimension
        position(i,j)= lowerbound+(upperbound-lowerbound)*rand();
    end
end



while current_iter<=iter+1
    %%
    error=cost(position,Hd,wp,ws,degree);
    [costss,Index]=sort(error,'ascend');
    position=position(Index,:);
    errors(current_iter)=error(Index(1));
    
    hickory=position(1,:);
    acorn=position(2:1+acorntree,:);
    normal=position(acorntree+2:end,:);
    normalrandom=randperm(normaltree,acorntree);
    best(current_iter,:)=hickory;
    
    %% MOVEMENT
    
    %ACORN TO HICKORY
    for i=1:acorntree
        if rand()>=predator
            acorn(i,:)=acorn(i,:)+dg*gliding_constant*(hickory(1,:)-acorn(i,:));
        else
            acorn(i,:)=lowerbound+(upperbound-lowerbound)*rand(1,dimension);
        end
    end
    %NORMAL TO ACORN
    normaltoacornmove=randi(normaltree);
    if rand()>=predator
        normal(normaltoacornmove,:)=normal(normaltoacornmove,:)+dg*gliding_constant*(acorn(randi(acorntree),:)-normal(normaltoacornmove,:));
    else
        normal(normaltoacornmove,:)=lowerbound+(upperbound-lowerbound)*rand(1,dimension);
    end
    %NORMAL TO HICKORY
    normaltohickorymove=randi(normaltree);
    if rand()>=predator
        normal(normaltohickorymove,:)=normal(normaltohickorymove,:)+dg*gliding_constant*(hickory(1,:)-normal(normaltohickorymove,:));
    else
        normal(normaltohickorymove,:)=lowerbound+(upperbound-lowerbound)*rand(1,dimension);
    end
    position=[hickory;acorn;normal];
    
    %%
    current_iter=current_iter+1;
    
    sc=sqrt(sum((sum(abs(hickory-acorn)).^2)));
    smin=(10*exp(-6))/(365).^(current_iter/(iter/2.5));
    
    if(sc<smin)
        season=summer;
        summercount=summercount+1;
        for i=1:normaltree
            sig=(gamma(2.5)*sin(pi*1.25))/(gamma(1.25)*2.5*2^((0.5)/2));
            levy=0.01*((rand()*sig)/(rand^(1/1.5)));
            normal(i,:)=lowerbound+levy*(upperbound-lowerbound);
        end
        
    else
        season=winter;
        wintercount=wintercount+1;
    end
    
    disp(['Iteration= ' num2str(current_iter) ' Best Cost SSA = ' num2str(costss(1))])
BestCostssa(current_iter)=costss(1);   
[h,w]=freqz(position(1,:),1,1024);
%     plot(w/pi,abs(h),'LineWidth',2);
%     hold on
%     plot(w0,Hd,'LineWidth',2);
%     plot(w1,h1,'LineWidth',2);
%     legend('Estimated','Desired','Equiripple')
%     title(['Degree= ' num2str(degree)])
%     drawnow;
%     hold off
    
end

 elapsedtime=toc;

    filename=['20ssa0625'];
    save(filename);
    
    %% End of Program
    
    clear all;
    close all;
    clc;
% [h,w]=freqz([BestSol.Position,fliplr(BestSol.Position(1:end-1))],1,2048,2);
% [h,w]=freqz(BestSol.Position,1,1024);
% figure;plot(w/pi,abs(h),'LineWidth',2)
end