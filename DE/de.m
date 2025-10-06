%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;
for degree=50:5:50
    tic
    rng('default')
wp=0.0625;
ws=0.125;
% degree=5;
%% Desired
    [hd,wd]=HDesired(wp,ws,1024);
    %% Equiripple
    eqfilt=designfilt('lowpassfir', 'FilterOrder', degree, 'PassbandFrequency', wp, 'StopbandFrequency', ws);
    [h1,w1]=freqz(eqfilt.Coefficients,1,1024);
    h1=abs(h1);
    w1=w1/pi;
%% Problem Definition

% CostFunction=@(x) Sphere(x);    % Cost Function

nVar=degree;            % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=-1;          % Lower Bound of Decision Variables
VarMax= 1;          % Upper Bound of Decision Variables

%% DE Parameters

MaxIt=1000;      % Maximum Number of Iterations

nPop=100;        % Population Size

beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor

pCR=0.2;        % Crossover Probability

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop

    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=mycostfunc(pop(i).Position,hd,wp,ws,degree);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
		y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=mycostfunc(NewSol.Position,hd,wp,ws,degree);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    
    [h,w]=freqz(BestSol.Position,1,1024);
                h=abs(h);
                w=w/pi;
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
     plotting(w,h,wd,hd,w1,h1,degree)
            
end

%% Show Results

% figure;
%plot(BestCost);
% semilogy(BestCost, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
 elapsedtime=toc;

    filename=['20de0625'];
    save(filename);
    
    %% End of Program
    
    clear all;
    close all;
    clc;
% [h,w]=freqz([BestSol.Position,fliplr(BestSol.Position(1:end-1))],1,2048,2);
% [h,w]=freqz(BestSol.Position,1,1024);
% figure;plot(w/pi,abs(h),'LineWidth',2)
end