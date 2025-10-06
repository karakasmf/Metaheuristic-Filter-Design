%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA114
% Project Title: Implementation of Artificial Bee Colony in MATLAB
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

% CostFunction=@(x) Sphere(x);        % Cost Function

nVar=degree;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=-1;         % Decision Variables Lower Bound
VarMax= 1;         % Decision Variables Upper Bound

%% ABC Settings

MaxIt=10000;              % Maximum Number of Iterations

nPop=100;               % Population Size (Colony Size)

nOnlooker=nPop;         % Number of Onlooker Bees

L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)

a=1;                    % Acceleration Coefficient Upper Bound

%% Initialization

% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];

% Initialize Population Array
pop=repmat(empty_bee,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Population
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=mycostfunc(pop(i).Position,hd,wp,ws,degree);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% ABC Main Loop

for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=mycostfunc(newbee.Position,hd,wp,ws,degree);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=mycostfunc(newbee.Position,hd,wp,ws,degree);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost=mycostfunc(pop(i).Position,hd,wp,ws,degree);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    BestCost(it)=BestSol.Cost;
    
    [h,w]=freqz(BestSol.Position,1,1024);
                h=abs(h);
                w=w/pi;
     plotting(w,h,wd,hd,w1,h1,degree)
end
    
%% Results
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
 elapsedtime=toc;

    filename=['120abc0625'];
    save(filename);
    
    %% End of Program
    
    clear all;
    close all;
    clc;
% [h,w]=freqz([BestSol.Position,fliplr(BestSol.Position(1:end-1))],1,2048,2);
% [h,w]=freqz(BestSol.Position,1,1024);
% figure;plot(w/pi,abs(h),'LineWidth',2)
end