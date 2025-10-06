%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA102
% Project Title: Implementation of Particle Swarm Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Particle Swarm Optimization in MATLAB (URL: https://yarpiz.com/50/ypea102-particle-swarm-optimization), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;
for degree=50:50
    tic
    rng('default')
wp=0.25;
ws=0.25;
% degree=5;
%% Desired
    [hd,wd]=HDesired(wp,ws,1024);
    %% Equiripple
%     eqfilt=designfilt('lowpassfir', 'FilterOrder', degree, 'PassbandFrequency', wp, 'StopbandFrequency', ws);
%     [h1,w1]=freqz(eqfilt.Coefficients,1,1024);
%     h1=abs(h1);
%     w1=w1/pi;
    
%% Problem Definition

% CostFunction=@(x) Sphere(x);        % Cost Function

nVar=degree;             % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

VarMin = -1;         % Lower Bound of Variables
VarMax = 1;         % Upper Bound of Variables


%% PSO Parameters

MaxIt = 1000;      % Maximum Number of Iterations

nPop = 100;        % Population Size (Swarm Size)

% PSO Parameters
w = 1;            % Inertia Weight
wdamp = 0.99;     % Inertia Weight Damping Ratio
c1 = 1.5;         % Personal Learning Coefficient
c2 = 2.0;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO, 
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1 = 2.05;
% phi2 = 2.05;
% phi = phi1+phi2;
% chi = 2/(phi-2+sqrt(phi^2-4*phi));
% w = chi;          % Inertia Weight
% wdamp = 1;        % Inertia Weight Damping Ratio
% c1 = chi*phi1;    % Personal Learning Coefficient
% c2 = chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax = 0.1*(VarMax-VarMin);
VelMin = -VelMax;

%% Initialization

empty_particle.Position = [];
empty_particle.Cost = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);

GlobalBest.Cost = inf;

for i = 1:nPop
    
    % Initialize Position
    particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    particle(i).Cost = mycostfunc(particle(i).Position,hd,wp,ws,degree);
    
    % Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest = particle(i).Best;
        
    end
    
end

BestCost = zeros(MaxIt, 1);

%% PSO Main Loop

for it = 1:MaxIt
    
    for i = 1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin);
        particle(i).Velocity = min(particle(i).Velocity, VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside = (particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside) = -particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        
        % Evaluation
        particle(i).Cost = mycostfunc(particle(i).Position,hd,wp,ws,degree);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest = particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it) = GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w = w*wdamp;
        [he,we]=freqz(GlobalBest.Position,1,1024);
                he=abs(he);
                we=we/pi;
     plotting(we,he,wd,hd,degree)
end

BestSol = GlobalBest;

 elapsedtime=toc;

    filename=['120pso0625'];
    save(filename);
    
    %% End of Program
    
    clear all;
    close all;
    clc;
% [h,w]=freqz([BestSol.Position,fliplr(BestSol.Position(1:end-1))],1,2048,2);
% [h,w]=freqz(BestSol.Position,1,1024);
% figure;plot(w/pi,abs(h),'LineWidth',2)
end
