function x = bwo(ws,wp,N,dimension,nPop,maxiter)

rng(rand)

% Filter Parameters

ub=1;
lb=-1;
pc=0.8;                                                                     % Percent of Crossover
nCross=round(pc*nPop/2)*2;                                                  % Number of Selected Parents
pMutation=0.4;
nMutation=round(pMutation*nPop);                                            % Number of Mutants
pCannibalism=0.5;                                                           % Percent of Canniibalism
nCannibalism=round(pCannibalism*dimension);  




