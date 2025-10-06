
%% Statrt of Program
clear all
close all
clc
for degree=20:5:20
    rng('default')
    
    tic
    
    
    
    wp=0.0624;
    ws=0.0626;
    wc=(wp+ws)/2;
    
    PopSize = 100;
    MaxIteration = 2500;
    
    % degree=30;
    DimNum=degree;
    
    x=designfilt('lowpassfir', 'FilterOrder', degree, 'PassbandFrequency', wp, 'StopbandFrequency', ws);
    [h1,w1]=freqz(x.Coefficients,1,1024);
    h1=abs(h1);
    w1=w1/pi;
    
    % % iseven=rem(degree,2);
    % % if iseven==0
    % %     DimNum = (degree+2)/2;
    % % else
    % %     DimNum = (degree+1)/2;
    % % end
    
    
    
    
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
    
    
    %% Algorithm Parameters
    SelMethod = 1;
    CrossMethod = 1;
    
    
    
    CrossPercent = 65;
    MutatPercent = 10;
    ElitPercent = 100 - CrossPercent - MutatPercent;
    
    CrossNum = round(CrossPercent/100*PopSize);
    
    if mod(CrossNum,2)~=0
        CrossNum = CrossNum - 1;
    end
    
    MutatNum = round(MutatPercent/100*PopSize);
    ElitNum = PopSize - CrossNum - MutatNum;
    
    %% Problem Satement
    lb = -1;
    ub =1;
    
    CostFuncName =@MyFucn_Fcn;
    
    %% Initial Population
    Population = rand(PopSize,DimNum) * (ub - lb) + lb;
    Cost = feval(CostFuncName,Population,Hd,wp,ws,degree);
    [Cost Indx] = sort(Cost);
    Population = Population(Indx,:);
    
    %% Main Loop
    MeanMat = [];
    MinMat = [];
    BestC=10000;
    Iter=1;
    % hWaitbar = waitbar(0.5, [num2str(degree) ' Coeff'],...
    %     'Name', 'Problem Çözülüyor','CreateCancelBtn','delete(gcbf)');
    itercost=[];
    % while BestC>0.5
    for Iter = 1:MaxIteration
        %% Elitism
        ElitPop = Population(1:ElitNum,:);
        
        %% Cross Over
        CrossPop = [];
        ParentIndexes = GASelectParents_Fcn(Cost,CrossNum,SelMethod);
        
        for ii = 1:CrossNum/2
            Par1Indx = ParentIndexes(ii*2-1);
            Par2Indx = ParentIndexes(ii*2);
            
            Par1 = Population(Par1Indx,:);
            Par2 = Population(Par2Indx,:);
            
            [Off1 Off2] = GAMyCrossOver_Fcn(Par1,Par2,CrossMethod);
            CrossPop = [CrossPop ; Off1 ; Off2];
        end
        
        %% Mutation
        MutatPop = rand(MutatNum,DimNum) * (ub - lb) + lb;
        
        %% New Population
        Population = [ElitPop ; CrossPop ; MutatPop];
        Cost = feval(CostFuncName,Population,Hd,wp,ws,degree);
        [Cost Indx] = sort(Cost);
        Population = Population(Indx,:);
        %% Algorithm Progress
        BestC = Cost(1);
        disp('----------------------------------------------')
        disp(['Iter= ',num2str(Iter),' BestCost= ', num2str(BestC)])
        
        %BestP = Pop(1,:)
        MinMat(Iter) = Cost(1);
        MeanMat(Iter) = mean(Cost);
        p=Population(1,:);
        
        % % % %     iseven=rem(degree,2);
        % % % % if iseven==0
        % % % %     [h,w]=freqz([p fliplr(p(1:end-1))],1,1024);
        % % % % %     [p fliplr(p(1:end-1))]
        % % % % else
        % % % %     [h,w]=freqz([p fliplr(p)],1,1024);
        % % % % %     [p fliplr(p)]
        % % % % end
        errors(Iter)=BestC;
        [h,w]=freqz(Population(1,:),1,1024);
        
%         plot(w/pi,abs(h),'LineWidth',2);
        
        %hold on
        %plot(w0,H0,'LineWidth',2);
        %plot(w1,h1,'LineWidth',2);
        %legend('Estimated','Desired','Equiripple')
        %title(['Degree= ' num2str(degree)])
%         drawnow;
        %hold off
        %     if ~ishandle(hWaitbar)
        %         % Stop the if cancel button was pressed
        %         disp('Stopped by user');
        %         break;
        %     end
        
        %plot(MinMat,'--r','linewidth',2);
        %hold on
        %     plot(MeanMat,'--k','linewidth',2);
        %     hold off
        %pause(.5)
        %     semilogy(Iter,MinMat(Iter),'r.')
        %hold on
        %semilogy(Iter,MeanMat(Iter),'b.')
        Iter=Iter+1;
        itercost(end+1)=BestC;
    end
    % ylim([0 5])
    %% Results
    BestSolution = Population(1,:);
    
    % % % iseven=rem(degree,2);
    % % % if iseven==0
    % % %    BestSolution1=[BestSolution fliplr(BestSolution(1:end-1))]
    % % % else
    % % %     BestSolution1=[BestSolution fliplr(BestSolution)]
    % % % end
    BestCost = Cost(1,:)
    [h,w]=freqz(BestSolution,1,1024);
    % figure;plot(w/pi,abs(h),'LineWidth',2)
    % hold on
    % plot(w0,H0,'LineWidth',2);
    % hold off
    % figure,plot(itercost)
    elapsedtime=toc;
    filename=['20ga0625'];
    save(filename);
    
    %% End of Program
    
    clear all;
    close all;
    clc;
end


function [Off1 Off2] = MyCrossOver_Fcn(Par1,Par2,CrossMethod)

switch CrossMethod
    case 1
        Beta1 = rand;
        Beta2 = rand;
        
        Off1 = Beta1*Par1 + (1-Beta1)*Par2;
        Off2 = Beta2*Par1 + (1-Beta2)*Par2;
    case 2
        
    case 3
        
end

function ParIndexes = SelectParents_Fcn(Cost,SelectionNum,SelMethod)
PopSize = size(Cost,1);

switch SelMethod
    case 1
        R = randperm(PopSize); ParIndexes = R(1:SelectionNum);
end

end