% Model_diffTau_test
% Published Wang and Pena 2013 Journal of Neuroscience

% Wang Y, Peña JL (2013) Direction selectivity mediated by adaptation in the owl’s inferior colliculus. J Neurosci 33:19167–19175. 
%
% Yunyan Wang 2013
% Wangyyj@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%% Adaptation Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



load ('AdaptDirSel_n44_20ms') 
load ('AdaptDirSel_n44_20ms_AzSoloOnly') 

adapt_n=1:length(AdaptDirSel);
adapt_n([12 22])=[]; %no DS info
for n=adapt_n
    A(n,:)=scaledata(AdaptDirSel{1,n}.mean(1:15,5)',min(AdaptDirSel{1,n}.mean(1:15,5)'),1);
    DS(n)=DirSel{1,n}.Stream.DirSel;
end
Adapt_start=1; 
AdaptCurveAll = A(Adapt_start:end,:);
isi=AdaptDirSel{1,11}.isi(1:15);

a_tau=1:1:1000;
t=1:550; %length of exponential
adaptRange=(2:14);
adaptIndex=floor(isi(adaptRange));
%Calculate exponentials
A_Curves=zeros(length(a_tau),length(t));
for n=1:length(a_tau)
    A_Curves(n,:)=1-exp(t/a_tau(n)*-1);
end

%fit all cells, find best tau
adapt_diff=zeros(length(AdaptCurveAll(:,1)),length(a_tau));
AdaptTau=zeros(1,length(AdaptCurveAll(:,1)));
for n=1:length(AdaptCurveAll(:,1))
    for r=1:length(a_tau)
        adapt_diff(n,r)= mean((A_Curves(r,adaptIndex) - (AdaptCurveAll(n,adaptRange))).^2);  
    end
    [temp smallestDiff]=min(adapt_diff(n,:));
    AdaptTau(n)=a_tau(smallestDiff);
end
degSupp=sum(1-AdaptCurveAll(:,adaptRange)');

figure('Color',[1 1 1]);
scatter(AdaptTau,degSupp)
xlabel('Best tau')
ylabel('Adaptation strength')

figure('Color',[1 1 1]);
hist(AdaptTau,20)
xlim([a_tau(1) a_tau(end)])
xlabel('tau (ms)')
ylabel('Incidence')
disp(sprintf('Median tau for Adapt fit = %d ms',median(AdaptTau)))

figure('Color',[1 1 1]);
hold on
plot([isi(adaptRange(1)) isi(adaptRange(end)) ],ones(1,2),'--k')
plot(isi(adaptRange),AdaptCurveAll(6,adaptRange),'r')
ylim([0 1.2])
xlabel('Inter-Click Interval (ms)')
ylabel('Normalized response')
legend('Response to first click','response to second click')
legend('boxoff')

% %Alt = normalize adapt curves
% for n=1:adapt_n
%     AdaptCurveAlln(n,:)=scaledata(AdaptCurveAll(n,:),0,1);
% end
% %fit all cells, find best tau
% for n=1:adapt_n
%     for r=1:length(a_tau)
%         adapt_diffn(n,r)= mean((A_Curves(r,adaptIndex) - (AdaptCurveAlln(n,adaptRange))).^2);  
%     end
%     [temp smallestDiffn]=min(adapt_diffn(n,:));
%     AdaptTauN(n)=a_tau(smallestDiffn);
% end
% degSuppN=sum(1-AdaptCurveAlln(:,adaptRange)');
% scatter(AdaptTauN,degSuppN)
% 
% figure('Color',[1 1 1]);
% hist(AdaptTauN,20)
% xlabel('tau (ms)')
% ylabel('frequency')
% median(AdaptTau)

%%%%%%%%%%%%%%%%%%%%%%%%% Motion Data %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot Mean Adapt Curve
figure('Color',[1 1 1]);
plot(mean(AdaptCurveAll))
title('Mean Adapt Data')

range=(2:13);
n_neurons=length(DS);
tau=1:4:400;
points=21; %Num of Speakers
L=points*20; %length of model (20pts/speaker)

n_neurons=1:n_neurons;
n_neurons([12 22])=[];
for n=n_neurons 
    %normalize
    temp=scaledata([DirSel{1,n}.Stream.LRpsth(2,:); DirSel{1,n}.Stream.RLpsth(2,:)],0,1);
    psth_t=DirSel{1,n}.Stream.LRpsth(1,:)-Lat(n);
    %interp and shift
    for m=1:2
        temp2(m,:)=interp1(psth_t,temp(m,:),1:floor(max(psth_t)));
    end
    PSTH{1,n}=temp2;
    clear ('temp','temp2');
end
Curves=zeros(length(tau),length(t));
for n=1:length(tau)
    Curves(n,:)=exp(t/tau(n)*-1);
end
RF_n=RF;
RF_nRL=fliplr(RF_n); %Flip normalized RF for RL summation calculation
RF_t=[1:21/420:22]; %21 speakers in array, interpolate points so dt=1ms
RF_t=RF_t(1:420);
%define where center is (speaker 8-14, center is at 11)
DSrange(1)=find(RF_t>=8,1,'first');
DSrange(2)=find(RF_t>=14,1,'first');
% Find adapt curve range 
adaptIndex=floor(isi(range));

for r=1:length(tau)
    for n=n_neurons
        RF_longLR(n,:)=interp1(1:21,RF_n(n,:),RF_t); %interpolate points for RF
        RF_longRL(n,:)=interp1(1:21,RF_nRL(n,:),RF_t);
        suppLR{n,r}=zeros(points,L);
        suppLR{n,r}=zeros(points,L);

        %first point LR
        if RF_n(n,1)>0
            suppLR{n,r}(1,1:L)= Curves(r,1:L)*RF_n(n,1);
        else
            suppLR{n,r}(1,1:L)= nan(1,L);
        end
        for m=2:points %dt
            if RF_n(n,m)>0
                temp=suppLR{n,r}(m-1,(m-1)*20); %find last value from last stim
                if ~isnan(temp)
                    suppLR{n,r}(m,(m-1)*20+1:L)= Curves(r,1:L-20*(m-1))*RF_n(n,m)+temp;
                else
                    suppLR{n,r}(m,(m-1)*20+1:L)=Curves(r,1:L-20*(m-1))*RF_n(n,m);
                end
            else
                suppLR{n,r}(m,(m-1)*20+1:L)= suppLR{n,r}(m-1,(m-1)*20:L-1); %if no activity, keep decaying
            end
        end
        %first point RL
        if RF_nRL(n,1)>0
            suppRL{n,r}(1,1:L)= Curves(r,1:L)*RF_nRL(n,1);
        else
            suppRL{n,r}(1,1:L)= nan(1,L);
        end
        for k=2:points
            if RF_nRL(n,k)>0
                temp2=suppRL{n,r}(k-1,(k-1)*20); %find last value from last stim
                if ~isnan(temp2)
                    suppRL{n,r}(k,(k-1)*20+1:L)= Curves(r,1:L-20*(k-1))*RF_nRL(n,k)+temp2;
                else
                    suppRL{n,r}(k,(k-1)*20+1:L)=Curves(r,1:L-20*(k-1))*RF_nRL(n,k);
                end
            else
                suppRL{n,r}(k,(k-1)*20+1:L)= suppRL{n,r}(k-1,(k-1)*20:L-1); %if no activity, keep decaying
            end
        end
        %summate
        for p=1:points
            EsuppLR{1,n}(r,1+(20*(p-1)):p*20)=suppLR{n,r}(p,1+(20*(p-1)):p*20);
            EsuppRL{1,n}(r,1+(20*(p-1)):p*20)=suppRL{n,r}(p,1+(20*(p-1)):p*20);
        end
    end
end
for n=n_neurons
    EsuppS{1,n}=1-scaledata([EsuppLR{1,n} EsuppRL{1,n}],0,1);
    for r=1:length(tau)
        %multiply with static RF
        RF_modelLR{1,n}(r,:)= EsuppS{1,n}(r,1:L).* RF_longLR(n,:);
        RF_modelRL{1,n}(r,:)= EsuppS{1,n}(r,L+1:L*2).* RF_longRL(n,:);
        tempa=RF_modelLR{1,n}(r,DSrange(1):DSrange(2));
        tempb=RF_modelRL{1,n}(r,DSrange(1):DSrange(2));
        DSmodel(r,n)=calcIndex_s(sum(tempa(~isnan(tempa))),sum(tempb(~isnan(tempb))));
    
        %Calculate difference between model and actual psth 
        psth_l=length(PSTH{1,n}(1,:));
        temp1=(RF_modelLR{1,n}(r,1:psth_l) - PSTH{1,n}(1,1:psth_l)).^2;
        CurveDiff{r,n}(1,:)=mean(temp1(~isnan(temp1)));
        temp2=(RF_modelLR{1,n}(r,1:psth_l) - PSTH{1,n}(2,1:psth_l)).^2;
        CurveDiff{r,n}(2,:)=mean(temp2(~isnan(temp2)));
        CurveDiffTot(r,n)=sum([CurveDiff{r,n}(1,:) CurveDiff{r,n}(2,:)]);
        
        clear('temp1', 'temp2','tempa','tempb','psth_l')
    end   
    DSDiff(:,n)=abs(DSmodel(:,n)-DS(n));
end
[tau_i bestCurveTau]=min(CurveDiffTot);
[tau_ds bestDSITau]=min(DSDiff);
for n=n_neurons
    bestDSTau(n)= DSmodel(bestCurveTau(n),n);
    TauFitMin(n,1)=find((CurveDiffTot(:,n)<=tau_i(n)+(max(CurveDiffTot(:,n))-tau_i(n))*0.05),1,'first');  %within 5% of best fit
    DSmodel5(n)=DSmodel(TauFitMin(n),n);
end

%plot suppression
figure('Color',[1 1 1]);
plot(1-EsuppS{1,4}(:,1:420)')
xlabel('Time (ms)')
ylabel('Suppression')

plotrange=[1 2 3 4 5 length(tau)];
for n=plotrange
    figure('Color',[1 1 1]);
    hold on
    scatter(DS,DSmodel(n,:))
    plot([-1 1],[-1 1],'--k')
    xlim([-1 1])
    ylim([-1 1]) 
    legend(sprintf('tau=%d ms',tau(n)))
    xlabel('DSI')
    ylabel('Model DSI')
end

%Data Fitting
fit_opts = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'on');
xr=[-1 1];

%Fit tau=1000
[fit_1000, goodness1000] = fit(DS(n_neurons)',DSmodel(length(tau),n_neurons)', 'poly1', fit_opts);
[b1000,stats1000] = robustfit(DS(n_neurons)',DSmodel(length(tau),n_neurons)');
regr_form1000(1:2)=b1000';
r2_1000=goodness1000.rsquare;

%Fit with optimum tau
[fit_Op, goodnessOp] = fit(DS(n_neurons)', bestDSTau(n_neurons)', 'poly1', fit_opts);
[bOp,statsOp] = robustfit(DS(n_neurons)',bestDSTau(n_neurons));
regr_formOp(1:2)=bOp';
r2_Op=goodnessOp.rsquare;

%Fit with tau within 5% of optimum tau
[fit_Op5, goodnessOp5] = fit(DS(n_neurons)', DSmodel5(n_neurons)', 'poly1', fit_opts);
[bOp5,statsOp5] = robustfit(DS(n_neurons)',DSmodel5(n_neurons)');
regr_formOp5(1:2)=bOp5';
r2_Op5=goodnessOp5.rsquare;

%plot DS prediction tau=1000
figure('Color',[1 1 1]);
hold on
scatter(DS(n_neurons), DSmodel(length(tau),n_neurons))
plot(xr,regr_form1000(1)+regr_form1000(2)*xr,'--r')
plot([-1 1],[-1 1],'--k')
xlim([-1 1])
ylim([-1 1])
xlabel('DSI')
ylabel('Model DSI')
title('Model DSI tau=1000ms')

%plot best predicted DS
figure('Color',[1 1 1]);
hold on
scatter(DS, bestDSTau,'r')
plot(xr,regr_formOp(1)+regr_formOp(2)*xr,'--r')
plot([-1 1],[-1 1],'--k')
xlim([-1 1])
ylim([-1 1])
xlabel('DSI')
ylabel('Model DSI')
title('Model DSI individual optimum tau')

%plot 5% within best predicted DS
figure('Color',[1 1 1]);
hold on
scatter(DS(n_neurons),DSmodel5(n_neurons),'r')
plot(xr,regr_formOp5(1)+regr_formOp5(2)*xr,'--r')
plot([-1 1],[-1 1],'--k')
xlim([-1 1])
ylim([-1 1])
xlabel('DSI')
ylabel('Model DSI')
title('Model DSI individual 5% of optimum tau')

%plot best taus for motion data
figure('Color',[1 1 1]);
hist(tau(bestCurveTau))
xlabel('tau (ms)')
ylabel('Incidence')
title('Tau of best fit')

%plot 5% within best taus for motion data
%plot 5% within best taus for motion data
figure('Color',[1 1 1]);
hist(tau(TauFitMin(n_neurons)'))
xlim([tau(1) tau(end)])
xlabel('tau (ms)')
ylabel('Incidence')
title('Tau within 10% of best fit')
title('Tau of best fit')
disp(sprintf('Median tau for Adapt fit = %d ms',median(tau(TauFitMin(n_neurons)))))

%plot best predicted DS
figure('Color',[1 1 1]);
hold on
plot([-1 1],[-1 1],'--k')
for n=n_neurons
    if tau(bestCurveTau(n))>tau(end)/2
        scatter(DS(n), bestDSTau(n),'r')
    else
        scatter(DS(n), bestDSTau(n),'b')
    end
end
xlabel('DSI')
ylabel('Model DSI')

%Comparison: adapt curve tau and DSI prediction tau
AdaptSeqTau(:,1)=(a_tau(AdaptTau(n_neurons)));
AdaptSeqTau(:,2)=(tau(TauFitMin(n_neurons)));
figure('Color',[1 1 1]);
scatter(AdaptSeqTau(:,1),AdaptSeqTau(:,2))
xlabel('Adaptation tau')
ylabel('DS tau')

figure('Color',[1 1 1]);
scatter(abs(DS(n_neurons)),tau(TauFitMin(n_neurons)))
xlabel('DS')
ylabel('Tau within 5 percent optimum')
