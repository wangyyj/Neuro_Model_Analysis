% calc_binaural_interaction
% Takes monaural stimuli and ILD curves to calculate binaural interaction
% REFERENCE
close all;
clear;
AnimalID = 'M117A';

% channel_no = input('Channel number?       ');
% % Enter monaural files (L/R)
% MonoL_file = input('Enter file # for monaural LEFT     ','s');
% MonoR_file = input('Enter file # for monaural RIGHT    ','s');
% MonoL_rate = calc_rate([AnimalID MonoL_file], channel_no);
% MonoR_rate = calc_rate([AnimalID MonoR_file], channel_no);
% 
% % Enter ILD files
%     %Enter another (y)? Done (n)
% ILD_file{1} = input('Enter ILD filename(#)       ','s');
% ILD_rate{1} = calc_rate([AnimalID ILD_file{1}], channel_no);
% nILD = 1;
% moreILD = input('Enter another (y)? Done (n)       ','s');
% while moreILD == 'y'
%     nILD = nILD+1;
%     ILD_file{nILD} = input('Enter ILD filename(#)       ','s');
%     ILD_rate{nILD} = calc_rate([AnimalID ILD_file{nILD}], channel_no);
%     moreILD = input('Enter another (y)? Done (n)       ','s');
% end

% Quick input, enter 5 files, plot
nILD=3;
channel_no = 4;
files_s = input('enter file names separate by space   ','s');
files_s = files_s(~isspace(files_s));
place = 1;
for n=1:5
    files{n}= files_s(place:place+3);
    place=place+4;
end
MonoL_file = num2str(files{1});
MonoR_file = num2str(files{2});
MonoL_rate = calc_rate([AnimalID MonoL_file], channel_no);
MonoR_rate = calc_rate([AnimalID MonoR_file], channel_no);
for n=1:3
    ILD_file{n} = num2str(files{2+n});
    ILD_rate{n} = calc_rate([AnimalID ILD_file{n}], channel_no);
end

% Find predicted firing rate by linearly combining L and R
for m = 1:nILD
    ILD_L = ILD_rate{1,m}.attn';
    ILD_R = fliplr(ILD_L);
    Lrate = MonoL_rate.tuning;
    RrateR = MonoR_rate.tuning;
    output.ABL(m) = ILD_L(ceil(length(ILD_L)/2)); %ABL
    r =1;
    for n = 1:length(ILD_L)
        [memL locL] = ismember(ILD_L(n), MonoL_rate.attn);
        [memR locR] = ismember(ILD_R(n), MonoR_rate.attn);
        if memL == 1 && memR ==1
            output.rate{1,m}(r,1) = ILD_rate{1,1}.ILD(n); %ILD
            output.rate{1,m}(r,2) = ILD_L(n); %Lattn
            output.rate{1,m}(r,3) = Lrate(locL);%Lresp
            output.rate{1,m}(r,4) = ILD_R(n);%Rattn
            output.rate{1,m}(r,5) = RrateR(locR);%Rresp
            output.rate{1,m}(r,6) = Lrate(locL) + RrateR(locR); %sum of L+R
            output.rate{1,m}(r,7) = ILD_rate{1,m}.tuning(n);%binaural resp
            r = r+1;
        end
    end
    ymax(m) = max([output.rate{1,m}(:,6)' output.rate{1,m}(:,7)' ILD_rate{1,m}.tuning]);
    output.compare(m,:) = output.rate{1,m}(:,6)'-output.rate{1,m}(:,7)';
end

%plot
init_fig_arial
subplot(nILD+1,1,1)
plot(MonoL_rate.attn,MonoL_rate.tuning,'g');
hold on
plot(MonoR_rate.attn,MonoR_rate.tuning,'r');
title('Monaural Response','fontweight','bold')
xlabel('Attn (dB)')
ylabel('Spikes/s')
legend('L','R','location','northeast')
legend('boxoff')
box off
ymax = max(ymax); %find max, set same ylim
for n = 1:nILD
    subplot(nILD+1,1,n+1)
    scatter(output.rate{1,n}(:,1),output.rate{1,n}(:,6),'r','filled')
    hold on
    plot(ILD_rate{1,n}.ILD,ILD_rate{1,n}.tuning,'k')
    box off
    ylim([0 ymax])
    xlabel('ILD (dB)')
    ylabel('Spikes/s')
    title(['ABL = ' num2str(output.ABL(n)) '   ' AnimalID ILD_file{n}],'fontweight','bold');
    % legend('Sum Monaural','Binaural','location','northwest')
    % legend('boxoff')
end

% Plot difference between predicted and
ILDs = output.rate{1,1}(:,1)';
init_fig_arial
plot(ILDs,output.compare')
hold on
plot([ILDs(1) ILDs(end)],[0 0],'k--')
l=legend(num2str(output.ABL(1)),num2str(output.ABL(2)),num2str(output.ABL(3)),'location','northwestoutside');
v = get(l,'title');
set(v,'string','ABL');
xlabel('ILD')
ylabel('(L+R) - (Binaural)')
box off





% Facilitating/suppressive (sign) 

% Categorize 

%Save analysis file
%savepath =
%'C:\Dropbox\WangLab\Data\M117A\AnalyzedData\BinauralInteraction';
%savefilename = [AnimalID 'U' Unit '_BI' '.mat'];

    
