% Calc_rate
% Calculates rate curve from spike file
% Modified from get_tuning by Chia-Jung Chang
% Created 6/2015 YYWang
function [rate] = calc_rate(filename, channel_no)
% rate.file     full file name
% rate.tuning   x-axis data
% rate.stimType stimulus type
% rate.attn     attenuation
% rate.dur
% rate.spkr
% rate.ILD
% rate.pre_stim
% rate.postStim
% rate.spont
% rate.rep
% rate.var      variable parameter  
% rate.y        y-axis data
% rate.se       standard deviation of the mean
% rate.sd       standard deviation


rate.file = filename;
x = feval(filename);  %ex. ;M117A0856'
rate.stimType = x.analysis_type;
rate.pre_stim = x.pre_stimulus_record_time;
rate.post_stim = x.post_stimulus_record_time;
rate.rep = max(x.data(:,2));

% Get test parameters
attn = x.stimulus_ch1(:, strcmp('attn(dB)', x.stimulus_tags_ch1));
dur = x.stimulus_ch1(:, strcmp('len(ms)', x.stimulus_tags_ch1));
spkr = x.stimulus_ch1(:, strcmp('spkr', x.stimulus_tags_ch1));
pre_stim = x.pre_stimulus_record_time;
post_stim = x.post_stimulus_record_time;
post_analyzed = 20; %period to analyze after sound ends

%ILD
if ~isempty(x.stimulus_tags_ch2)
    rate.ILD =  x.stimulus_ch1(:, strcmp(' ILD dB', x.stimulus_tags_ch1));
    rate.ABI = rate.ILD(floor(length(rate.ILD)/2)+1);
    rate.attnILD(1,:) = x.stimulus_ch1(:, strcmp('attn(dB)', x.stimulus_tags_ch1));
    rate.attnILD(2,:) = x.stimulus_ch2(:, strcmp('attn(dB)', x.stimulus_tags_ch2));
end

% Determine variable
if max(attn) == min(attn)
    rate.attn = attn(1);
    attn_change = 0;
else
    rate.attn = attn;
    rate.y = attn;
    rate.var = 'Attn';
    attn_change = 1;
end

if max(dur) == min(dur)
    rate.dur = dur(1);
    dur_change = 0;
else
    rate.y = dur;
    rate.var = 'Duration';
    dur_change = 1;
end

if max(spkr) == min(spkr)
    rate.spkr = spkr(1);
    spkr_change = 0;
else
    rate.y = spkr;
    rate.var = 'Location' 
    spkr_change = 1;
    
    %rate.ILD
    %rate.freq
end
% Reformat data to be continuous

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified from E.Remington open_m_datafile
data=x.data;
stim = max(length(x.user_stimulus_desc_ch1), max(x.data(:,1)));
if x.data_format < 5 %non-continuous, not time corrected
    index_nn = find(data(:,4)> -1);
    data(index_nn,4)=data(index_nn,4)*1.024;
else %continuous and already time-corrected. Need to reformat
    %timing information to be trial-based
    if size(data,2) == 6
        behavior_data = data;
        data(:,5:6) = [];
    elseif size(data,2) == 5
        data(:,5) = [];
    end
    %     stimulus_deliveries_indices = find(data(:,3) == 7);
    stimulus_deliveries = data(data(:,3) == 7,:);
    for i = 1:length(stimulus_deliveries)
        %Convert from continuous to trial based timing
        data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
            stimulus_deliveries(i,2),4) = data(data(:,1) == stimulus_deliveries(i,1) & data(:,2) == ...
            stimulus_deliveries(i,2),4) - (stimulus_deliveries(i,4) - x.pre_stimulus_record_time*1000);
    end
end
reps=max(data(:,2));
rep_numbers = cell(1,stim);
if any(data(:,3)==7)
    for i = 1:max(data(:,1))
        behav_reps(i) = sum(data(:,1)==i & data(:,3)==7);
        rep_numbers{i} = data(data(:,1) == i & data(:,3) == 7,2);
    end
    data(data(:,3)>6,:) = [];
else
    behav_reps = ones(1,stim)*reps;
end
data((data(:,3)~=channel_no | data(:,4)==-1),:)=[]; %get rid of -1 and other channels

data(:,4)=round(data(:,4)/1000); %convert to ms
% if isempty(tWindow)
    tmin = pre_stim + 15;
    tmax = pre_stim + max(dur) + post_analyzed;
% else
%     tmin = tWindow{1}(1);
%     tmax = tWindow{1}(2);
% end
% tmax = pre_stim + max(stim_dur) + post_stim - 50; %users prefer that
% analysis window defaults to stimulus length
spont_rate = size(data((data(:,4) < pre_stim & data(:,4) > 0),:),1)/(sum(behav_reps)*pre_stim) * 1000; %spikes/sec
rate.spont = spont_rate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate rate, modified from zero_calc_rate
% accomodate varying duration stimuli
stim_dur = dur;
stim_overhang = tmax - (pre_stim + max(stim_dur));
if length(stim_dur)==1 | stim_overhang < 0 | tmin > min(stim_dur)+pre_stim
    reduced_data=data(find(data(:,4)>tmin & data(:,4)<tmax),:);
    for i=1:stim
        old_len = size(reduced_data,1);
        driven_spikes = find(reduced_data(:,1)==i);
        driven_data(i) = length(driven_spikes); %not right, too short, save last thing
        reduced_data(driven_spikes,:)=[];
        new_len = size(reduced_data,1);
        rates(i)=(old_len - new_len)/(reps*(tmax-tmin))*1000;%spikes/sec
    end
else
    reduced_data = data(find(data(:,4)>tmin),:);
    for i=1:stim
        old_len = size(reduced_data,1);
        driven_spikes = find(reduced_data(:,1)==i & reduced_data(:,4)<(pre_stim + stim_dur(i)+stim_overhang));
        driven_data_matrix = reduced_data(find(reduced_data(:,1)==i & reduced_data(:,4)<(pre_stim + stim_dur(i)+stim_overhang)),:);
        driven_data(i) = length(driven_spikes); 
        reduced_data(driven_spikes,:)=[];
        new_len = size(reduced_data,1);     
        rates(i)=(old_len - new_len)/(reps*(pre_stim + stim_dur(i) + stim_overhang - tmin))*1000; %spikes/sec
        for n=1:reps
            rate_rep(i,n) = length(find(driven_data_matrix(:,2)==n))/(pre_stim + stim_dur(i) + stim_overhang - tmin)*1000;
        end
    end
end %if length(stim_dur)
%Calc Stats
rate.se = zeros(1,stim);
for n = 1:stim
    rate.se(n) = std(rate_rep(n,:))/sqrt(reps);
    rate.sd(n) = std(rate_rep(n,:));
end

% r = mean(rate_rep');   %check things match up

rate.tuning = rates;



  