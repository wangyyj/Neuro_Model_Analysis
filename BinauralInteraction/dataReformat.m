% Time-corrects data to be trial-based
% YYWang 6/2015
function data = dataReformat(data, preStim)
pre_stimulus_record_time = preStim;


%edited from E.Remington open_m_datafile
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

