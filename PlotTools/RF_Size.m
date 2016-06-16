% RF_Size.mat  plots bar graph to show changes in RF size. 
% data: two columns; col1= Seq RF pixel #; col2=RevRF RF pixel #
% Plots bar, scatter and lines connecting scatter from 2 columns

dotaxis=ones(length(data),2);
dotaxis(:,2)=2;

%bar graph
bar(1,mean(data(:,1)),'w');
hold on
bar(2,mean(data(:,2)),'w');

%scatter data points
scatter(dotaxis(:,1),data(:,1),'r')
scatter(dotaxis(:,2),data(:,2),'r')

std1=std(data(:,1));
std2=std(data(:,2));


%errorbar
errorbar(1,mean(data(:,1)),std1,'xk')
errorbar(2,mean(data(:,2)),std2,'xk')


%lines
for count=1:length(data)
 plot(data(count,:),'-.k')
 hold on
end


ylabel('Receptive Field Size (pixel)')
title ('Reduction in Receptive Field Size with Overlapping Stimuli')
