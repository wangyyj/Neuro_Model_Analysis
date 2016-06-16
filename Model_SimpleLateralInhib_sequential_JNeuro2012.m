% SimpleLateralInhib_sequential
% Inhibitory weights are gaussian, population histogram is chi^2
% Surround inhibitory weights are integrated across population
% Inhibitory weights are scaled by tuning curves - increasing widths with
% laterality

% Published Wang et al. 2012 Journal of Neuroscience
% Wang Y, Shanbhag SJ, Fischer BJ, Peña JL. (2012) Population-Wide Bias of Surround suppression in auditory spatial receptive fields of the owl's midbrain. J Neurosci. 32(31):10470-8. 
%
% Yunyan Wang 2012
% Wangyyj@gmail.com

clear all;

%% Network

%Number of neurons
N=501;
%N=1001;
%Centers (deg)
center=chi2rnd(5,1,N)'; %Chi2 cell population
center=center*5.5-20;
center=sort(center);
%load(    )%not randomized center every time script is run

%Weight matrix
% %Gaussian Weights for anatomical distance, default sigma=0.4
% gau=normpdf(-1:2/(N-1):1,0,0.4);
% %gau=normpdf(-3:6/N:3)-((normpdf(-20:40/N:20))-0.01); %Gaussian with notch
% gau=gau-min(gau);
% gau(250)=0;
% gau=[gau(1:250) fliplr(gau(1:250))];
a_x=-250:250;
%a_x=-500:500;
a_halfwidth=70;
a_sigma=a_halfwidth/2.35482;
gau=normpdf(a_x,0,a_sigma);
gau=gau-min(gau);
gau(250)=0;
gau=gau*(1/max(gau));

[temp mid]=max(gau);
W=zeros(N,N+2*mid);

place=1;
for n=1:N
    W(n,place:place+N-1)=gau;
    place=place+1;
end
W2=W(:,mid:(mid+N)-1);

% Sequential - halfwidth ranges from 23-33 from 0 to 40
x=0:ceil(max(center)-min(center)); 
for n=1:N
    theta=abs(center(n));
    halfwidth(n)=0.25*theta+23; %increasing halfwidth(seq)
    %halfwidth(n)=23; %no change in tuning width(SWN)
    sigma=halfwidth(n)/2.35482;
    y=normpdf(x,0,sigma);
    y=y*1/max(y);
    tuningCurve(n,:)=y;
end

%Calculate activation to centers, each center(n) is preferred direction
for n=1:N 
    for m=1:N
        y=tuningCurve(n,:);
        delta=abs(center(n)-center(m));
        activation(n,m)=interp1(x,y,delta);
    end
end

% Modify suppression over anatomic space by tuning curve info
% When speaker at center(400) is activated, deg of activation of  the surround
%   is activation(:,400)
for n=1:N
    W3(n,:)=W2(n,:).*activation(:,n)';
end
W3(W3==NaN)=0;
% %Calculate frontal vs peri Suppression
for n=1:N
    frontal(n)=sum(W3(n,1:n));
    peri(n)=sum(W3(n,n+1:end));
end
figure('Color',[1 1 1]);
plot(center,frontal);
hold on;
plot(center,peri, 'r');
title('Frontal vs Peripheral Inhibitory Drive')
xlabel('Preferred Tuning in Az (deg)')
ylabel('Inhibitory Drive')
legend('Frontal Suppression', 'Peripheral Suppression')

%plot percent difference
dif=100*(frontal-peri)./frontal;
figure('Color',[1 1 1]);
plot(center,dif);
title('% Difference between Frontal and Peripheral Inhibitory Drive')
xlabel('Preferred Tuning in Az (deg)')
ylabel('% Difference in Inhibitory Drive')
xlim([-20 70])
ylim([0 110])


% %groups 0-60 at 15deg interval
% place=0;
% for n=1:5
%     deg(n)=find(center>place,1,'first');
%     place=place+15;
% end
% %mean suppression for groups
% place=1;
% for n=1:4
%     groups(n,1)= mean(frontal(deg(place):deg(place+1)));
%     groups(n,2)= mean(peri(deg(place):deg(place+1)));
%     place=place+1;
% end    
% %Plot mean suppression for groups
% norm=max(max(groups));
% groups=groups/norm;
% figure; plot(groups(:,1));
% hold on; plot(groups(:,2),'r')
% title('% Mean Suppression for Groups')
% xlabel('Groups)')
% ylabel('Normalized Total Suppression')
% legend('Frontal Suppression', 'Peripheral Suppression')
% 
% %half-width integration dist for each cell on each side
% 
% for n=1:N   %
%     dist(n,1)= find(W3(n,:)>max(gau)/2,1,'first');
%     dist(n,2)= find(W3(n,:)>max(gau)/2,1,'last');
%     dist(n,3)= center(dist(n,2))-center(dist(n,1)); %total
%     dist(n,4)= center(n)-center(dist(n,1)); %frontal
%     dist(n,5)= center(dist(n,2))-center(n); %peripheral
% end
% 
% figure;
% plot(center,dist(:,3))
% hold on
% plot(center,dist(:,4),'r')
% plot(center,dist(:,5),'g')
% xlim([0,90])
% legend('Total', 'Frontal', 'Peripheral');
% xlabel('Preferred Azimuth')
% ylabel('Integration Distance (degs)')

% %Sequential Stimulation
% for n=1:N
%     short_dist(n,1)=length(find(center<center(n)& center>center(n)-10));
%     short_dist(n,2)=length(find(center>center(n)& center<center(n)+10));
%     short_dist(n,3)=sum(short_dist(n,1:2));   %Tot # cells in 10deg surround
%     short_dist(n,4)=sum(gau(250-short_dist(n,1):250)); 
%     short_dist(n,5)=sum(gau(251:251+short_dist(n,2)));
% end
% short_dist(:,6)= short_dist(:,4)+short_dist(:,5); %Tot Suppression in 10deg Surround
% short_dist(:,7)= sum(W2'); %Tot Suppression with entire surround activated
% RF_dist(:,1)=300-short_dist(:,6); %Sequential
% RF_dist(:,2)=300-short_dist(:,7); %SWN
% RF_dist= RF_dist/max(max(RF_dist));
% figure;
% plot(center,RF_dist(:,1))
% hold on
% plot(center,RF_dist(:,2),'r')
% title('RF size with Sequential / SWN Stimulation')
% legend('Sequential Stimulation','SWN Stimulation')
% xlabel('Preferred Azimuth (degs)')
% ylabel('Relative Receptive Field Size')
% xlim([0 80])


% % Find average 1/2 width for groups
% place=1;
% for n=1:4
%     halfwidth_groups(n,1)=mean(halfwidth(deg(place:place+1),1)); %Frontal
%     halfwidth_groups(n,2)=mean(halfwidth(deg(place:place+1),2)); %Peripheral   
%     halfwidth_groups(n,3)=sum(halfwidth_groups(n,1:2)); %sum
%     place=place+1;
% end
% %Plot half-width distance
% figure; 
% plot(halfwidth_groups(:,1));
% hold on; 
% plot(halfwidth_groups(:,2),'r')
% plot(halfwidth_groups(:,3),'k')
% title('% Integration Widths (1/2 max of suppression)')
% xlabel('Groups)')
% ylabel('Integration Distance (degs)')
% legend('Frontal Suppression', 'Peripheral Suppression')

% %%temp loop for data analysis - 1/2 width linear suppression
% for n=1:3
%     widths(n,1)=length(find(data(n,1:11)>max(data(n,1:11))/2));
%     widths(n,2)=length(find(data(n,11:end)>max(data(n,11:end))/2));
%     widths(n,3)=sum(widths(n,1:2));
% end
%         
    
    
    
    
    
    
    

