clear
%______________________________________
% run lines below before each analysis
global g_zero
tab_index = g_zero.open_tabs;
cur_tab = g_zero.tabs{1,tab_index};
y_data= cur_tab.y_data{1,2};
x_data= cur_tab.x_data{1,2};
%______________________________________

%Plot SRF
set(0,'DefaultAxesFontSize',9)
zero_plot_srf2;
cpos = get(c,'Position');
cpos(3) = 0.3*cpos(3);
axpos = get(gca,'Position');
set(c,'Position',cpos)
set(gca,'Position',axpos)
ylabel(c,'Spikes (Hz)')
axis tight

%plot RASTER
zero_plot_raster_offline;
set(F, 'Position', [100, 100, 370, 220]);
% ILD only
if length(y_data)<10
    set(gca,'YTick',[1:9],'YTickLabel',[-20:5:20])
else
    set(gca,'YTick',[1:2:13],'YTickLabel',[-30:10:30])
end
ylabel('ILD (dB)')
xlim([100 600])

%plot Rate
y_data= cur_tab.y_data{1,2};
x_data= cur_tab.x_data{1,2};
spontR= cur_tab.spont_rate;
init_fig_arial
plot(x_data,y_data,'k')
hold on
plot([1 length(y_data)],[spontR spontR],'Color',[0.7 0.7 0.7],'LineStyle','--')
box off
%ILD only
set(F, 'Position', [100, 100, 370, 220]);
xlabel('ILD (dB)')
ylabel('Response (spk/s)')
if length(y_data)<10
    set(gca, 'Xtick',[1:9],'XTickLabel',[-20:5:20])
else
    set(gca, 'Xtick',[1:2:13],'XTickLabel',[-30:10:30])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot raster

title('M6x1586   Stimulus: Twit-peep Call')
set(gca,'YTick',[])
ylabel('Calls','FontSize', 13)
set(gca,'YTick',[1:2:20],'YTickLabel',[1:2:20])
set(gca,'YTick',[1:4:40],'YTickLabel',[1:4:40])

