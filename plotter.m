clear 
clc
close all

load("Copy_of_monte_carlo_data.mat")

subset1 = results(1:250);
subset2 = results(251:500);
subset3 = results(501:750);
subset4 = results(751:1000);

data = zeros(250,4);

for i=1:250
    data(i,4) = subset1(i).ss_pct_mpc_improvement;
    data(i,3) = subset2(i).ss_pct_mpc_improvement;
    data(i,2) = subset3(i).ss_pct_mpc_improvement;
    data(i,1) = subset4(i).ss_pct_mpc_improvement;
end

figure_configuration_IEEE_standard

bh = boxplot(data,'Labels',{'25%','50%','75%','100%'},'BoxStyle','outline','Jitter',0.25)
xlabel('Network Connectivity Level')
ylabel(sprintf('MPC Cost Improvement \nover Model-Free'))
ytickformat('percentage')
set(bh,'LineWidth', 2);
yt = [0,0.05,0.1];
ytl = compose('%.0f%%', yt);
set(gca,'YTick',yt);
%%
clc
pct_shift_mpc = zeros(250*20,4);
pct_shift_mf = zeros(250*20,4);

for i=1:250
    pct_shift_mpc(20*(i-1)+1:20*(i-1)+20,4) = subset1(i).pct_shift_mpc;
    pct_shift_mpc(20*(i-1)+1:20*(i-1)+20,3) = subset2(i).pct_shift_mpc;
    pct_shift_mpc(20*(i-1)+1:20*(i-1)+20,2) = subset3(i).pct_shift_mpc;
    pct_shift_mpc(20*(i-1)+1:20*(i-1)+20,1) = subset4(i).pct_shift_mpc;
    pct_shift_mf(20*(i-1)+1:20*(i-1)+20,4) = subset1(i).pct_shift_mf;
    pct_shift_mf(20*(i-1)+1:20*(i-1)+20,3) = subset2(i).pct_shift_mf;
    pct_shift_mf(20*(i-1)+1:20*(i-1)+20,2) = subset3(i).pct_shift_mf;
    pct_shift_mf(20*(i-1)+1:20*(i-1)+20,1) = subset4(i).pct_shift_mf;

end

% data = {pct_shift_mpc,pct_shift_mf};
% fig = boxplot(pct_shift_mpc);
% set(fig,'LineWidth', 2);
% % fig.axis.YAxis.Scale ="log";
% % ytickformat(fix.axis,'percentage')
% ylim([-1 93])
% % 
% yt = get(gca, 'YTick');
% yt = [0, 5, 10, max(pct_shift_mpc,[],'all')]
% ytl = compose('%.0f%%', yt)
% set(gca, 'YTick',yt, 'YTickLabel',ytl)
% 
% breakyaxis([11,92],0.02)




% Step 1: Prepare the data
interleavedData = zeros(size(pct_shift_mpc, 1), size(pct_shift_mpc, 2) * 2);
for i = 1:size(pct_shift_mpc, 2)
    interleavedData(:, (i-1)*2 + 1) = pct_shift_mpc(:, i);
    interleavedData(:, (i-1)*2 + 2) = pct_shift_mf(:, i);
end

% Step 2: Plot with boxplot
f = figure; % Create a new figure
% Calculate positions for each box. Adjust the spacing as needed.
positions = 1:size(interleavedData, 2);
spacing = 1; % Change this value to adjust the space between pairs
% for i = 2:2:size(interleavedData, 2)
%     positions(i:end) = positions(i:end) + spacing;
% end
positions = [1,2,4,5,7,8,10,11]

bh1 = boxplot(gca,pct_shift_mpc,'Positions',[1,3.5,6,8.5],'Colors','b','Symbol','+b','Widths',0.9,'Jitter',0.25);
hold on
bh2 = boxplot(gca,pct_shift_mf,'Positions',[2,4.5,7,9.5],'Colors','r','Widths',0.8,'Jitter',0.25);

hold on
% bh = boxplot(interleavedData, 'Positions', positions, 'ColorGroup', repmat([1, 2], 1, size(pct_shift_mpc, 2)), 'Colors', ['b', 'r']);

% Customizing the appearance
set(bh1, 'LineWidth', 1.5); % Set the line width
set(bh2, 'LineWidth', 1.5);
% legend('MPC', 'MF'); % Add a legend to distinguish between MPC and MF

% Step 3: Customize labels
% Adjust the x-tick labels to reflect the new grouping
set(gca, 'XTick', [1.5 4, 6.5, 9], 'XTickLabel', {'25%', '50%', '75%', '100%'});
% 
ylim([-1 93.5]); % Adjust as necessary
% 
% % Formatting Y-axis labels as percentage
% yt = get(gca, 'YTick');
% ytl = compose('%.0f%%', yt);
% set(gca, 'YTickLabel', ytl);
% 
% yt = [0,5,10,92.5]
% ytl = compose('%.0f%%', yt);
% set(gca, 'YTickLabel', ytl);
% % Add your breakyaxis call here
xlh = xlabel('Network Connectivity Level')
% xlh.Position(2) = xlh.Position(2)-0.5
ylh = ylabel('Opinion Shift')
ylh.Position(1)=ylh.Position(1)-0.25;
breakinfo = breakyaxis([11,90.5],0.02);
% asdf = ylabel(breakinfo.lowAxes,"Opinion Shift")
% yt = get(gca, 'YTick');
% yt = [0,5,10,92.5]
% ytl = compose('%.0f%%', yt);
% set(gca, 'YTickLabel', ytl);
%%
yt = get(breakinfo.lowAxes,'YTick');
ytl = compose('%.0f%%', yt);
set(breakinfo.lowAxes,'YTickLabel',ytl);
% yt = get(breakinfo.highAxes,'YTick');
yt = [91 92];
ytl = compose('%.0f%%', yt);
set(breakinfo.highAxes,'YTick',yt);
set(breakinfo.highAxes,'YTickLabel',ytl,'fontsize',4*8);

h1 = line(NaN, NaN, 'Marker', 'none', 'Color', 'b');
h2 = line(NaN, NaN, 'Marker', 'none', 'Color', 'r');

% Creating the legend
legend([h1, h2], {'MPC', 'Model-Free'}, 'AutoUpdate','off');
% ylabel('Opinion Shift')
ax = get(gcf,'children');
ind = find(isgraphics(ax,'Legend'));
set(gcf,'children',ax([ind:end,1:ind-1]))