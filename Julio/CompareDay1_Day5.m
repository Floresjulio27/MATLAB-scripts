zap;

%% Load data
Day_1 = load('Fig1_S2_Day 1_Data.mat'); data_Day_1 = Day_1.data;
Day_5 = load('Fig1_S2_Day 5_Data.mat'); data_Day_5 = Day_5.data;

summaryFigureOverlay = figure;
sgtitle('Stimulus-evoked responses: Day_1 vs Day_5')
%% Plotting 

%set colors 
gfp_base_color = [0 0 0];
gfp_Day_5_color = [0.55, 0.45, 0.70];
cbv_base_color = [0.30, 0.30, 0.30];
cbv_Day_5_color = [0.70, 0.40, 0.10];

% Contra stim PFC fiber GFP (Right hemisphere, left whisker stim)
ax1 = subplot(2,2,1);

% GFP Day_1
x = data_Day_1.Contra.mean_timeVector(:)';
y = data_Day_1.Contra.P_NEGFP_Mean(:)';
e = data_Day_1.Contra.P_NEGFP_SEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_base_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', gfp_base_color, 'LineWidth', 2);

% GFP Day_5 
x2 = data_Day_5.Contra.mean_timeVector(:)';
y2 = data_Day_5.Contra.P_NEGFP_Mean(:)';
e2 = data_Day_5.Contra.P_NEGFP_SEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day_5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day_5_color, 'LineWidth', 2);

title('Contra stim PFC fiber GFP')
ylabel('\DeltaF/F SST Ca2+ activity (%)')
xlabel('Time (s)')
legend([p1 p2],'Day 1', 'Day 5')
axis square; 
axis tight;
xlim([-5 15]);
ylim([-0.4 1.5]);

% Contra stim PFC fiber GFP (Right hemisphere, left whisker stim)
ax2 = subplot(2,2,2);

%CBV Day_1
x3 = data_Day_1.Contra.mean_timeVector(:)';
y3 = data_Day_1.Contra.P_NECBV_Mean(:)';
e3 = data_Day_1.Contra.P_NECBV_SEM(:)';
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], cbv_base_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p3 = plot(x3, y3, 'Color', cbv_base_color, 'LineWidth', 2);

%CBV Day_5
x4 = data_Day_5.Contra.mean_timeVector(:)';
y4 = data_Day_5.Contra.P_NECBV_Mean(:)';
e4 = data_Day_5.Contra.P_NECBV_SEM(:)';
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], cbv_Day_5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p4 = plot(x4, y4, 'Color', cbv_Day_5_color, 'LineWidth', 2);

title('Contra stim PFC fiber CBV')
ylabel('\DeltaF/F CBV (%)')
xlabel('Time (s)')
legend([p3 p4],'Day 1', 'Day 5')
axis square; 
axis tight;
xlim([-5 15]);
ylim([-0.5 4.0])


% %CBV
% ax1 = subplot(2,3,1);
% plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean + data.Contra.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.Contra.mean_timeVector,data.Contra.P_AChCBV_Mean - data.Contra.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% title('Contra stim S1BF fiber')
% ylabel('\DeltaF/F S1BF CBV (%)')
% ax1.YLim = [-1 2];
% 
% %GFP
% yyaxis right
% plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean + data.Contra.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.Contra.mean_timeVector,data.Contra.P_AChGFP_Mean - data.Contra.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% ylabel('\DeltaF/F S1BF SST Ca2+ acitvity (%)')
% 
% 
% 
% %%%% other
% plot(data_baseline.Contra.mean_timeVector, data_baseline.Contra.P_NEGFP_Mean, '-', 'Color', gfp_base_color, 'LineWidth', 2)
% plot(data_baseline.Contra.mean_timeVector, data_baseline.Contra.P_NEGFP_Mean + data_baseline.Contra.P_NEGFP_yCI95, '-', 'Color', gfp_base_color, 'LineWidth', 0.5)
% plot(data_baseline.Contra.mean_timeVector, data_baseline.Contra.P_NEGFP_Mean - data_baseline.Contra.P_NEGFP_yCI95, '-', 'Color', gfp_base_color, 'LineWidth', 0.5)
% plot(data_did4.Contra.mean_timeVector, data_did4.Contra.P_NEGFP_Mean, '--', 'Color', gfp_did4_color, 'LineWidth', 2)
% plot(data_did4.Contra.mean_timeVector, data_did4.Contra.P_NEGFP_Mean + data_did4.Contra.P_NEGFP_yCI95, '--', 'Color', gfp_did4_color, 'LineWidth', 0.5)
% plot(data_did4.Contra.mean_timeVector, data_did4.Contra.P_NEGFP_Mean - data_did4.Contra.P_NEGFP_yCI95, '--', 'Color', gfp_did4_color, 'LineWidth', 0.5)
% ylabel('\DeltaF/F SST Ca2+ activity (%)')
% xlabel('Time (s)'); xlim([-5 15]); axis square
