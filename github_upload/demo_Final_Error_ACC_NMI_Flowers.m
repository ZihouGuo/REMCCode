clear
% close all
%% ��������
load Error_NonKtSVD_Soft_Flower1;
X= 1:length(out.error1);
Error = max(out.error1,out.error2);
ACC   = out.MeanACC;
NMI   = out.MeanNMI;


%% ��ͼ����������y������
fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
%�������
yyaxis left
plot(X, Error,'b--', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','b');   
h  = ylabel('Objective function value')
set(h,'fontsize',20,'color','b');
% xlabel('Iterations')
%���ÿ̶�
% axis([1 5 98.8 99.4]);
% set(gca,'YTick',[98.8 99 99.2 99.4]);
%�����Ҳ�
yyaxis right
plot(X, ACC,'m--', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','r');
hold on
plot(X, NMI,'r--', 'LineWidth',1.5, 'MarkerSize',6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','r');

h = ylabel('NMI/ACC')
set(h,'fontsize',20,'color','b');
h = title('(e) Flower')
set(h,'fontsize',20);
h = xlabel('Iterations');
set(h,'fontsize',20);
h = legend('Error', 'ACC', 'NMI');
set(h,'fontsize',20);
xlim([1 length(ACC)])
%���ÿ̶�
% axis([1 5 0.984 0.99]);
% set(gca,'YTick',[0.984 0.986 0.988 0.99]);
%������
grid on
grid minor
