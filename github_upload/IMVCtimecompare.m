clear;
clc;
figure('position',[150,100,900,550])%确定图片的位置和大小，[x y width height]
%准备数据
%Y=[8.92,10.23,13.93,12.61,17.07;5.67,6.70,8.24,10.03,9.79;6.78,7.97,0,10.17,11.40;6.98,8.68,0,9.85,11.64;5.58,7.24,6.91,7.65,9.75];
Y=[8.4,8.6;6.21,8.04;7.32,10.1];
X=1:3;
 %画出5组柱状图，宽度1
h=bar(X,Y,1);      
 %修改横坐标名称、字体
set(gca,'XTickLabel',{'Reuters','Cifar10','Cifar100'},'FontSize',14,'FontName','Times New Roman');
% 设置柱子颜色,颜色为RGB三原色，每个值在0~1之间即可
set(h(1),'FaceColor',[130,205,208]/255)     
set(h(2),'FaceColor',[54,49,143]/255)    
% set(h(1),'FaceColor',[228,0,127]/255)    
%  set(h(2),'FaceColor',[234,186,202]/255)    
%  set(h(2),'FaceColor',[135,20,104]/255)   
ylim([0,17.5]);      %y轴刻度
%修改x,y轴标签
ylabel('\fontname{Times New Roman}\fontsize{18}Time (log_2)');
%xlabel('\fontname{宋体}\fontsize{14}不同组'); 
%修改图例
legend({'\fontname{Times New Roman}IMVC-CBG','\fontname{Times New Roman}TDASC'},'FontSize',13,'Location','NorthWest');
