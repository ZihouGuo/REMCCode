clear;
clc;
figure('position',[150,100,900,550])%ȷ��ͼƬ��λ�úʹ�С��[x y width height]
%׼������
%Y=[8.92,10.23,13.93,12.61,17.07;5.67,6.70,8.24,10.03,9.79;6.78,7.97,0,10.17,11.40;6.98,8.68,0,9.85,11.64;5.58,7.24,6.91,7.65,9.75];
Y=[8.4,8.6;6.21,8.04;7.32,10.1];
X=1:3;
 %����5����״ͼ�����1
h=bar(X,Y,1);      
 %�޸ĺ��������ơ�����
set(gca,'XTickLabel',{'Reuters','Cifar10','Cifar100'},'FontSize',14,'FontName','Times New Roman');
% ����������ɫ,��ɫΪRGB��ԭɫ��ÿ��ֵ��0~1֮�伴��
set(h(1),'FaceColor',[130,205,208]/255)     
set(h(2),'FaceColor',[54,49,143]/255)    
% set(h(1),'FaceColor',[228,0,127]/255)    
%  set(h(2),'FaceColor',[234,186,202]/255)    
%  set(h(2),'FaceColor',[135,20,104]/255)   
ylim([0,17.5]);      %y��̶�
%�޸�x,y���ǩ
ylabel('\fontname{Times New Roman}\fontsize{18}Time (log_2)');
%xlabel('\fontname{����}\fontsize{14}��ͬ��'); 
%�޸�ͼ��
legend({'\fontname{Times New Roman}IMVC-CBG','\fontname{Times New Roman}TDASC'},'FontSize',13,'Location','NorthWest');
