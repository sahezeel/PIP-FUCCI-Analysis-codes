clear
clc
close all

%%
%================User Input======================
data_links = {'\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\HCC KD 06\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\HCC KD 07\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\HCC KD 08\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\HCC KD 09\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\HCC KD 10\';};

%saves the merged mean of median (7th graph) to a specific location ONLY
mergedpath = '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\Merged ALL\HCC KD 06 7 8 9 10\';

%xticks
x_ticks = 2;

%window size for moving average filter to smooth the graphs
win = 7; %individual graph smoothning
win2 = 7; % median smoothening

%The color of the confidence intervals is decided by setting this value
%between 0 and 1. 0 is blank with no color and 1 is completely rich
patchSatVal = 0.1;
% winM = 5

%%
L = length(data_links);
Max_Num_frames = 0;
Gr = [];
Gr3 = [];
Re = [];
Re3 = [];
for i=1:L
    data_link = data_links{i};
    FindCrossOver(data_link);
    data = [data_link,'MEDIANS.mat'];
    load(data)%FRAMES, MEDIAN_G, MEDIAN_R
    frames{i} = FRAMES;
%     greens{i} = MEAN_G;
%     reds{i} = MEAN_R;
    greens{i} = MEDIAN_G;
    reds{i} = MEDIAN_R;
    
    sGr = size(Gr);
    sG = size(G);
    b = max(sGr(2),sG(2));
    Gr = [[Gr,zeros(abs([0,b]-sGr))];[G,zeros(abs([0,b]-sG))]];
    
    sGr3 = size(Gr3);
    sG3 = size(G3);
    b = max(sGr3(2),sG3(2));
    Gr3 = [[Gr3,zeros(abs([0,b]-sGr3))];[G3,zeros(abs([0,b]-sG3))]];
    
    sRe = size(Re);
    sR = size(R);
    b = max(sRe(2),sR(2));
    Re = [[Re,zeros(abs([0,b]-sRe))];[R,zeros(abs([0,b]-sR))]];
    
    sRe3 = size(Re3);
    sR3 = size(R3);
    b = max(sRe3(2),sR3(2));
    Re3 = [[Re3,zeros(abs([0,b]-sRe3))];[R3,zeros(abs([0,b]-sR3))]];
    
%     Gr = [Gr; G];
%     Gr3 = [Gr3; G3];
%     Re = [Re; R];
%     Re3 = [Re3; R3];

    if length(FRAMES)>Max_Num_frames
        Max_Num_frames = length(FRAMES);
        frames_with_Max_Num = FRAMES;
    end
    clear FRAMES MEDIAN_G MEDIAN_R
end
clear sG sGr sG sGr2 sG2 sRe sR sRe2 sR2 b;
frames = frames';
greens = greens';
reds = reds';

c3 = findXticks(frames_with_Max_Num,x_ticks);
x = (-c3*x_ticks:x_ticks:c3*x_ticks);


for i =1:L
    interm1 = zeros(Max_Num_frames,1);
    interm2 = zeros(Max_Num_frames,1);
    interm_g = greens{i};
    interm_r = reds{i};
    %interm_f = frames{i};
    %plot(interm_f,interm_g,'g',interm_f,interm_r,'r')
    interm_g = interm_g';
    interm_r = interm_r';
    interm1((length(interm1)-length(interm_g))/2+1:end-(length(interm1)-length(interm_g))/2) = interm_g;
    interm2((length(interm2)-length(interm_g))/2+1:end-(length(interm2)-length(interm_g))/2) = interm_r;
    for j=1:Max_Num_frames
        r(i,j) = interm2(j,1);
        g(i,j) = interm1(j,1);
    end
    %plot(frames_with_Max_Num,r(i,:),'r',frames_with_Max_Num,g(i,:),'g')
    clear interm1 interm2
end
median_r = MEDIAN(r);
median_g = MEDIAN(g);

%plot green with outliers for all ROIs 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(Gr,1)
    xplot = frames_with_Max_Num;
    yplot = Gr(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [144/255,238/255,144/255] ,'LineWidth',1)
end
% plot(Max_Num_frames,smooth(MEAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',2)
plot(frames_with_Max_Num,smooth(median_g,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [mergedpath,'GreenChannel_all_ROIs_w_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off


%plot green without outliers for all ROIs 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(Gr3,1)
    xplot = frames_with_Max_Num;
    yplot = Gr3(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [144/255,238/255,144/255] ,'LineWidth',1)
end
% plot(frames_with_Max_Num,smooth(MEAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',2)
plot(frames_with_Max_Num,smooth(median_g,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [mergedpath,'GreenChannel_all_ROIs_wo_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

%plot red with outliers for all ROIs 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(Re,1)
    xplot = frames_with_Max_Num;
    yplot = Re(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [220/255,20/255,60/255] ,'LineWidth',1)
end
% plot(frames_with_Max_Num,smooth(MEAN_G,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',2)
plot(frames_with_Max_Num,smooth(median_r,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [mergedpath,'RedChannel_all_ROIs_w_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off


%plot red without outliers for all ROIs 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(Re3,1)
    xplot = frames_with_Max_Num;
    yplot = Re3(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [220/255,20/255,60/255] ,'LineWidth',1)
end
% plot(frames_with_Max_Num,smooth(MEAN_G,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',2)
plot(frames_with_Max_Num,smooth(median_r,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [mergedpath,'RedChannel_all_ROIs_wo_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

% Median of Median graph. This is the final graph we want

FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(frames_with_Max_Num,g,{@median,@std},'lineprops',{'-','color',[0/255, 128/255, 0/255], 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);

shadedErrorBar(frames_with_Max_Num,r,{@median,@std},'lineprops',{'-','color', [220/255,20/255,60/255], 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);
% hold on

% plot(frames_with_Max_Num,median_r,'r',frames_with_Max_Num,median_g,'g')
% plot(frames_with_Max_Num,smooth(median_r,winM,'moving'),'color',[220/255,20/255,60/255],'LineWidth',6)
% hold on
% plot(frames_with_Max_Num,smooth(median_g,winM,'moving'),'color',[0/255,128/255,0/255],'LineWidth',6)

xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [mergedpath, 'Median of Medians.fig'];
savefig(FigH,savefile) %SA added this to save file as .fig
hold off;

close all
save ([mergedpath,'Entire workspace all ROIs.mat']);


%%
function MEDIAN_R = MEDIAN(R)
sR = size(R);%86,471
for j3=1:sR(2)
    a = R(:,j3);
    an = nonzeros(a);
    if ~isempty(an)
        MEDIAN_R(j3) = median(an);
    else
        MEDIAN_R(j3) = 0;
    end
end
end

function c3 = findXticks(FRAMES,x_ticks)
c = (length(FRAMES) - 1)/2;
c1 = 0;
c2 = -1;
while c2<0
c1 = c1 + 1;
c2 = c1*x_ticks - c;
end
c3 = c1 - 1;
end