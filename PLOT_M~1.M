clear
clc
close all

%%
%================User Input======================
data_links = {'\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\16\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\17\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\18\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\19\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\20\';};

%saves the merged mean of median (7th graph) to a specific location ONLY
mergedpath = '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\Merged Mean\KD_Palblo 16-20\';

%xticks
x_ticks = 2;

%window size for moving average filter to smooth the graphs
win = 7; %individual graph smoothning
win2 = 7; % median smoothening

%The color of the confidence intervals is decided by setting this value
%between 0 and 1. 0 is blank with no color and 1 is completely rich
patchSatVal = 0.1;
% winM = 5

% number of colors for Jet color map
n = length(data_links);

%%
L = length(data_links);
Max_Num_frames = 0;
Gr = {};
Re = {};
for i=1:L
    data_link = data_links{i};
    PlotRawGreenRed(data_link);
    data = [data_link,'MEDIANS_raw.mat'];
    load(data)%FRAMES, MEDIAN_G, MEDIAN_R
    frames{i} = FRAMES;
%     greens{i} = MEAN_G;
%     reds{i} = MEAN_R;
    greens{i} = MEDIAN_G;
    reds{i} = MEDIAN_R;
    Gr{i} = G;
    Re{i} = R;
    if length(FRAMES)>Max_Num_frames
        Max_Num_frames = length(FRAMES);
        frames_with_Max_Num = FRAMES;
    end
    clear FRAMES MEDIAN_G MEDIAN_R
    pos_of_fwdslashes = strfind(convertCharsToStrings(data_link),'\');
    leg(((i-1)*4)+4) = convertCharsToStrings(data_link(pos_of_fwdslashes(end-1)+1:pos_of_fwdslashes(end)-1));
end

%-------plot green with median and shaded error bars for all ROIs---------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
ii = (1:n)';
c = num2cell(jet(length(ii)));
for i = 1:L
    shadedErrorBar(frames{i},Gr{1,i},{@median,@std},'lineprops',{'-', 'MarkerFaceColor', cell2mat(c(i,:)), 'color', cell2mat(c(i,:))},'patchSaturation',patchSatVal);
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(cellstr(leg));
savefile = [mergedpath,'GreenChannel_shaded_all_ROIs.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------% 

%--------plot red with median and shaded error bars for all ROIs----------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
ii = (1:n)';
c = num2cell(jet(length(ii)));
for i = 1:L
    shadedErrorBar(frames{i},Re{1,i},{@median,@std},'lineprops',{'-', 'MarkerFaceColor', cell2mat(c(i,:)), 'color', cell2mat(c(i,:))},'patchSaturation',patchSatVal);
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(cellstr(leg));
savefile = [mergedpath,'RedChannel_shaded_all_ROIs.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------%
