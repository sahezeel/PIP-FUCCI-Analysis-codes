clear
clc
close all

%%
data_links.name{1} = 'A2 CT 1-5';
data_links.treatment{1} = {'\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\A2 KD 55\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\A2 KD 56\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\A2 KD 57\';
    '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\A2 KD 59\';};

% data_links.name{2} = 'KEN 16-20';
% data_links.treatment{2} = {'Z:\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\04232021 H441 A2 W D K IND RO top RO Palblo bot 5 pos 8 well\2021-04-23\Analysis using red drop codes\A2_KEN_16\';
%     'Z:\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\04232021 H441 A2 W D K IND RO top RO Palblo bot 5 pos 8 well\2021-04-23\Analysis using red drop codes\A2_KEN_17\';
%     'Z:\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\04232021 H441 A2 W D K IND RO top RO Palblo bot 5 pos 8 well\2021-04-23\Analysis using red drop codes\A2_KEN_18\';
%     'Z:\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\04232021 H441 A2 W D K IND RO top RO Palblo bot 5 pos 8 well\2021-04-23\Analysis using red drop codes\A2_KEN_19\';
%     'Z:\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\04232021 H441 A2 W D K IND RO top RO Palblo bot 5 pos 8 well\2021-04-23\Analysis using red drop codes\A2_KEN_20\';};

% data_links.name{3} = 'ken';
% data_links.treatment{3} = {'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Analysis using Zhouping codes\Stuck in G1\H1975 KD 23\';
%     'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Analysis using Zhouping codes\Stuck in G1\H1975 KD 24\';};

% This is 
% data_links(4).name(4) = 'ken';
% data_links.treatment{4} = {'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Analysis using Zhouping codes\Stuck in G1\H1975 KD 21\';
%     'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Analysis using Zhouping codes\Stuck in G1\H1975 KD 22\';};


%saves the merged mean of median (7th graph) to a specific location ONLY
mergedpath = '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\UM Core and Indika Confocal Data\082020 HCC H44 LC3 a2 H1975 PIP Ct KD RNAiMax\Analysis\Analysis using New codes 2021 May\Merged ALL\Stuck in G1\H441 55 56 57 59 stuck in g1\';

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
L = length(data_links.treatment);
Gr = {};
Re = {};
for i=1:L
    data_link = data_links.treatment{i};
    name = data_links.name{i};
    find_median_Mitosis_stuck_in_G1(data_link, mergedpath, name);
    sname = ['MEDIANS_stuck_in_G1_', name, '.mat'];
    data = [mergedpath,sname];
    load(data)%FRAMES, MEDIAN_G, MEDIAN_R
    fr{i} = frames_with_Max_Num;
%     greenMean{i} = mean_g;
%     redsMean{i} = mean_r;
    greenMed{i} = median_g;
    redMed{i} = median_r;
    Gr{i} = g;
    Re{i} = r;
%     x{i} = x;
%     if length(FRAMES)>Max_Num_frames
%         Max_Num_frames = length(FRAMES);
%         frames_with_Max_Num = FRAMES;
%     end
    clear median_g median_r r g frames_with_Max_Num;
%     pos_of_fwdslashes = strfind(convertCharsToStrings(data_link),'\');
%     leg(((i-1)*4)+4) = convertCharsToStrings(data_link(pos_of_fwdslashes(end-1)+1:pos_of_fwdslashes(end)-1));
end

%%
%-------plot green with median and shaded error bars for all ROIs---------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
dispLeg = [];
ii = (1:L)';
c = num2cell(jet(length(ii)));
for i = 1:L
    Hg(i) = shadedErrorBar(fr{i},Gr{1,i},{@median,@std},'lineprops',{'-', 'MarkerFaceColor', cell2mat(c(i,:)), 'color', cell2mat(c(i,:))},'patchSaturation',patchSatVal);
    dispLeg = [dispLeg, Hg(i).mainLine];
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(dispLeg, cellstr(data_links.name));
savefile = [mergedpath,'GreenChannel_shaded_all_Treatments.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------% 


%--------plot red with median and shaded error bars for all ROIs----------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
dispLeg = [];
ii = (1:L)';
c = num2cell(jet(length(ii)));
for i = 1:L
    Hr(i) = shadedErrorBar(fr{i},Re{1,i},{@median,@std},'lineprops',{'-', 'MarkerFaceColor', cell2mat(c(i,:)), 'color', cell2mat(c(i,:))},'patchSaturation',patchSatVal);
    dispLeg = [dispLeg, Hr(i).mainLine];
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(dispLeg, cellstr(data_links.name));
savefile = [mergedpath,'RedChannel_shaded_all_Treatments.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------%





%----------plot green with only median of medians for all ROIs------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
ii = (1:L)';
c = num2cell(jet(length(ii)));
for i = 1:L
   gH(i) = plot(fr{i},smooth(greenMed{i},win2,'moving'),'color',cell2mat(c(i,:)),'LineWidth',6);
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(cellstr(data_links.name));
savefile = [mergedpath,'GreenChannel_median_of_medians_all_Treatments.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------% 


%-----------plot red with only median of medians for all ROIs-------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
ii = (1:L)';
c = num2cell(jet(length(ii)));
for i = 1:L
   rH(i) = plot(fr{i},smooth(redMed{i},win2,'moving'),'color',cell2mat(c(i,:)),'LineWidth',6);
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
legend(cellstr(data_links.name));
savefile = [mergedpath,'RedChannel_median_of_medians_all_Treatments.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------%

