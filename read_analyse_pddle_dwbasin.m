function [pd_data] = read_analyse_pddle_dwbasin(data_path,plots_path,nheader,tmeas_s,fsamp_Hz)
%READ_ANALYSE_WG_DWBASIN 
% works similarly to read_analyse_wg_dwbasin.m but the plotting is much simpler
% it just compares the pd_data.angle of the pddata stucture between two repetition (Yf must be the same)
% it is storngly advised that data for this are placed in different folder that the WG data. Just pick the two runs you want to compare and put 
% them in a new folder for data path

%   Inputs: 
%   (1) data_path          e.g., 'E:/wave_gauge/' etc. [string]
%   (2) plots_path         e.g., 'E:/plots/' etc. [string]
%   (3) nheader            number of header lines [int OR int array]
%   (4) tmeas_s            duration of data (s) [int OR int array]
%   (5) fsamp_Hz           sampling frequency (Hz) [int OR int array]

% Andrew Smith - Konstantinos Chasapis 2023

%% Section 1: Call to data path, file ID, figure declarations

% 1.1 Call to data path
cd(data_path);

% 1.2 File ID
files = dir('surf_el__*.txt');

%% Section 2: Reading data
pd_data = [];

for j = 1:length(files)
    [t,pddata]=  fReadDWB_PDLDIAGS_AWS_KC(files(j).name,nheader(j),tmeas_s(j),fsamp_Hz(j));
    pd_data(j).t = t(:);
    pd_data(j).angle = pddata.angle;
 
    disp(['Data from file ' num2str(j) ' of ' num2str(length(files)) ' read in.']);
end

%% Section 3: Analysis
% .*180/pi rad to deg



%% Section 4: Plots
count = 1;
cd(plots_path)  
count = count + 1;
figure(count); 
fig_Width       = 45;
fig_AspectRatio = [3.5 1];
left_Margin     = 2.4;
bottom_Margin   = 1.5;
right_Margin    = 1.4;
top_Margin      = 0.5;
haxes           = tightPlots_Axes(1, 1,  fig_Width, ...
                  fig_AspectRatio,	...
                  [0 0], ...
                  [bottom_Margin top_Margin], ...
                  [left_Margin right_Margin], 'centimeters');
                                        
set(gcf, 'color', 'w')

subplot(1,3,1)
pdnum = 18;
 for j = 1:length(files)
     plot(pd_data(j).t(:,1),pd_data(j).angle(:,pdnum),'linewidth',2,'DisplayName',['wave paddle: ',num2str(pdnum)]); grid on; hold on;
     legend show
     legend('Interpreter','latex');
 end
% axes(haxes);
xlabel('time (s)'); ylabel('angle (rad)','interpreter', 'latex'); xlim([0 64]);
title('Wave paddle angle');

subplot(1,3,2)
pdnum = 28;
 for j = 1:length(files)
     plot(pd_data(j).t(:,1),pd_data(j).angle(:,pdnum) ,'linewidth',2,'DisplayName',['wave paddle: ',num2str(pdnum)]); grid on; hold on;
     legend show
     legend('Interpreter','latex');

 end
xlabel('time (s)'); ylabel('angle (rad)','interpreter', 'latex'); xlim([0 64]);
title('Wave paddle angle');

subplot(1,3,3)
pdnum = 38;
 for j = 1:length(files)
     plot(pd_data(j).t(:,1),pd_data(j).angle(:,pdnum),'linewidth',2,'DisplayName',['wave paddle: ',num2str(pdnum)]); grid on; hold on;
     legend show
     legend('Interpreter','latex');
 end
xlabel('time (s)'); ylabel('angle (rad)','interpreter', 'latex'); xlim([0 64]);
title('Wave paddle angle');


print('-dpng','-r600',['Figure_paddles_', files(1).name(10:end-4),'_Vs_',files(2).name(10:end-4)]);
savefig(['Figure_paddles_', files(1).name(10:end-4),'_Vs_',files(2).name(10:end-4),'.fig'])
close all
cd(data_path)






































