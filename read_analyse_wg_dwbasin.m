function [wg_data] = read_analyse_wg_dwbasin(data_path,plots_path,nheader,ngauges,tmeas_s,fsamp_Hz,pw_avg_s)
%READ_ANALYSE_WG_DWBASIN Reads and analyses wave gauge data, plots, and
%saves figures

%   Inputs: 
%   (1) data_path          e.g., 'E:/wave_gauge/' etc. [string]
%   (2) plots_path         e.g., 'E:/plots/' etc. [string]
%   (3) nheader            number of header lines [int OR int array]
%   (4) ngauges            number of wave gauges [int OR int array]
%   (5) tmeas_s            duration of data (s) [int OR int array]
%   (6) fsamp_Hz           sampling frequency (Hz) [int OR int array]
%   (7) pw_avg_s           pre-wave time for avg removal (s) [dbl OR array]
%   (8) gauge_id           which gauge for plots [int OR int arrayt] ------> optional see comments in script

%% Section 1: Call to data path, file ID, figure declarations

% 1.1 Call to data path
cd(data_path)

% 1.2 File ID
files = dir('surf_el__*.txt');

% Store some default declarations here for reference

% Run the below...
% data_path = '/media/aws/ICL_2TB/ICL_Work/Data/AWS_KC_2023_data/wg_data/run_4_1_04/surf_elev/'; cd(data_path);
% plots_path = '/media/aws/ICL_2TB/ICL_Work/Data/AWS_KC_2023_data/plots/';
% files = dir('surf_el__*.txt');
% ngauges = 16.*ones(length(files),1);
% nheader = 5.*ones(length(files),1);
% tmeas_s = 64.0.*ones(length(files),1);
% fsamp_Hz = 128.0.*ones(length(files),1);
% pw_avg_s = 4.0.*ones(length(files),1);
% gauge_id = 3.*ones(length(files),1);
%[~] = read_analyse_wg_dwbasin(data_path,plots_path,nheader,ngauges,tmeas_s,fsamp_Hz,pw_avg_s,gauge_id);

%% Section 2: Reading data
wg_data = [];

for j = 1:length(files)
    [t,eta,depth] = fReadWG_AWS_KC(files(j).name,ngauges(j),nheader(j),tmeas_s(j),fsamp_Hz(j),pw_avg_s(j));
    wg_data(j).t = t(:);
    wg_data(j).eta = eta(:,:);
    wg_data(j).depth = depth(:);
    clear t eta depth
    disp(['Data from file ' num2str(j) ' of ' num2str(length(files)) ' read in.']);
end

%% Section 3: Analysis of data

for j = 1:length(files)
        for k = 1:ngauges(j)
        % 3.1 eta^2
        wg_data(j).eta_sq(:,k) = wg_data(j).eta(:,k).^2;
        % 3.2 rms(eta)
        wg_data(j).eta_rms(k) = rms(wg_data(j).eta(:,k));
        % 3.3 rms(eta) using a moving average of tavg [s]
        tavg = 1/8; 
        wg_data(j).eta_mrms(:,k) = sqrt(movmean(wg_data(j).eta_sq(:,k),fsamp_Hz(j)*tavg));
        % 3.4 wave spectrum
        spec_out = spectf(wg_data(j).eta(:,k),1.0/fsamp_Hz(j),8);
        wg_data(j).f(:,k) = spec_out(:,1);                                                        % [Hz]
        wg_data(j).S_eta(:,k) = spec_out(:,2);                                                    % [m^2/Hz]
        clear spec_out
        % 3.5 integrated wave spectrum
        wg_data(j).S_eta_integ(k) = nansum(nansum(wg_data(j).f(:,k).*wg_data(j).S_eta(:,k)));           % [m^2]    
        end
        cd(data_path)
        run_string = data_path(end-20:end-11);
        ts_string = files(j).name(10:end-4);
        save(['wg_data_' ts_string '.mat'],'wg_data');    
        disp(['File : wg_data_' ts_string '.mat']);  
end

%% Section 4: Plotting of data

repeatIndxs = [1:25:length(files)];
numRepeats = length(repeatIndxs);
count = 0;
% n =  1: numRepeats not used anymore in this version 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure Set 1: Eta, Spectra, Along-Crest Profile for WG 1-8  -------------------------------------------(plot every 25th yf, e.g. 1, 26, etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(files)
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
    
    h1 = plot(wg_data(j).t,wg_data(j).eta(:,1),'linewidth',2); grid on; hold on;
    h2 = plot(wg_data(j).t,wg_data(j).eta(:,2),'linewidth',2); grid on; hold on;
    h3 = plot(wg_data(j).t,wg_data(j).eta(:,3),'linewidth',2); grid on; hold on;
    h4 = plot(wg_data(j).t,wg_data(j).eta(:,4),'linewidth',2); grid on; hold on;
    h5 = plot(wg_data(j).t,wg_data(j).eta(:,5),'linewidth',2); grid on; hold on;
    h6 = plot(wg_data(j).t,wg_data(j).eta(:,6),'linewidth',2); grid on; hold on;
    h7 = plot(wg_data(j).t,wg_data(j).eta(:,7),'linewidth',2); grid on; hold on;
    h8 = plot(wg_data(j).t,wg_data(j).eta(:,8),'linewidth',2); grid on; hold on;
    
    xlabel('Time (s)'); ylabel('$\eta$ (m)','interpreter', 'latex')
    title('Timeseries: WG 1-8 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8],{'WG 1','WG 2','WG 3','WG 4','WG 5','WG 6','WG 7','WG 8'},'Location','NorthEast','Interpreter','latex'); 
    ylim([-0.2 0.3]);
    
    subplot(1,3,2)
    
    h1 = loglog(wg_data(j).f(:,1),wg_data(j).S_eta(:,1),'linewidth',2); grid on; hold on;
    h2 = loglog(wg_data(j).f(:,2),wg_data(j).S_eta(:,2),'linewidth',2); grid on; hold on;
    h3 = loglog(wg_data(j).f(:,3),wg_data(j).S_eta(:,3),'linewidth',2); grid on; hold on;
    h4 = loglog(wg_data(j).f(:,4),wg_data(j).S_eta(:,4),'linewidth',2); grid on; hold on;
    h5 = loglog(wg_data(j).f(:,5),wg_data(j).S_eta(:,5),'linewidth',2); grid on; hold on;
    h6 = loglog(wg_data(j).f(:,6),wg_data(j).S_eta(:,6),'linewidth',2); grid on; hold on;
    h7 = loglog(wg_data(j).f(:,7),wg_data(j).S_eta(:,7),'linewidth',2); grid on; hold on;
    h8 = loglog(wg_data(j).f(:,8),wg_data(j).S_eta(:,8),'linewidth',2); grid on; hold on;
    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 1-8 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8],{'WG 1','WG 2','WG 3','WG 4','WG 5','WG 6','WG 7','WG 8'},'Location','NorthEast','Interpreter','latex'); 
    ylim([1e-12 1e2]);
    
    subplot(1,3,3)
    gaugePos = [-0.4:0.2:1.0];
    xx = gaugePos;
    yy = wg_data(j).S_eta_integ(1:8);
    
    h1 = scatter(xx(1),yy(1),100,'o','linewidth',2); grid on; hold on; 
    h2 = scatter(xx(2),yy(2),100,'o','linewidth',2); grid on; hold on; 
    h3 = scatter(xx(3),yy(3),100,'o','linewidth',2); grid on; hold on; 
    h4 = scatter(xx(4),yy(4),100,'o','linewidth',2); grid on; hold on; 
    h5 = scatter(xx(5),yy(5),100,'o','linewidth',2); grid on; hold on;  
    h6 = scatter(xx(6),yy(6),100,'o','linewidth',2); grid on; hold on; 
    h7 = scatter(xx(7),yy(7),100,'o','linewidth',2); grid on; hold on; 
    h8 = scatter(xx(8),yy(8),100,'o','linewidth',2); grid on; hold on; 
    
    p = polyfit(xx(2:6),yy(2:6),2);
    xxf = linspace(-3,3,100);
    yyf = polyval(p,xxf); 
    
    h9 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Position (m)'); ylabel('$\int S_{\eta}$ (m$^{2}$)','interpreter', 'latex');
    title('Symmetry: WG 1-8 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],{num2str(round(yy(1),4)),num2str(round(yy(2),4)),num2str(round(yy(3),4)),num2str(round(yy(4),4)),num2str(round(yy(5),4)),num2str(round(yy(6),4)),num2str(round(yy(7),4)),num2str(round(yy(8),4)),'Fit'},'Location','NorthEast','Interpreter','latex'); 
    xlim([-2 4]);
    set(gca,'XTick',[-2:2:2]);
    ylim([0 0.03]);
    
    % Save and clear figures
    run_string = data_path(end-20:end-11);
    ts_string = files(j).name(10:end-4);
    
    figure(count);
    cd(plots_path)    
    print('-dpng','-r600',['Figure_' num2str(count) '_' ts_string])
    savefig(['Figure_' num2str(count) '_' ts_string '.fig'])
    close all
end
cd(data_path)

%% Figure Set 2: Eta, Spectra, Along-Crest Profile for WG 9-16 ------------------------------------------------(e.g. plot every 25th yf, e.g. 1, 26, etc.)

for j = 1:length(files)
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
    
    h1 = plot(wg_data(j).t,wg_data(j).eta(:,9),'linewidth',2); grid on; hold on;
    h2 = plot(wg_data(j).t,wg_data(j).eta(:,10),'linewidth',2); grid on; hold on;
    h3 = plot(wg_data(j).t,wg_data(j).eta(:,11),'linewidth',2); grid on; hold on;
    h4 = plot(wg_data(j).t,wg_data(j).eta(:,12),'linewidth',2); grid on; hold on;
    h5 = plot(wg_data(j).t,wg_data(j).eta(:,13),'linewidth',2); grid on; hold on;
    h6 = plot(wg_data(j).t,wg_data(j).eta(:,14),'linewidth',2); grid on; hold on;
    h7 = plot(wg_data(j).t,wg_data(j).eta(:,15),'linewidth',2); grid on; hold on;
    h8 = plot(wg_data(j).t,wg_data(j).eta(:,16),'linewidth',2); grid on; hold on;
    
    xlabel('Time (s)'); ylabel('$\eta$ (m)','Interpreter','latex')
    title('Timeseries: WG 9-16 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8],{'WG 9','WG 10','WG 11','WG 12','WG 13','WG 14','WG 15','WG 16'},'Location','NorthEast','Interpreter','latex'); 
    ylim([-0.2 0.3]);
    
    subplot(1,3,2)
    
    h1 = loglog(wg_data(j).f(:,9),wg_data(j).S_eta(:,9),'linewidth',2); grid on; hold on;
    h2 = loglog(wg_data(j).f(:,10),wg_data(j).S_eta(:,10),'linewidth',2); grid on; hold on;
    h3 = loglog(wg_data(j).f(:,11),wg_data(j).S_eta(:,11),'linewidth',2); grid on; hold on;
    h4 = loglog(wg_data(j).f(:,12),wg_data(j).S_eta(:,12),'linewidth',2); grid on; hold on;
    h5 = loglog(wg_data(j).f(:,13),wg_data(j).S_eta(:,13),'linewidth',2); grid on; hold on;
    h6 = loglog(wg_data(j).f(:,14),wg_data(j).S_eta(:,14),'linewidth',2); grid on; hold on;
    h7 = loglog(wg_data(j).f(:,15),wg_data(j).S_eta(:,15),'linewidth',2); grid on; hold on;
    h8 = loglog(wg_data(j).f(:,16),wg_data(j).S_eta(:,16),'linewidth',2); grid on; hold on;
    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 9-16 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8],{'WG 9','WG 10','WG 11','WG 12','WG 13','WG 14','WG 15','WG 16'},'Location','NorthEast','Interpreter','latex'); 
    ylim([1e-12 1e2]);
    
    subplot(1,3,3)
    gaugePos = [-0.4:0.2:1.0];
    xx = gaugePos;
    yy = wg_data(j).S_eta_integ(9:16);
    
    h1 = scatter(xx(1),yy(1),100,'o','linewidth',2); grid on; hold on; 
    h2 = scatter(xx(2),yy(2),100,'o','linewidth',2); grid on; hold on; 
    h3 = scatter(xx(3),yy(3),100,'o','linewidth',2); grid on; hold on; 
    h4 = scatter(xx(4),yy(4),100,'o','linewidth',2); grid on; hold on; 
    h5 = scatter(xx(5),yy(5),100,'o','linewidth',2); grid on; hold on;  
    h6 = scatter(xx(6),yy(6),100,'o','linewidth',2); grid on; hold on; 
    h7 = scatter(xx(7),yy(7),100,'o','linewidth',2); grid on; hold on; 
    h8 = scatter(xx(8),yy(8),100,'o','linewidth',2); grid on; hold on; 
    
    p = polyfit(xx(2:6),yy(2:6),2);
    xxf = linspace(-3,3,100);
    yyf = polyval(p,xxf); 
    
    h9 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Position (m)'); ylabel('$\int S_{\eta}$ (m$^{2}$)','interpreter', 'latex');
    title('Symmetry: WG 9-16 Series');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],{num2str(round(yy(1),4)),num2str(round(yy(2),4)),num2str(round(yy(3),4)),num2str(round(yy(4),4)),num2str(round(yy(5),4)),num2str(round(yy(6),4)),num2str(round(yy(7),4)),num2str(round(yy(8),4)),'Fit'},'Location','NorthEast','Interpreter','latex'); 
    xlim([-2 4]);
    set(gca,'XTick',[-2:2:2]);
    ylim([0 0.03]);
    
    % Save and clear figures
    run_string = data_path(end-20:end-11);
    ts_string = files(j).name(10:end-4);
    
    figure(count);
    cd(plots_path)    
    print('-dpng','-r600',['Figure_' num2str(count) '_' ts_string])
    savefig(['Figure_' num2str(count) '_' ts_string '.fig'])
    close all
end
cd(data_path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure Set 3: Change in Total Energy with Wave and Along-Crest for a single row (ONE ROW - ONE RUN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1: length(files)

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
    
    gaugePos = [0:0.2:1.0];
    xx = gaugePos;
    
    yy = -1.*(wg_data(j).S_eta_integ(11:16) - wg_data(j).S_eta_integ(3:8));
    
    h1 = scatter(xx(1),yy(1),100,'o','linewidth',2); grid on; hold on; 
    h2 = scatter(xx(2),yy(2),100,'o','linewidth',2); grid on; hold on; 
    h3 = scatter(xx(3),yy(3),100,'o','linewidth',2); grid on; hold on;  
    h4 = scatter(xx(4),yy(4),100,'o','linewidth',2); grid on; hold on; 
    h5 = scatter(xx(5),yy(5),100,'o','linewidth',2); grid on; hold on; 
    h6 = scatter(xx(6),yy(6),100,'o','linewidth',2); grid on; hold on; 
    
    p = polyfit(xx,yy,1);
    xxf = linspace(0,1,100);
    yyf = polyval(p,xxf); 
    
    h7 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Gauge Position (m)'); ylabel('$\Delta \int S_{\eta}$ (m$^{2}$)','interpreter', 'latex')
    title('Along-X (Along Wave Dir)');
    l1 = legend([h1,h2,h3,h4,h5,h6,h7],{'WG 11 - WG 3','WG 12 - WG 4','WG 13 - WG 5','WG 14 - WG 6','WG 15 - WG 7','WG 16 - WG 8','Fit'},'Location','NorthEast'); 
    xlim([-0.2 1.2]);
    set(gca,'XTick',[0:0.2:1]);
%     ylim([0 0.005]);
            
    subplot(1,3,2)
    
    initialGaugePos = [0:0.2:0.8];
    xx = []; 
    xx = initialGaugePos;
    yy = [];
    yy(1) = -1.*(wg_data(j).S_eta_integ(4) - wg_data(j).S_eta_integ(3));
    yy(2) = -1.*(wg_data(j).S_eta_integ(5) - wg_data(j).S_eta_integ(4));
    yy(3) = -1.*(wg_data(j).S_eta_integ(6) - wg_data(j).S_eta_integ(5));
    yy(4) = -1.*(wg_data(j).S_eta_integ(7) - wg_data(j).S_eta_integ(6));
    yy(5) = -1.*(wg_data(j).S_eta_integ(8) - wg_data(j).S_eta_integ(7));
    
    h1 = scatter(xx(1),yy(1),100,'o','linewidth',2); grid on; hold on; 
    h2 = scatter(xx(2),yy(2),100,'o','linewidth',2); grid on; hold on; 
    h3 = scatter(xx(3),yy(3),100,'o','linewidth',2); grid on; hold on;  
    h4 = scatter(xx(4),yy(4),100,'o','linewidth',2); grid on; hold on; 
    h5 = scatter(xx(5),yy(5),100,'o','linewidth',2); grid on; hold on; 
    
    p = polyfit(xx,yy,1);
    xxf = linspace(0,1,100);
    yyf = polyval(p,xxf); 
    
    h6 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Initial Gauge Position (m)'); ylabel('$\Delta \int S_{\eta}$ (m$^{2}$)','interpreter', 'latex')
    title('Along-Y (Along Crest Dir)');
    l1 = legend([h1,h2,h3,h4,h5,h6],{'WG 4 - WG 3','WG 5 - WG 4','WG 6 - WG 5','WG 7 - WG 6','WG 8 - WG 7','Fit'},'Location','NorthEast'); 
%     ylim([0 0.0005]);
    xlim([-0.2 1.2]);
    set(gca,'XTick',[0:0.2:1]);
    
    subplot(1,3,3)
    
    h1 = loglog(wg_data(j).f(:,11),(wg_data(j).S_eta(:,11) ./ wg_data(1).S_eta(:,3)),'linewidth',2); grid on; hold on;
    h2 = loglog(wg_data(j).f(:,8),(wg_data(j).S_eta(:,8) ./ wg_data(1).S_eta(:,3)),'linewidth',2); grid on; hold on;
    h3 = loglog(wg_data(j).f(:,16),(wg_data(j).S_eta(:,16) ./ wg_data(1).S_eta(:,11)),'linewidth',2); grid on; hold on;
    
    xlabel('Frequency (Hz)'); ylabel('$\Delta S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Energy Loss by Frequency');
    l1 = legend([h1,h2,h3],{'WG 11/WG 3','WG 8/WG 3','WG 16/WG 11'},'Location','NorthEast'); 
    ylim([1e-6 1e6]);
    
    % Save and clear figures
    run_string = data_path(end-20:end-11);
    ts_string = files(j).name(10:end-4);
    
    figure(count);
    cd(plots_path)    
    print('-dpng','-r600',['Figure_' num2str(count) '_' ts_string])
    savefig(['Figure_' num2str(count) '_' ts_string '.fig'])
    close all
    cd(data_path)
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Many different runs combined - Figure Set 3: Eta, Spectra as y-focusing Position Changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = count + 1;
j = count;
figure(j); 
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
legends = {'WG 3 @ $y_{f}$ = -0.2','WG 4 @ $y_{f}$ = 0.0','WG 5 @ $y_{f}$ = +0.2','WG 6 @ $y_{f}$ = +0.4','WG 7 @ $y_{f}$ = +0.6','WG 8 @ $y_{f}$ = +0.8'};

if length(files) > 2
    subplot(1,2,1)
     for j = 1:length(files)
        h1(j) = plot(wg_data(j).t(:,1),wg_data(j).eta(:,j+2),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
     end
     legend show
     legend('Interpreter','latex');
    
    xlabel('Time (s)'); ylabel('$\eta$ (m)','Interpreter','latex')
    title('Timeseries: WG 1-8 Series');
    ylim([-0.2 0.3]);
            
    subplot(1,2,2)
    
     for j = 1:length(files)
        h1(j) = loglog(wg_data(j).f(:,1),wg_data(j).S_eta(:,j+2),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
     end
     legend show
     legend('Interpreter','latex');
    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 1-8 Series');
    ylim([1e-12 1e2]);

else
    subplot(1,3,1)
     for j = 1:length(files)
        h1(j) = plot(wg_data(j).t(:,1),wg_data(j).eta(:,j+2),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
     end
     legend show
     legend('Interpreter','latex');
    
    xlabel('Time (s)'); ylabel('$\eta$ (m)','Interpreter','latex')
    title('Timeseries: WG 1-8 Series');
    ylim([-0.2 0.3]);
            
    subplot(1,3,2)
    
     for j = 1:length(files)
        h1(j) = loglog(wg_data(j).f(:,1),wg_data(j).S_eta(:,j+2),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
     end
     legend show
     legend('Interpreter','latex');
    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 1-8 Series');
    ylim([1e-12 1e2]);

end

% energy dissipation between WG with different yf 
if length(files) > 2
    subplot(1,3,3)
    initialGauge = [3:1:7];
    xx = initialGauge;
    yy = [];

    m = 1; 
    for j = 1 : length(files) - 1
        yy(m) = wg_data(j+1).S_eta_integ(j+3) - wg_data(j).S_eta_integ(j+2);
        m = m + 1;
    end
%     clear legends
%     legends = zeros(1,length(yy));
%     legends(1) = num2str(round(yy(1),4));
%     for j = 1 : length(yy)
%         legends(j) = num2str(round(yy(j),4));
%     end

    clear m ;
    for j = 1 : length(yy)
        h1 = scatter(xx(j),yy(j),100,'o','linewidth',2,'DisplayName',char(num2str(round(yy(j),4)))); grid on; hold on;  
    end
    legend show
    legend('Interpreter','latex');
    xxf = linspace(1,8,100);
    yyf = 0.*ones(length(xxf),1);
    
    h6 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Initial Gauge'); ylabel('$\Delta \int S_{\eta}$ (m$^{2}$)','interpreter', 'latex');
    title('Total Energy Differences: WG 1-8 Series');
    set(gca,'XTick',[1:1:8]);
    ylim([-0.1 0.1]);
    set(gca,'YTick',[-0.1:0.05:0.1]);
    xlim([0 9]);
end
% Save and clear figures
run_string = data_path(end-20:end-11);
ts_string = files(j).name(10:end-4);

figure(count);
cd(plots_path)    
print('-dpng','-r600',['Figure_' num2str(count) '_' ts_string])
savefig(['Figure_' num2str(count) '_' ts_string '.fig'])
close all
cd(data_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same but for WG 8-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

legends = {'WG 11 @ $y_{f}$ = -0.2','WG 12 @ $y_{f}$ = 0.0','WG 13 @ $y_{f}$ = +0.2','WG 14 @ $y_{f}$ = +0.4','WG 15 @ $y_{f}$ = +0.6','WG 16 @ $y_{f}$ = +0.8'};

if length(files) > 2
    subplot(1,2,1)
    for j = 1 : length(files)
        h1(j) = plot(wg_data(j).t,wg_data(j).eta(:,j+10),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
    end
    legend show
    legend('Interpreter','latex');
    xlabel('Time (s)'); ylabel('$\eta$ (m)')
    title('Timeseries: WG 9-16 Series');
    
    ylim([-0.2 0.3]);
            
    subplot(1,2,2)
    for j = 1 : length(files)
        h1 = loglog(wg_data(j).f(:,11),wg_data(j).S_eta(:,j+10),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
    end
    legend show
    legend('Interpreter','latex');    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 9-16 Series');
    ylim([1e-12 1e2]);


else
    subplot(1,3,1)
    for j = 1 : length(files)
        h1(j) = plot(wg_data(j).t,wg_data(j).eta(:,j+10),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
    end
    legend show
    legend('Interpreter','latex'); 
    xlabel('Time (s)'); ylabel('$\eta$ (m)')
    title('Timeseries: WG 9-16 Series');
    
    ylim([-0.2 0.3]);
            
    subplot(1,3,2)
    for j = 1 : length(files)
        h1 = loglog(wg_data(j).f(:,11),wg_data(j).S_eta(:,j+10),'linewidth',2,'DisplayName',char(legends(j))); grid on; hold on;
    end
    legend show
    legend('Interpreter','latex');  
    
    xlabel('Frequency (Hz)'); ylabel('$S_{\eta}$ (m$^{2}$ Hz$^{-1}$)','interpreter', 'latex');
    title('Spectra: WG 9-16 Series');
    ylim([1e-12 1e2]);
end


if length(files) > 2
    subplot(1,3,3)
    initialGauge = [11:1:15];
    xx = initialGauge;
    yy = [];
    m = 1;
    for j = length(files) - 1
        yy(m) = wg_data(j+1).S_eta_integ(j+11) - wg_data(j).S_eta_integ(j+10);
        m = m + 1;
    end
    clear m;
    legends = num2str(round(yy(1),4));
    for j = 2 : length(yy)
        legends = [legends(j-1),num2str(round(yy(m),4))];
    end
    for j = length(yy)
        h1 = scatter(xx(j),yy(j),100,'o','linewidth',2,'DisplayName',char(legends(j))); grid on; hold on; 
    end
    xxf = linspace(9,16,100);
    yyf = 0.*ones(length(xxf),1);
    
    h6 = plot(xxf,yyf,'k','linewidth',2);
    
    xlabel('Initial Gauge'); ylabel('$\Delta \int S_{\eta}$ (m$^{2}$)','interpreter', 'latex');
    title('Total Energy Differences: WG 9-16 Series');
    legend show
    legend('Interpreter','latex'); 
    set(gca,'XTick',[9:1:16]);
    ylim([-0.1 0.1]);
    set(gca,'YTick',[-0.1:0.05:0.1]);
    xlim([8 17]);
end
% Save and clear figures
run_string = data_path(end-20:end-11);
ts_string = files(j).name(10:end-4);

figure(count);
cd(plots_path)    
print('-dpng','-r600',['Figure_' num2str(count) '_' ts_string])
savefig(['Figure_' num2str(count) '_' ts_string '.fig'])
close all
cd(data_path)








end

