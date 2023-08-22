function[hydro_data] = read_analyse_hydro(data_path,plots_path)
% Input - output: similar to read_analyse_micro.m

%% Section 1: Load data 
files = numel(data_path); % number of micro files
fieldname = cell(1,files); 
for i = files:-1:1
    cd(data_path{i});
    fieldname{i} = char(data_path{1,i}(1,51:end));  % names of fields
    data(i).a = dlmread('hy_001.txt');
    hydro_data(i).name = fieldname{i};
    hydro_data(i).time = data(i).a(1:1:end,1);
    hydro_data(i).p = data(i).a(1:1:end,2);
end


%% Section 2: Analysis
Sdb = -210; % sensitivity in dB (-146,-200 alternatives)
Fmin = 3; % low frequency threshold for high pass filter (wave filtering)
fs = 1/(hydro_data(1).time(3)-hydro_data(1).time(2));  % high frequency threshold for high pass filter (wave filtering)
window = [1 64];
nfft = 2^nextpow2(fs/Fmin);
for i = files:-1:1
    hydro_data(i).pfl = movmedian(hydro_data(i).p,30);   % ----SMOOTHING OPTIONS  1----
%     [hydro_data(i).pfl] = FUNCTION_Filtering_Components([Fmin fs],hydro_data(i).p,fs);  % 2 - 

    hydro_data(i).pwindow = hydro_data(i).p(hydro_data(i).time <= window(2) & hydro_data(i).time>= window(1));
end
for i = files:-1:1
    [~,hydro_data(i).F, hydro_data(i).tSpectrogram, hydro_data(i).SppSpectrogram] = spectrogram(hydro_data(i).pfl,hanning(fs/Fmin),nfft/2,nfft,fs,'psd');
    hydro_data(i).SppSpectrogram = (-Sdb + 10 * log10(hydro_data(i).SppSpectrogram));  %If working with the signal in volts
end


%% Section 3: Plots
cd(plots_path);
% A = figure ();
% fig_Width       = 45;
% fig_AspectRatio = [3.5 1];
% left_Margin     = 2.4;
% bottom_Margin   = 1.1;
% right_Margin    = 1.4;
% top_Margin      = 0.5;
% haxes           = tightPlots_Axes(1, 1,  fig_Width, ...
%                   fig_AspectRatio,	...
%                   [0 0], ...
%                   [bottom_Margin top_Margin], ...
%                   [left_Margin right_Margin], 'centimeters');
%                                         
% set(gcf, 'color', 'w')
% 
% fig_title = 'pressure time series';
% % set(A,'units','normalized','outerposition',[0.1 0. 0.7 0.7]);
% % set(A,'PaperSize',[70 70]);
% % set(A,'PaperUnits', 'centimeters');
% % set(A,'Color',[1 1 1]);
% for i = 1:1
%     p(1) = plot(hydro_data(i).time,hydro_data(i).p*10^-9,'color',fClr(1));hold on
%     p(2) = plot(hydro_data(i).time,hydro_data(i).pfl*10^-9,'color',fClr(2),'linewidth',0.3); % Multiply by10^-9 to convert to kilopascals
% end
% legend(p,'Original Signal',sprintf('High-pass filter'),...
%     'interpreter','latex','Location','southoutside','FontSize',22);
% xlim([window(1) window(end)]);
% xlabel('Time [s]','interpreter','latex','fontsize',22);
% ylabel('Pressure [$kPa$]','interpreter','latex','fontsize',22);
% set(gca,'fontsize',22);
% grid on; box on;
% %--------------------------------------------------------------------------
% cd(plots_path);
% % exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
% print(A,'-djpeg','-r500',[fig_title,'.jpg']);
% saveas(A,[fig_title,'.fig']);
% %--------------------------------------------------------------------------

% spectrogram
A = figure ();
fig_title ='spectrogram_f';
set(A,'units','normalized','outerposition',[0.1 0. 0.45 0.8]);
set(A,'PaperSize',[70 70]);
set(A,'PaperUnits', 'centimeters');
set(A,'Color',[1 1 1]); clear p

multi = 1000; % this is for KiloHertz
for i = 1:1
    pcolor(window(1) + hydro_data(i).tSpectrogram,hydro_data(i).F/multi,hydro_data(i).SppSpectrogram); 
end
shading flat;
% lighting phong;
% lightangle(-40,30); % -90,30
% camlight('left');
% shading interp;
% EdgeColor = 'none';
c = colorbar('Location','southoutside'); h = colormap(jet);
grid on;
grid minor;
set(gca, 'YScale', 'log')
ylabel(c,'$E_{var}$ $[dB$ $re.$ $\mu Pa^2/Hz]$','Interpreter','latex','fontsize',18);
% caxis([0 max(max(SppSpectrogram))]);
ylabel('Frequency [$kHz$]','interpreter','latex','fontsize',18)
yticks([0.05 0.1 0.5 1 5 10 50 100 250]);
ylim([Fmin*(1/multi) (fs/2)*(1/multi)]);
xlabel('Time [$s$]','interpreter','latex','fontsize',18);
% xlim([window(1) window(end)]);
% xlim([31 33]);
set(gca,'fontsize',18);
grid on; box on;
%--------------------------------------------------------------------------
cd(plots_path);
% exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
print(A,'-djpeg','-r500',[fig_title,'.jpg']);
saveas(A,[fig_title,'.fig']);
%--------------------------------------------------------------------------
    
    
    
    
    
    



end








 