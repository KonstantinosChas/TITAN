function[micro_data] = read_analyse_micro(data_path,plots_path)
% micro_time,micro_surf,mrwaw,micro_spectra,micro_s_tot
% loads data from FReadMicro, data must be places in parent folder and names of specific folders (mi_datatime) should be specified in data_path
% check bblcam_dgnstics how to use the fucntion
% Input: path to data and path to save plots
% Output: micro data struct for each folder containing - names, time, surface elevation, spectra, total energy (variance)
% Konstantinos Chasapis 2023

%% Section 1: load data from fReadMICRO

legends = {'yf =  0    m T1','yf = -0.4 m T1','yf = -0.8 m T1','yf =  0    m T2','yf = -0.4 m T2','yf = -0.8 m T2'}';  % this is random based on my present cases
files = numel(data_path); % number of micro files
fieldname = cell(1,files); 
for i = 1:files
    fieldname{i} = char(data_path{1,i}(1,51:end));  % names of fields

    [t.(fieldname{i}),m.(fieldname{i}),mraw.(fieldname{i}),S_m.(fieldname{i}),dist_SWL.(fieldname{i})] = fReadMICRO_AWS_KC(data_path{i},'mi_001.txt',1,2,64,128,4);
end
for i = 1:files
    S_tot.(fieldname{i}) = sum(S_m.(fieldname{i})); % total energy per wave run
end
for i = files:-1:1
    micro_data(i).name = fieldname{i};
    micro_data(i).t = t.(fieldname{i});
    micro_data(i).surfel = m.(fieldname{i});
    micro_data(i).S = S_m.(fieldname{i});
    micro_data(i).Stot = S_tot.(fieldname{i});
    micro_data(i).raw = mraw.(fieldname{i});
end

%% Section 2: Analysis
% spike removal
% first derivative inspection 
for i = files: -1:1
    m_smooth(i).surfel_md = movmedian(micro_data(i).raw,12);
    m_smooth(i).deriv1 = diff(micro_data(i).raw);

    m_smooth(i).indices = find(m_smooth(i).deriv1<0.02 & m_smooth(i).deriv1>-0.02);
    m_smooth(i).surfel_fl = micro_data(i).raw(m_smooth(i).indices);
    m_smooth(i).t_fl = micro_data(i).t(m_smooth(i).indices);
    
    m_smooth(i).deriv11 = diff(m_smooth(i).surfel_fl);
    m_smooth(i).indices = find(m_smooth(i).deriv11<0.005 & m_smooth(i).deriv11>-0.005);
    m_smooth(i).surfel_fl11 = m_smooth(i).surfel_fl(m_smooth(i).indices);
    m_smooth(i).t_fl11 = m_smooth(i).t_fl(m_smooth(i).indices);
end



cd(plots_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------------Plots micro----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_comments = '_';
cc= [rand(1,numel(fieldname))',rand(1,numel(fieldname))',rand(1,numel(fieldname))']; % random colors for plots

% cc = jet(files);  % colorbar
A = figure();
sizefont = 32;
subplot(1,2,1);
set(A,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);   
set(A,'PaperSize',[90 90]);
set(A,'PaperUnits', 'centimeters');
set(A,'Color',[1 1 1]);
markers = {'-','-.',':','-','-.',':'};
fig_title = ['microsonic_smoothing',plot_comments];
title('surface elevation')
for i = 1:1
%    plot(micro_data(i).t,dist_SWL.(fieldname{i})-m_smooth(i).surfel_md,'linewidth',2,'color',cc(i+1,:,:),'LineStyle',markers{i+1}); hold on; % 

    plot(micro_data(i).t,dist_SWL.(fieldname{i})- micro_data(i).raw,'color',[0.9290 0.6940 0.1250],'linewidth',0.2,'LineStyle',markers{i+2}); hold on; % 
    plot(m_smooth(i).t_fl,dist_SWL.(fieldname{i}) - m_smooth(i).surfel_fl,'color',fClr(1),'linewidth',4,'LineStyle',markers{i}); hold on; %  
    plot(m_smooth(i).t_fl11,dist_SWL.(fieldname{i}) - m_smooth(i).surfel_fl11,'color',fClr(2),'linewidth',2,'LineStyle',markers{i}); hold on; %  
  
end
grid on; box on;
ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
set(gca,'Fontsize',sizefont);
xticks(0:10:64); 
% ylim([-0.3 0.4]);
yticks(-0.40:0.05:0.40);
lgd = legend('raw','smoothed once','smoothed twice'); 
fontsize(lgd,14,'points');

subplot(1,2,2);
title('first derivative');
for i = 1:1
    plot(micro_data(i).t(1:end-1),m_smooth(i).deriv1,'color',fClr(1),'linewidth',1.5,'LineStyle',markers{i}); hold on; %   
    plot(m_smooth(i).t_fl(1:end-1),m_smooth(i).deriv11,'color',fClr(2),'linewidth',1.5,'LineStyle',markers{i}); hold on; %   

end
grid on; box on;
ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
set(gca,'Fontsize',sizefont);
xticks(0:10:64); 
% ylim([-0.3 0.4]);
% yticks(-0.30:0.05:0.40);
lgd = legend('deriv of raw','deriv after 1st smooth'); 
fontsize(lgd,14,'points');
%--------------------------------------------------------------------------
cd(plots_path);
% exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
print(A,'-djpeg','-r500',[fig_title,'.jpg']);
saveas(A,[fig_title,'.fig']);
%--------------------------------------------------------------------------

%% -- good plots - uncomment and use
% cc = jet(files);  % colorbar
% A = figure();
% sizefont = 32;
% subplot(1,2,1);
% set(A,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);   
% set(A,'PaperSize',[90 90]);
% set(A,'PaperUnits', 'centimeters');
% set(A,'Color',[1 1 1]);
% markers = {'-','-.',':','-','-.',':'};
% fig_title = ['microsonic_yf_comp_',plot_comments];
% for i = 1:files/2
%     plot(t.(fieldname{i})(1:wind:end),m.(fieldname{i})(1:wind:end),'linewidth',1.5,'color',cc(i,:,:),'LineStyle',markers{i}); hold on; %   
% end
% grid on; box on;
% ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
% xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
% set(gca,'Fontsize',sizefont);
% xticks(0:10:64); 
% ylim([-0.3 0.4]);
% yticks(-0.30:0.05:0.40);
% lgd = legend(legends{1},legends{2},legends{3}); 
% fontsize(lgd,14,'points');
% 
% 
% disp('do you want to diplay two microsonic runs? y/n')
% % in = input('','s');
% in = 'y';
% if strcmp(in,'y')
% subplot(1,2,2);
% set(A,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);  
% set(A,'PaperSize',[90 90]);
% set(A,'PaperUnits', 'centimeters');
% set(A,'Color',[1 1 1]);
% for i = (files/2)+1:files 
%     plot(t.(fieldname{i})(1:wind:end),m.(fieldname{i})(1:wind:end),'linewidth',1.5,'color',cc(i,:,:),'LineStyle',markers{i}); hold on;
% end
% lgd = legend(legends{4},legends{5},legends{6}); 
% fontsize(lgd,14,'points');
% 
% grid on; box on;
% ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
% xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
% set(gca,'Fontsize',sizefont); 
% xticks(0:10:64); 
% ylim([-0.3 0.4]);
% yticks(-0.30:0.05:0.40); 
% else
%     disp('only one subplot this time');
% end
% 
% %--------------------------------------------------------------------------
% cd(plots_path);
% % exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
% print(A,'-djpeg','-r500',[fig_title,'.jpg']);
% saveas(A,[fig_title,'.fig']);
% %% --------------------------------------------------------------------------
% 
% % Surf_elev comparison between runs same yf
% cc = jet(files);
% A = figure();
% sizefont = 32;
% left = 0.25;  
% bottom = 0.23;
% width = 0.68; 
% height = 0.72;
% 
% subplot(1,2,1);
% set(A,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);   
% set(A,'PaperSize',[90 90]);
% set(A,'PaperUnits', 'centimeters');
% set(A,'Color',[1 1 1]);
% markers = {'-','-.',':','-','-.',':'};
% fig_title = ['microsonic_rep_comp',plot_comments];
% 
% plot(t.(fieldname{1})(1:wind:end),m.(fieldname{1})(1:wind:end),'linewidth',1.5,'color',cc(1,:,:),'LineStyle',markers{1}); hold on;
% plot(t.(fieldname{3})(1:wind:end),m.(fieldname{4})(1:wind:end),'linewidth',1.5,'color',cc(6,:,:),'LineStyle',markers{3});
% 
% grid on; box on;
% ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
% xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
% set(gca,'Fontsize',sizefont);
% xticks(0:10:64); 
% ylim([-0.3 0.4]);
% yticks(-0.30:0.05:0.40);
% lgd = legend(legends{1},legends{4}); fontsize(lgd,14,'points');
% %-------------------------
% subplot(1,2,2);
% set(A,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);  
% set(A,'PaperSize',[90 90]);
% set(A,'PaperUnits', 'centimeters');
% set(A,'Color',[1 1 1]);
% 
% plot(t.(fieldname{2})(1:wind:end),m.(fieldname{2})(1:wind:end),'linewidth',1.5,'color',cc(1,:,:),'LineStyle',markers{2}); hold on;
% plot(t.(fieldname{4})(1:wind:end),m.(fieldname{5})(1:wind:end),'linewidth',1.5,'color',cc(6,:,:),'LineStyle',markers{6});
% 
% grid on; box on;
% ylabel('$\eta$ (m)','Fontsize', sizefont,'interpreter', 'latex');
% xlabel('time (s)','Fontsize', sizefont,'interpreter', 'latex');
% set(gca,'Fontsize',sizefont); 
% xticks(0:10:64); 
% ylim([-0.3 0.4]);
% yticks(-0.30:0.05:0.40); 
% lgd = legend(legends{2},legends{5}); fontsize(lgd,14,'points');
% %--------------------------------------------------------------------------
% % exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
% print(A,'-djpeg','-r500',[fig_title,'.jpg']);
% saveas(A,[fig_title,'.fig']);
% %--------------------------------------------------------------------------
% 
% % spectra 
% A = figure();
% sizefont = 32;
% set(A,'units','normalized','outerposition',[0.1 0.1 0.5 0.8]);   
% set(A,'PaperSize',[90 90]);
% set(A,'PaperUnits', 'centimeters');
% set(A,'Color',[1 1 1]);
% fig_title = ['microsonic_spectra',plot_comments];
% for i = 1:files
%     loglog(S_m.(fieldname{i})(:,1),S_m.(fieldname{i})(:,2),'linewidth',1,'color',cc(i,:,:),'LineStyle',markers{i}); hold on;
% end
% grid on;
% set(gca,'XScale','log','YScale','log')
% set(gca,'Fontsize',sizefont); 
% ylabel('$S_n \mathrm{(m^2/Hz)} $','Fontsize', sizefont,'interpreter', 'latex');
% xlabel('Frequency (Hz)','Fontsize', sizefont,'interpreter', 'latex');
% 
% lgd = legend(legends{:}); fontsize(lgd,14,'points');
% %--------------------------------------------------------------------------
% % exportgraphics(P1,[fig_title,'.emf'],'ContentType','vector');  % other formats pdf, eps, tiff
% print(A,'-djpeg','-r500',[fig_title,'.jpg']);
% saveas(A,[fig_title,'.fig']);
% %--------------------------------------------------------------------------
close all;



end