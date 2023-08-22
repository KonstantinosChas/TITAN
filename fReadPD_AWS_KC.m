function [t,dwb_pdl_diag] = fReadDWB_PDLDIAGS_AWS_KC(file,nheader,tEnd,sampleRate)
% Reads surface element data produced by Edinburgh Design Wave Gauge
% Software, new version
% [t,eta,depth] = fReadWG_AWS_KC(file,ngauge,nheader,tEnd,sampleRate,tSWL)
% inputs:
%       file    - filename
%       ngauge  - number of gauges
%       nheader - number of headerlines
%       tEnd    - [sec], duration of the sampling
%       sampleRate - [Hz], this is used to check the file length in
%                    order to spot potential sampling problems
%       tSWL - [sec], optional, specifies the number of seconds at the
%              beginning of the record to be used as the averaged still
%              water level measurement
% outputs:
%       t       - time vector
%       eta     - surface elevation
%       depth   - averaged still water level measurement, over tSWL

% Structure array output units
% Angle [rad]
% Drive torque [Nm]
% Dynamic torque [Nm]
% Drive current [A]
% Demand [Nm]
% Offset angle [rad]

% AWS, Jul 2023

% read in data
fid = fopen(file,'r');
form = ['%s',repmat('%f',[1,336]),'%*[^\n]'];
txtData = textscan(fid,form,'headerlines',nheader,'delimiter','\t','treatAsEmpty','NULL','EmptyValue',0);
fclose(fid);
% report date and time when the txt file is written
fid = fopen(file,'r');
firstLine = textscan(fid,'%s %s %s %s %s %s',1,'delimiter','\t'); 
date = firstLine{3}{1};
time = firstLine{5}{1};
expname = firstLine{6}{1};
filename = textscan(file,'%s','delimiter','/'); filename = filename{1}{end};
fprintf(1,'\n EXPERIMENT: %s (%s) \n PERFORMED ON: %s \t AT: %s \n',filename,expname,date,time);

fclose(fid);

% check file
nt = length(1/sampleRate:1/sampleRate:tEnd);
if nt ~= length(txtData{2})
    warning(['Data is not of expected length. Expected nt=',num2str(nt),', actual nt=',num2str(length(txtData{1}))])
    ifTrunc = input('Do you want to truncate the end of the file? Y/n?','s');
    if strcmpi(ifTrunc,'n')
        nt = length(txtData{2});
    elseif strcmp(ifTrunc, 'y')
        fprintf('recorded truncated to expected length! \n');
    end
        
    t = linspace(1/sampleRate,nt/sampleRate,nt)';
else
    t = linspace(1/sampleRate,tEnd,nt)';
end

% create output structure array
dwb_pdl_diag = [];
numericData = zeros(length(t),length(txtData)-1);

for j = 2:length(txtData)
    numericData(:,j) = txtData{j}(1:end);
end
numericData(:,1) = [];

dwb_pdl_diag.t = t;
dwb_pdl_diag.timestamp = txtData{1}(1:end);
dwb_pdl_diag.angle = numericData(:,1:6:end);
dwb_pdl_diag.drive_torque = numericData(:,2:6:end);
dwb_pdl_diag.dynamic_torque = numericData(:,3:6:end);
dwb_pdl_diag.drive_current = numericData(:,4:6:end);
dwb_pdl_diag.demand = numericData(:,5:6:end);
dwb_pdl_diag.offset_angle = numericData(:,6:6:end);

ts_string = file(10:end-4);
save(['dwb_pdl_diag_' ts_string '.mat'],'dwb_pdl_diag');    
disp(['File : dwb_pdl_diag_' ts_string '.mat']);  
end

