function [t,eta,depth] = fReadWG_AWS_KC(file,ngauge,nheader,tEnd,sampleRate,tSWL)
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

% AWS, Jul 2023

% read in data
fid = fopen(file,'r');
form = ['%s',repmat('%*f',[1,336]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f',repmat('%*f',[1,3]),'%f','%*[^\n]'];
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

% construct eta & t
eta = zeros(nt,ngauge);
for i = 1:ngauge
    eta(:,i) = txtData{i+1}(1:nt);
end

depth = zeros(1,ngauge);

% if depth data is at the begining of the file
if nargin > 5
    fprintf('(d) \n \n')
    for i = 1:ngauge
%         swl = eta(t<tSWL,i);
%         depth(i) = mean(eta(swl(~isnan(swl)),i));
        depth(i) = mean(eta(t<tSWL,i));
        % steven added this to incoorporate cases where depth data is
        % obtained otherwise
        if isnan(depth(i))
            depth(i) = 0;
            fprintf(1,'no depth value is interpreted, please incoorporate this another way.')
        end
        eta(:,i) = eta(:,i) - depth(i);
    end
end

