function [t,m,mraw,S_m,dist_SWL] = fReadMICRO_AWS_KC(path,file,ngauge,nheader,tEnd,sampleRate,tSWL)
% Reads surface element data produced by microsonic
% Software, new version
% [t,eta,depth] = fReadWG_AWS_KC(file,ngauge,nheader,tEnd,sampleRate,tSWL)
% inputs:
%       file    - filename
%       ngauge  - number of microsonics
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
cd(path);
fid = fopen(file,'r');
form = ['%f','%f'];
txtData = textscan(fid,form,'headerlines',nheader,'delimiter','\t','treatAsEmpty','NULL','EmptyValue',0);
fclose(fid);

% check file
nt = length(1/sampleRate:1/sampleRate:tEnd);
if nt ~= length(txtData{2})
%     warning(['Data is not of expected length. Expected nt=',num2str(nt),', actual nt=',num2str(length(txtData{1}))])
%     ifTrunc = input('Do you want to truncate the end of the file? Y/n?','s');
%     if strcmpi(ifTrunc,'n')
        nt = length(txtData{2});
%     elseif strcmp(ifTrunc, 'y')
%         fprintf('recorded truncated to expected length! \n');
%     end
        
    t = linspace(1/sampleRate,nt/sampleRate,nt)';
else
    t = linspace(1/sampleRate,tEnd,nt)';
end

% construct V (Volts) & t
V = zeros(nt,ngauge);
for i = 1:ngauge
    V(:,i) = txtData{i+1}(1:nt);
end

% construct m from V (conversion)
mraw = 0.0529.*V + 0.0497;


% filters: 
% moving median - smoothing signal

% deriv1mraw = diff(mraw);
% indices = find(deriv1mraw<0.005 & deriv1mraw>-0.005);
% mraw = mraw(indices);
% t = t(indices);
m = movmedian(mraw,5); 
% m = medfilt1(mraw,50);

% m = sgolayfilt(mraw, 5, 61); golay
% fc = 7.3;
% Tc = 1./fc;
% [m,~] = lowpasf(mraw,1.0/128.0,Tc); low pass

dist_SWL = zeros(1,ngauge);

% if depth data is at the begining of the file
if nargin > 5
    fprintf('(d) \n \n')
    for i = 1:ngauge
%         swl = eta(t<tSWL,i);
%         depth(i) = mean(eta(swl(~isnan(swl)),i));
        dist_SWL(i) = mean(m(t<tSWL,i));
        % steven added this to incoorporate cases where depth data is
        % obtained otherwise
        if isnan(dist_SWL(i))
            dist_SWL(i) = 0;
            fprintf(1,'no depth value is interpreted, please incoorporate this another way.')
        end
        m(:,i) = dist_SWL(i) - m(:,i);  % subtract offset zero elevation value from data
    end
end

S_m = spectf(m,1./sampleRate,8);

end
