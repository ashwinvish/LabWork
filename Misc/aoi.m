% modified version of aoi . . . for NI PCI-6289
% currently works with an xPC Target model
% written by Clay Didier
% csdidier@mit.edu
%
%Usage:
%   AOI ( DATA, CHANNELS, FREQUENCY ) is a routine that interfaces the
%   NI-DAC board throught the xpc target system, which outputs specified 
%   data on up to 4 analog channels, and receives input on 4 analog inputs.
%
%   DATA is the data matrix to be outputted by the NI PCI-6289.  It is an M
%   row by N column matrix, where M is the number of samples to be output
%   (hence defining the total execution time), and N must be a number from
%   1-4 (up to 4 channels).
%
%   CHANNELS is a row vector containing the enumeration of input channels
%   to be sampled.  The length of this vector must match the column length
%   of the DATA entry.  Channels must be specified in increasing order.
%   The four channels are specified starting at index 0, with 0 being the
%   first, and 3 being the last.
%
%   FREQUENCY is the sampling rate of the input/output signals in cycles/s.
%
%   EXAMPLE:
%	output_data_ch1 = [ 1 2 3 4 5 ]';
%       output_data_ch2 = [ 6 7 8 9 0 ]';
%       input_channels = [ 0 1 3 5 ];
%       frequency = 100;
%       input_data = aoi( [ output_data_ch1 output_data_ch2 ], input_channels, frequency); 
%

% Based on:
% wrapper for RT testpulse
% this is engineered to be an aoi replacement
% Written by Neville Sanjana, 6/2005
%
% Directions: Change target computer IP (and port, if necessary) to get
% target object
%
% Based on:
%$Revision 0.992b 08/21/2002
%written by: Alan Chen
%contact: alanchen@mit.edu
%
%

function current=aoi(varargin)

if isequal(size(varargin,2),3);
    DATA=varargin{1};
    CHANNELS=varargin{2};
    FREQUENCY=varargin{3};
    READCHANNELS=CHANNELS;
else 
    if isequal(size(varargin,2),4);
        DATA=varargin{1};
        CHANNELS=varargin{2};
        FREQUENCY=varargin{3};
        READCHANNELS=varargin{4};
    else
        error('Wrong number of input arguments!');
    end;
end; 


%BASIC ERROR CHECKING:
%first we need to make certain the dimensions of DATA and CHANNELS match
dataSize=size(DATA);
chanSize=size(CHANNELS);
readchanSize=size(READCHANNELS);
numberOfChannels=chanSize(2);
numberOfReadChannels=readchanSize(2);
numberOfSamples=dataSize(1);

%BASIC ERROR CHECKING:
%check to make sure the number of channels is 4 or less
if isequal(dataSize(2),chanSize(2));
else
    error('DATA dimensions must correspond with CHANNELS dimensions');
end;

%BASIC ERROR CHECKING:
%check to make sure the number of channels is 4 or less
if (chanSize(2)<=4);
else
    error('Too many channels specified . . . only up to 4 allowed!');
end

%BASIC ERROR CHECKING:
%check to make sure the channels are in range
if (ismember(CHANNELS,[0 1 2 3]));
else
    error('The channels must be member of {0, 1, 2, 3}');
end
if (ismember(READCHANNELS,[0 1 2 3]));
else
    error('The channels read must be member of {0, 1, 2, 3}');
end

%BASIC ERROR CHECKING:
%check to make sure the channels are in increasing order
if isequal(CHANNELS,sort(CHANNELS));
else
    error('The channels must be listed in increasing order!');
end

%TARGET SETUP:
%connect to xPC target and load bigmamav7a if necessary


targetIP = '10.10.10.11';
tg = xpctarget.xpc('TCPIP', targetIP, '22222');
if isequal(get(tg,'Application'),'bigmamav7o'); 
else
    load(tg,'bigmamav7o');
end

%TARGET SETUP:
% get needed parameter Ids
programSwitchParameter=getparamid(tg,'program index','Value');
channelIndexParameter=getparamid(tg,'aoi/xpcfromfilebyclay','P3');
dataReadWidthParameter=getparamid(tg,'aoi/xpcfromfilebyclay','P2');
ADchannelIndexParameter=getparamid(tg,'PCI6289AD','P1');
DAchannelIndexParameter=getparamid(tg,'PCI6289DAclaymod','P1');
%gainIndexParameter=getparamid(tg,'Gain','Gain');
%sensitivityIndexParameter=getparamid(tg,'Sensitivity','Gain');

%TARGET SETUP:
%create the channelIndex vector
channelIndex=ones(1,4).*(-1);

for i=1:numberOfChannels;
    channelIndex(CHANNELS(i)+1)=i;
end;

%readchannelIndex=ones(1,4).*(-1);
%for i=1:numberOfReadChannels;
%   readchannelIndex(READCHANNELS(i)+1)=READCHANNELS(i)+1;
%end;
readchannelIndex=[1 2 3 4];


%TARGET SETUP:
%set the parameters 
aoiProgramIndex = 1;
tg.setparam(programSwitchParameter, aoiProgramIndex);
tg.setparam(channelIndexParameter, channelIndex);
tg.setparam(dataReadWidthParameter, numberOfChannels);
tg.setparam(ADchannelIndexParameter, readchannelIndex);
tg.setparam(DAchannelIndexParameter, channelIndex);
%tg.setparam(gainIndexParameter, 1.0);
%tg.setparam(sensitivityIndexParameter, 1.0);


%TARGET SETUP:
% declare variables and setup xpc target
samplePeriod = 1/FREQUENCY;
stopTime = numberOfSamples * samplePeriod;
tg.StopTime=stopTime;
tg.SampleTime = samplePeriod;

%DATA SETUP AND TRANSFER:
% loop through channels and
% modify DATA variable to include all four
% channels  in order of 0, 1, 2, 3.
%tic;
DATA=DATA';
%disp('Data Conversion Done');
%toc

%DATA SETUP AND TRANSFER:
% save revised DATA variable to local file in byte format
xpcbytes2file('aoidata.dat',DATA);
%disp('File "aoidata.dat" saved to local disk');
%toc

%DATA SETUP AND TRANSFER:
% transfer data file to xpc target using FTP
xpcftp = xpctarget.ftp(tg);
xpcftp.put('aoidata.dat');
%disp('File "aoidata.dat" transferred to target');

%toc

%DATA SETUP AND TRANSFER:
% remove local data file
delete aoidata.dat;
%disp('Local file "aoidata.dat" removed');
%toc


%XPC TARGET EXECUTION:
tg.start;
%disp('XPC Target started');
%toc
while(strcmp(tg.status,'running'))
    pause(0.1);
end

output = tg.OutputLog;

%current=output(:,READCHANNELS(1)+1);
%current=output(2:end,READCHANNELS(1)+1);  %tg.OutputLog appears to have one more sample than DATA
%tested with cosine and analog inputs start one sample behind analog output
current=[];

%for i=2:numberOfReadChannels;
for i=1:numberOfReadChannels;
%current=[current output(:,READCHANNELS(i)+1)];
current=[current output(2:end,READCHANNELS(i)+1)];
end;



