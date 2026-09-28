% 1. Split asdfGlobal into asdfWindowed
% 2. Recalculate sig_te and sig_del
% 3. Verify the analysis using "old data", which contains Dij, Pij for
% global asdf, then perform the analysis on the windowed ASDFs

% Shortcut: use global delay matrix from old data
%
%
% %% Windowing procedure
% Here we take the simplest method of dividing the recording time into a
% fixed number of windows. This is accomplished by dividing the asdf file
% into windows, such that c-webs do not cross into into the next window,
% i.e., we break the asdfGlobal file into separate asdfWindowed files, each
% of length NT/nWindows. We will examine more complex methods where c-webs
% may overrun their windows and c-webs corresponding to different windows
% may occur simultaneously.
%
% %% Data Formatting
% Spike times are formatted as uint32 (range: [0,2^32-1]) for SSC3 data 
% since the 1-hr recordings are binned into 1-ms bins, spike times occur
% within the range [1,3.6e+06-1].
%
% Neuronal proximity is formatted as uint16 since the dimensions of the
% 512-electrode array used in SSC3 is approximately 2 mm by 1 mm, the
% neuron cell body positions are stored in microns, and hence the neurons
% detected should fall within the range [0,2^16-1].
%
% \rvwg2019

% Full analysis, from scratch:
%% load the data set
dataset = '17_05-1';
load(strcat('DataSet',dataset))

N = data.nNeurons;
NT = data.recordinglength;

temp = data.timescale; %binning timescale, in ms
btime = str2double(temp(1:find(temp==' ')-1));

temp = data.samplingfrequency; %sampling rate, in kHz
fs = str2double(temp(1:find(temp==' ')-1));
clear temp

Xij = data.x;   %2d location of neuron cell bodies
Yij = data.y;

%% prepare analysis parameters
nWindows = 10;

%% prepare global asdf file
asdfGlobal = cell(N+2,1);
asdfGlobal{end-1} = btime;
asdfGlobal{end} = [N,NT];
numbursts = zeros(N,1); %stores the number of times we have multiple spikes
                        %in a single bin per neuron

% if NT<2^8
%     p = 'uint8';
% elseif NT<2^16
%     p = 'uint16';
% elseif NT<2^32
%     p = 'uint32';
% else
%     p = 'double';
% end

for n=1:N
    A = data.spikes{n};
    asdfGlobal{n} = unique(cast(A,'uint32'));
    
    numbursts(n) = numel(A)-numel(asdfGlobal{n});
    
    if numel(find(diff(asdfGlobal{n})==0))~=0
        error('Spike times are not increasing regularly.')
    end
end
clear data A

%% prepare windowed asdf files

windowBorders = (0:NT/nWindows:NT)';

w = 1;  %window ID:
asdfWindowed = cellfun(@(x) x(x>windowBorders(w) & ...
    x<=windowBorders(w+1)),asdfGlobal(1:N),'un',0);
asdfWindowed{N+1} = btime;
asdfWindowed{N+2} = [N,NT/nWindows];

%% prepare neuronal proximity
Rij = zeros(N);   %neuron proximity matrix (in micrometers)

for i=1:N
    for j=1:i-1
        Rij(j,i) = sqrt((Xij(j)-Xij(i))^2+(Yij(j)-Yij(i))^2);
    end
end
clear Xij Yij
Rij = cast(Rij'+Rij,'uint16');

%% avalanche shapes
tic
[aShapes,aLengths,aSizes] = asdf2Shapes(asdfGlobal);
toc

%% TE analysis

    

