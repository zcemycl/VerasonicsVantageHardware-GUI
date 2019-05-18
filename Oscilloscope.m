% Aim: To produce a GUI to control multiple channels in terms of 
%      receive function
%
% File name: Oscilloscope.m
%
% Last update:
% 03/23/2018
clear
% ---------------------------------------------------------------------- %
% Parameters that will not be changed
% ---------------------------------------------------------------------- %
% literals needed for verasonics (all in normal units)
NUMBER_CYCLES                       = 3;
PULSE_REPETITION                    = 1e-3;             % [s]
PARAMETRIC_HALF_CYCLE_ON_TIME_FRAC  = 0.67;
PARAMETRIC_START_POLARITY           = 1;
MAXIMUM_VOLTAGE                     = 30;               % [V]
FREQUENCY                           = 5e6;              % [Hz]
TXAPOD                              = ones(1,32);
TXDELAY                             = zeros(1,32);

% define receive defaults used to set the receive buffer
RECEIVE_ACQ_START_TIME              = 0;                % [s] THIS CRASHES WINDOWS IF SET > 10e-6
RECEIVE_ACQ_END_TIME                = 500e-6;           % [s]

% number of traces to store in output data
NUM_TRACES_TO_STORE                 = 1000;

% number of frames
NUM_FRAMES                          = 1;

% set the default sampling mode
%   'NS200BW'   Nyquist sampling (200% bandwidth) of demodFrequency.
%   'NS200BWI'  2-1 interleaved sampling (requires 2 acquisitions).
%   'BS100BW'   100% bandwidth sampling of demodFrequency.
%   'BS67BW'    67% bandwidth sampling of demodFrequency.
%   'BS50BW'    50% bandwidth sampling of demodFrequency.
%   'custom'    sample rate set by decimSampleRate
SAMPLING_MODE                       = 'custom';

% set the sampling frequency (must be 250/n, n = 4, 5, 6, 7, ..., 20)
SAMPLING_FREQ_MHZ                   = 12.5;             % [MHz]
% ---------------------------------------------------------------------- %
% Parameters under control
% ---------------------------------------------------------------------- %
% Number of averages
number_of_averages                  = 10;
% rcv amplification 
rcvapod                             = ones(1,32);
% rcv on/off 
initialize                          = [1, zeros(1,31)];
% used to change the plot time axis
% Startpoint of the time axis
Time_Window_Offset                  = 0;
% Range of the time axis
Time_Window_Length                  = 5e-6;

% voltage-axis
amplitude                           = 30 ;
% auto (0) or fixed (1)
yaxisstate                          = 0;

% ---------------------------------------------------------------------- %
% Upper and Lower limit of each parameters
% ---------------------------------------------------------------------- %
% y-axis control
amplitudemax                        = 10000;
amplitudemin                        = 1;
% Number of Averages
noamax                              = 1000; 
noamin                              = 1;
% Delay 
delmax                              = 100e-6; 
delmin                              = 0;
% Seconds per Division 
spdmax                              = RECEIVE_ACQ_END_TIME; 
spdmin                              = 1e-9;
% Sampling frequency
fremax                              = 20 ; 
fremin                              = 1;
% Amplification
ampmax                              = 4 ; 
ampmin                              = -4;

% ---------------------------------------------------------------------- %
% filename to save VSX variables
FILENAME = [tempdir, 'Verasonics_Leo_SGIV'];
% ---------------------------------------------------------------------- %
% Resource Parameters (Define your connector type and nature
% ---------------------------------------------------------------------- %
% specify system parameters (verasonics connector board)
Resource.Parameters.Connector = 3;                  % the middle region of the connector board(LEMO)
Resource.Parameters.numTransmit = 32;               % no. of transmit channels (2 brds).
Resource.Parameters.numRcvChannels = 32;            % no. of receive channels (2 brds).
Resource.Parameters.speedOfSound = 1485;            % speed of sound in m/sec
Resource.Parameters.verbose = 2;                    % Warning levels
Resource.Parameters.initializeOnly = 0; 
Resource.Parameters.simulateMode = 0;               % runs script in simulate mode
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.GUI = 'gui_oscilloscope';
% ---------------------------------------------------------------------- %
% Transducer Details
% ---------------------------------------------------------------------- %
% Specify Trans structure array.
Trans.name = 'custom';
Trans.id = 100;
Trans.frequency = FREQUENCY*1e-6; 
Trans.type = 0; 
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans.numelements = 32 ; 
Trans.ElementPos = zeros(Trans.numelements, 5);
Trans.ElementSens = 1;
Trans.maxHighVoltage = 30 ;
Trans.Bandwidth = [0, 40];
Trans.impedance = 50;
Trans.Connector = (1:Trans.numelements)' ; 
Trans.connType = 6; 

% Calculation for TGC and receive
% compute number of samples recorded in RECEIVE_ACQ_END_TIME sampling at
% SAMPLING_FREQ_MHZ, adding a small additional factor
num_samples = round(1.1 * (RECEIVE_ACQ_END_TIME - RECEIVE_ACQ_START_TIME) * SAMPLING_FREQ_MHZ * 1e6);

% force to be an even number
if rem(num_samples, 2)
    num_samples = num_samples + 1;
end

% convert acquisition times to depth in wavelengths, dividing by 2 to
% account for round trip time
acq_start = FREQUENCY * RECEIVE_ACQ_START_TIME / 2;
acq_end   = FREQUENCY * RECEIVE_ACQ_END_TIME / 2;
% ---------------------------------------------------------------------- %
% Specify Resource Buffers
% ---------------------------------------------------------------------- %
Resource.RcvBuffer(1).datatype = 'int16'; % only data type in here
Resource.RcvBuffer(1).rowsPerFrame = num_samples;   % this allows for 1/4 maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = NUM_FRAMES;       % allocate 4 frames.
% ---------------------------------------------------------------------- %
% Transmit Waveform TW object Details
% ---------------------------------------------------------------------- %
% specify the transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [FREQUENCY *1e-6, PARAMETRIC_HALF_CYCLE_ON_TIME_FRAC, ...
    2*NUMBER_CYCLES, PARAMETRIC_START_POLARITY];    % A, B, C, D

% bind the waveform (TW) the transmit object (TX)
TX(1).waveform = 1;                                 % binds TW(1) to TX(1)
TX(1).Apod = zeros(1, Trans.numelements);
TX(1).focus = 0;
TX(1).Delay = zeros(1, Trans.numelements);
% ---------------------------------------------------------------------- %
% specify TGC Waveform structure
% ---------------------------------------------------------------------- %
TGC(1).CntrlPts = 0*ones(8,1);
TGC(1).rangeMax = acq_end;  % need to recalculate later
TGC(1).Waveform = computeTGCWaveform(TGC);
% ---------------------------------------------------------------------- %
% Specify Receive structure array 
% 1st for replace every n frames, where n the number of averages
% 2nd for accumulation of n-1 frames 
% ---------------------------------------------------------------------- %
for n = 1:2 
    Receive(n).Apod       = initialize.*rcvapod ;
    Receive(n).startDepth = acq_start; % start of data acquiaition in wavelength
    Receive(n).endDepth   = acq_end;  % end of data acq, wavevlengths
    Receive(n).TGC        =  1;
    if n == 1
        Receive(n).mode       =  0;         % for replace
    else
        Receive(n).mode       =  1;         % for accumulation 
    end
    Receive(n).bufnum     =  1;
    Receive(n).framenum   =  1;
    Receive(n).acqNum     =  1;
    Receive(n).ADCRate    =  SAMPLING_FREQ_MHZ;
    Receive(n).sampleMode =  SAMPLING_MODE;
end 
% ---------------------------------------------------------------------- %
% specify process events
% ---------------------------------------------------------------------- %
% Display settings process
% To plot the curves into the handle from the GUI axis
Process(1).classname = 'External';
Process(1).method = 'DisplayRF';
Process(1).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...       % process the most recent frame.
                         'dstbuffer','none'};

% ---------------------------------------------------------------------- %
% specify sequence events
% ---------------------------------------------------------------------- %
numAccum = number_of_averages; 
% SeqControl
% number of counts initialized later by loopTst
SeqControl(1).command = 'loopCnt';
SeqControl(1).argument = numAccum - 1;

% Conditionless jump back (unlike loopTst)
SeqControl(2).command = 'jump';
SeqControl(2).argument = 1;

% if loopCnt.argument > 0, it instructs a jump back to the event 
% specified by the loopTst.argument
SeqControl(3).command = 'loopTst';
SeqControl(3).argument = 3;

% Transfer data
SeqControl(4).command = 'transferToHost';

% Trigger in to synchronize the signal
SeqControl(5).command = 'triggerIn';
SeqControl(5).condition = 'Trigger_1_Rising';
SeqControl(5).argument = 0;
% ---------------------------------------------------------------------- %
% Event
% ---------------------------------------------------------------------- %
n = 1;
% 1: start to acquire the data using 1st rcv, and replace those previous
% ones
Event(n).info = 'Receive data -- Replace';
Event(n).tx = 0;
Event(n).rcv = 1;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 5;
n = n+1;
% 2: set the loop count for accumulation 
Event(n).info = 'Set loop count for number of accumulates.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;
% 3: accumulate the data using 2nd rcv
Event(n).info = 'Receive data -- Accumulates';
Event(n).tx = 0;
Event(n).rcv = 2;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 5;
n = n+1;
% 4: initilize the loop count, 
%    **jump back to event 3 and reduce the count number by 1 
%    if the count is greater than 0
%    **continue to the next event 5 if the count = 0
Event(n).info = 'Test loop count - if nz, jmp back to start of accumulates.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;
n = n+1;
% 5: Transfer all data from 1st and 2nd rcv
Event(n).info = 'Transfer data to host';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;
n = n+1;
% 6: Display using the external function for plotting
Event(n).info = 'External function for Display RF';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 1;
n = n+1;
% 7: Repeat, starting from event 1
Event(n).info = 'Jump to Event 1';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;


% ---------------------------------------------------------------------- %
% External functions
% ---------------------------------------------------------------------- %
% External function  -- DisplayRF
EF(1).Function = text2cell('%EF#1');
% ---------------------------------------------------------------------- %

% create buffer to store the time traces, and current index
rf_data_array = zeros(num_samples, NUM_TRACES_TO_STORE, 'int16');
rf_data_array_index = 1;

% ---------------------------------------------------------------------- %
% save all the structures to a .mat file.
save(FILENAME);

% run VSX
filename = FILENAME;
VSX;
return

% ---------------------------------------------------------------------- %
% Defining External functions
% ---------------------------------------------------------------------- %
%EF#1
DisplayRF(RData)

% This assignin allows the data could be saved by the save button
% during execution
assignin('base','RData',RData);

% read data from workspace
yaxisstate          = evalin('base', 'yaxisstate');
amplitude           = evalin('base', 'amplitude');
receive_data        = evalin('base', 'Receive');
rf_data_array       = evalin('base', 'rf_data_array');
rf_data_array_index = evalin('base', 'rf_data_array_index');
num_traces_to_store = evalin('base', 'NUM_TRACES_TO_STORE');
dt                  = 1 ./ (1e6 * evalin('base', 'SAMPLING_FREQ_MHZ'));
acq_start_time      = evalin('base', 'RECEIVE_ACQ_START_TIME');
Time_Window_Offset  = evalin('base', 'Time_Window_Offset');
Time_Window_Length  = evalin('base', 'Time_Window_Length');
initialize          = evalin('base', 'initialize');
number_of_averages  = evalin('base', 'number_of_averages');

% copy the data to the data array
rf_data_array(:, rf_data_array_index) = RData(:, 4);
rf_data_array_index = rf_data_array_index + 1;

% loop if buffer is full
if rf_data_array_index > num_traces_to_store
    rf_data_array_index = 1;
end

% replace values in workspace
assignin('base', 'rf_data_array', rf_data_array); 
assignin('base', 'rf_data_array_index', rf_data_array_index);

% create plot axis
time_axis = (acq_start_time + (0:receive_data(1).endSample - 1) * dt) * 1e6;
assignin('base','time_axis',time_axis);

acq_axis  = 1:num_traces_to_store;

% obtain the GUI handle
obj = findobj('tag','UI');
% the last number of the handles from the GUI is from the plotting axes
num_children = length(obj.Children);
receive_data = evalin('base', 'Receive');
ax = obj.Children(num_children);

% Plot the signals if the channels are on, 
% erase the signals plotting if off
legend_str = {};
if initialize(1)
    plot(ax, time_axis, RData(1:receive_data(1).endSample,1)./(number_of_averages), 'r'); 
    hold(ax, 'on');
    legend_str = [legend_str, {'Channel 1'}];
end
if initialize(2)
    plot(ax, time_axis, RData(1:receive_data(1).endSample,2)./(number_of_averages), 'k');
    hold(ax, 'on');
    legend_str = [legend_str, {'Channel 2'}];
end
if initialize(3)
    plot(ax, time_axis, RData(1:receive_data(1).endSample,3)./(number_of_averages), 'g');
    hold(ax, 'on');
    legend_str = [legend_str, {'Channel 3'}];
end
if initialize(4)
    plot(ax, time_axis, RData(1:receive_data(1).endSample,4)./(number_of_averages), 'b'); 
    hold(ax, 'on');
    legend_str = [legend_str, {'Channel 4'}];
end
% Grid
grid on
grid minor
% Legend to distinguish the signals from each channel
legend('show');
legend('boxoff');
% X-Y axis label
ylabel('Voltage(a.u.)');
xlabel('Time(us)');
legend(legend_str);
hold(ax, 'off');

% ======================================================================= %
% Voltage axis
% ======================================================================= %
% Follow the amplitude of the received signals to adjust the amplitude axis
axis(obj.Children(num_children),     'tight');
% If the amplitude axis state toggle button is pressed, 
% it fixes the axis to the amplitude value from the workspace
if yaxisstate  == 1
    set(obj.Children(num_children), 'YLim', [-amplitude, amplitude]);
end

% ======================================================================= %
% Time axis
% ======================================================================= %
% Follow the values of Time_Window_Offset and Time_Window_Length,
% to adjust the startpoint and range of the time axis
set(obj.Children(num_children), 'XLim', ...
    [Time_Window_Offset, Time_Window_Offset + Time_Window_Length] * 1e6);

drawnow limitrate;
%EF#1
