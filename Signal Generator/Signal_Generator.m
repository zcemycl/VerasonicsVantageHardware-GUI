% Aim: To produce a GUI to control multiple channels in terms of 
%      ultrasound emission.
%
% File name: Signal_Generator.m
%
% Last update:
% 02/07/2018

clear;

% ---------------------------------------------------------------------- %
% DEFAULT PARAMETERS
% ---------------------------------------------------------------------- %

% starting values for signal output
PULSE_REPETITION                    = 1e-3;     % [s]
TIME_DELAY                          = 0;        % [s]
FREQUENCY                           = 6e6;      % [Hz]
NUMBER_CYCLES                       = 3;
MAXIMUM_VOLTAGE                     = 30;       % [V]
STARTING_VOLTAGE                    = 10;       % [V]
TIME_AXIS_MAXIMUM                   = 5e-6;     % [s]
POLARITY                            = -1;       % it is a positive polarity

% literals needed for verasonics
PARAMETRIC_HALF_CYCLE_ON_TIME_FRAC  = 0.67;


% upper limit and lower limit of each edit text
% Number of cycles
ncpmax                              = 20 ; 
ncpmin                              = 1;
% Pulse repetition period (us)
prrmax                              = 4190000; 
prrmin                              = 10;
% Frequency (MHz)
fremax                              = 41; 
fremin                              = 0.5;
% Voltage (V)
volmax                              = 30; 
volmin                              = 0;
% Delay (us)
delmax                              = 45; 
delmin                              = 0;

% ---------------------------------------------------------------------- %
% Parameters under control
% ---------------------------------------------------------------------- %

% start with all elements turned off
initialize = zeros(1,32);

% calculate delay in wavelengths
delay = TIME_DELAY * FREQUENCY; % [wavelengths]

% assign additional values to workspace (modified by GUI)
pulse_repetition = PULSE_REPETITION * 1e6; % [us]
frequency = FREQUENCY * 1e-6; % [MHz]
number_cycles = NUMBER_CYCLES;
xmax = TIME_AXIS_MAXIMUM * 1e6;
polarity = POLARITY;
% set initial apodisation and delays
txapod = zeros(1, 32);
txdelay = delay .* ones(1, 32);

% ---------------------------------------------------------------------- %

% filename to save VSX variables
FILENAME = [tempdir, 'Verasonics_Leo_SGIV'];

% ---------------------------------------------------------------------- %
% Resource Parameters (Define your connector type and nature
% ---------------------------------------------------------------------- %

% specify system parameters (verasonics connector board)
Resource.Parameters.Connector       = 3;
Resource.Parameters.numTransmit     = 32;      % no. of transmit channels (2 brds).
Resource.Parameters.numRcvChannels  = 32;    % no. of receive channels (2 brds).
Resource.Parameters.speedOfSound    = 1485;    % speed of sound in m/sec
Resource.Parameters.verbose         = 2;
Resource.Parameters.initializeOnly  = 0;
Resource.Parameters.simulateMode    = 0;       % runs script in simulate mode
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.GUI = 'gui_signal_generator';

% ---------------------------------------------------------------------- %
% Transducer Details
% ---------------------------------------------------------------------- %

% Specify Trans structure array.
Trans.name = 'custom';
Trans.id = 100;
Trans.frequency = frequency; 
Trans.type = 0; 
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans.numelements = 32; 
Trans.ElementPos = zeros(Trans.numelements, 5);
Trans.ElementSens = 1;
Trans.maxHighVoltage = MAXIMUM_VOLTAGE;
Trans.Bandwidth = [0, 40];
Trans.impedance = 50;
Trans.Connector = (1:Trans.numelements)' ; 
Trans.connType = 6; 

% ---------------------------------------------------------------------- %
% Transmit Waveform TW object Details
% ---------------------------------------------------------------------- %

% specify the transmit waveform structure. 
TW(1).type = 'parametric';
TW(1).Parameters = [frequency, PARAMETRIC_HALF_CYCLE_ON_TIME_FRAC, 2*number_cycles, polarity];   % A, B, C, D

% bind the waveform (TW) the transmit object (TX)
TX(1).waveform = 1;                 % binds TW(1) to TX(1)

% set the apodization to zero
TX(1).Apod = txapod;

% turn off focusing
TX(1).focus = 0;

% set the time delay (given in fraction of a period)
TX(1).Delay = txdelay;

% ---------------------------------------------------------------------- %
% specify sequence events
% ---------------------------------------------------------------------- %

% specify the pulse repetition rate 
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = pulse_repetition;

% loop back to the SeqControl(1)
SeqControl(2).command = 'jump';
SeqControl(2).argument = 1;

% trigger out
SeqControl(3).command = 'triggerOut';
   
% event 1: transmit with trigger out
Event(1).info = 'Transmission.';
Event(1).tx = 1;                % use 1st TX structure.
Event(1).rcv = 0;               % no receiving
Event(1).recon = 0;             % no reconstruction.
Event(1).process = 0;           % no processing
Event(1).seqControl = [3, 1];   % wait for specified time

% after transmit are acquired do a jump back to the beginning
Event(2).info = 'Jump back to Event 1.'; % when VSX executes it returns
Event(2).tx = 0;         % no TX structure.
Event(2).rcv = 0;        % no Rcv structure.
Event(2).recon = 0;      % no reconstruction.
Event(2).process = 0;    % no processing
Event(2).seqControl = 2; % jump back to Event 1.

% save all the structures to a .mat file.
save(FILENAME);

% run VSX
filename = FILENAME;
VSX;