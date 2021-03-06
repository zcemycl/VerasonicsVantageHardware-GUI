function gui_oscilloscope
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% VSX_GUI  This gui is opened by the main VSX program and allows control of various
% acquisition and processing parameters.

% Last update: 03/23/2018

% The base drive voltage
% this will be multiplied with the apod arrays to give the actual voltage 
% of each element
hv = evalin('base','MAXIMUM_VOLTAGE');
% an utility function for control of base voltage
% ** can be found in Vantage/System/P_Files
setTpcProfileHighVoltage(hv,1);
% =========================================================================
% Parameters under control
% =========================================================================
Time_Window_Length = evalin('base','Time_Window_Length'); 
Time_Window_Offset = evalin('base','Time_Window_Offset');
if Time_Window_Offset > Time_Window_Length
    disp('Time_Window_Offest cannot be greater than Time_Window_Length');
end
number_of_averages = evalin('base','number_of_averages');
initialize = evalin('base','initialize');
rcvapod = evalin('base','rcvapod');
% voltage-axis
amplitude = evalin('base','amplitude') ;% from the signal generator [V]
% auto (0) or fixed (1)
yaxisstate = evalin('base','yaxisstate');
% =========================================================================
% Upper and Lower Limit
% =========================================================================
% y-axis amplitude
amplitudemax = evalin('base','amplitudemax');
amplitudemin = evalin('base','amplitudemin');
% Number of Averages
noamax = evalin('base','noamax'); noamin = evalin('base','noamin');
% Time_Window_Offest 
delmax = evalin('base','delmax'); delmin = evalin('base','delmin');
% Time Window Length 
spdmax0 = evalin('base','spdmax'); spdmin0 = evalin('base','spdmin');
% Amplification
ampmax = evalin('base','ampmax'); ampmin = evalin('base','ampmin');
% =========================================================================
% CHECK VERASONICS STATUS
% =========================================================================
% close any previously opened GUI windows
delete(findobj('tag', 'UI'));

% get status of hardware
VDAS = 0;
if evalin('base', 'exist(''VDAS'',''var'')')
    VDAS = evalin('base', 'VDAS');
end

% don't allow the GUI to run without hardware
if VDAS ~= 1
    error('Verasonics hardware not connected.');
end

% determine simulateMode value
simMode = 0; 
if evalin('base','exist(''Resource'',''var'')')
    if evalin('base','isfield(Resource,''Parameters'')')
        if evalin('base','isfield(Resource.Parameters,''simulateMode'')')
            simMode = evalin('base','Resource.Parameters.simulateMode');
        end
    end
end

% don't allow the GUI to run in simulation mode
if simMode ~= 0
    error('GUI can''t be run in simulate mode');
end

% =========================================================================
% CREATE FIGURE WINDOW
% =========================================================================
% Close any previously opened GUI windows.
delete(findobj('tag','UI'));
% Initialize and hide the GUI as it is being constructed.
ScrnSize = get(0,'ScreenSize');
Bkgrnd = [1 1 1];
f = figure('Visible','off',...  %'Units','normalized',...
    'Position',[300,120,700,550],... %'Position',[0.7,0.25,0.25,0.50],...
    'Name','Oscilloscope',...
    'NumberTitle','off',...
    'MenuBar','none', ...
    'Resize','on', ...
    'tag','UI');
set(f,'CloseRequestFcn',{@closefunc});
set(f,'DefaultUicontrolBackgroundColor',Bkgrnd)

% Make the GUI visible, unless the call has requested that it be hidden.
set(f, 'Visible', 'on');
% ====================================================================== %
% Self-customized Parts
% ====================================================================== %
% UI objects setup
% Axis to plot transmit waveform
ax = axes('Units','normalized',...
    'Position',[0.07,0.48,0.47,0.5],...
    'tag','Axes',...
    'Box','on',...
    'Nextplot','replacechildren');
ylabel('Voltage(a.u.)');
xlabel('Time(us)');
grid on
grid minor
legend('show');
legend('boxoff');

% ======================================================================= %
% UIs appearance settings
% Frame: to organise UIs that belong to the same parameter control
% Text:  to notify the user of what the edit box is controlling
% Edit:  to edit the number of the parameter
% Limit: the lower and upper limits that are placed on the left or right
%        of the edit, to show the limits of the respective edit.
% ======================================================================= %
% General Settings
% ======================================================================= %
% Frame
gen_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.545,0.48,0.45,0.5]);
% Text
gen_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','General Settings',...
    'Units','normalized',...
    'Position',[0.57,0.88,0.25,0.08],...
    'FontUnits','normalized',...
    'FontSize',0.5);
% ======================================================================= %
% Time Window Length (spd)
% ======================================================================= %
% Frame
spdfram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.553,0.8,0.435,0.1]);
% Text
spdtext = uicontrol('Parent',f,...
    'Style','text',...
    'string','Time Window Length(us)',...
    'Units','normalized',...
    'Position',[0.555,0.82,0.1875,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4);
% Limit
spdup = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(spdmax0/1e-6),...
    'Units','normalized',...
    'Position',[0.935,0.8425,0.05,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
spdlow = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(spdmin0/1e-6),...
    'Units','normalized',...
    'Position',[0.74,0.8425,0.04,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
% Edit
spdedit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(Time_Window_Length/1e-6),...
    'Units','normalized',...
    'Position',[0.785,0.84,0.14925,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@spdedit_callback});
% ======================================================================= %
% Number of Averages (noa)
% ======================================================================= %
% Frame
noafram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.553,0.7,0.435,0.1]);
% Text
noatext = uicontrol('Parent',f,...
    'Style','text',...
    'string','Number of Averages',...
    'Units','normalized',...
    'Position',[0.57,0.72,0.18,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4);
% Limit
noaup = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(noamax),...
    'Units','normalized',...
    'Position',[0.945,0.7425,0.04,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
noalow = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(noamin),...
    'Units','normalized',...
    'Position',[0.745,0.7425,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
% Edit
noaedit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(number_of_averages),...
    'Units','normalized',...
    'Position',[0.785,0.74,0.15,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@noaedit_callback});
% ======================================================================= %
% Time_Window_Offest 
% ======================================================================= %
% Frame
dtfram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.553,0.6,0.435,0.1]);
% Text
dttext = uicontrol('Parent',f,...
    'Style','text',...
    'string','Time Window Offset(us)',...
    'Units','normalized',...
    'Position',[0.57,0.62,0.18,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4);
% Limit
dtup = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(delmax/1e-6),...
    'Units','normalized',...
    'Position',[0.945,0.6425,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
dtlow = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(delmin),...
    'Units','normalized',...
    'Position',[0.745,0.6425,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
% Edit
dtedit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(Time_Window_Offset/1e-6),...
    'Units','normalized',...
    'Position',[0.785,0.64,0.15,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@dtedit_callback});
% ======================================================================= %
% auto-or-fix voltage axis
% ======================================================================= %
% Frame
yaxisfram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.553,0.5,0.435,0.1]);
% If the toggle is unpressed, the voltage axis is in auto mode
% If the toggle is pressed, the voltage axis is in fixed mode
if yaxisstate == 0
    str_y = 'AUTO';
else
    str_y = 'FIXED';
end
% Text besides the toggle button
yaxistxt = uicontrol('Parent',f,...
    'Style','text',...
    'string','Voltage axis (a.u.)',...
    'Units','normalized',...
    'Position',[0.57,0.52,0.18,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4);
% Text in the toggle button representing the state
ytogglebutton = uicontrol('Parent',f,...
    'Style','togglebutton',...
    'string',str_y,...
    'Value',yaxisstate,...
    'Units','normalized',...
    'Position',[0.8,0.51,0.12,0.035],...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Callback',{@ytoggle_callback});
% Limit
yup = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(amplitudemax),...
    'Units','normalized',...
    'Position',[0.945,0.5525,0.04,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
ylow = uicontrol('Parent',f,...
    'Style','text',...
    'string',num2str(amplitudemin),...
    'Units','normalized',...
    'Position',[0.745,0.5525,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
% Edit to input the value that will be used in fixed mode
yedit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(amplitude),...
    'Units','normalized',...
    'Position',[0.785,0.55,0.15,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@yedit_callback});
% ======================================================================= %
% Save button
% ======================================================================= %
save_gui_button = uicontrol('Parent',f,...
    'Style','pushbutton',...
    'Units','normalized',...
    'string','Save: GUI',...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Position',[0.6,0.42,0.15,0.05],...
    'Callback',{@save_gui_callback});
save_data_button = uicontrol('Parent',f,...
    'Style','pushbutton',...
    'Units','normalized',...
    'string','Save: Data',...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Position',[0.8,0.42,0.15,0.05],...
    'Callback',{@save_data_callback});
% ======================================================================= %
% Individual Settings 
% ======================================================================= %
% Frame
ind_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.01,0.05,0.96,0.3]);
% Text
ind_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Individual Settings',...
    'Units','normalized',...
    'Position',[0.02,0.25,0.3,0.08],...
    'FontUnits','normalized',...
    'FontSize',0.5);
% ======================================================================= %
% Channel 1 control
% ======================================================================= %
fram_1 = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.02,0.07,0.22,0.2]);
ch_1_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.03,0.19,0.195,0.07]);
ch_1_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Channel 1',...
    'Units','normalized',...
    'Position',[0.035,0.195,0.15,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.6);
% ======================================================================= %
% On/Off togglebutton 1
% ======================================================================= %
if initialize(1) == 0
    str_1 = 'OFF';
else 
    str_1 = 'ON';
end
onoff1_toggle = uicontrol('Parent',f,...
    'Style','togglebutton',...
    'string',str_1,...
    'Value',initialize(1),...
    'Units','normalized',...
    'Position',[0.17,0.2,0.04,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Callback',{@onoff_1_callback});
% ======================================================================= %
% Amplification 1
% ======================================================================= %
V_1_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Amplification:',...
    'Units','normalized',...
    'Position',[0.025,0.12,0.16,0.06],...
    'FontUnits','normalized',...
    'FontSize',0.4);
voltup1 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmax),...
    'Units','normalized',...
    'Position',[0.175,0.12,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
voltlow1 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmin),...
    'Units','normalized',...
    'Position',[0.057,0.12,0.05,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
V_1_edit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(rcvapod(1)),...
    'Units','normalized',...
    'Position',[0.09,0.12,0.085,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@v1edit_callback});
% ======================================================================= %
% Channel 2 control
% ======================================================================= %
fram_2 = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.02+0.24,0.07,0.22,0.2]);
ch_2_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.03+0.24,0.19,0.195,0.07]);
ch_2_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Channel 2',...
    'Units','normalized',...
    'Position',[0.035+0.24,0.195,0.15,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.6);
% ======================================================================= %
% On/Off togglebutton 2
% ======================================================================= %
if initialize(2) == 0
    str_2 = 'OFF';
else 
    str_2 = 'ON';
end
onoff2_toggle = uicontrol('Parent',f,...
    'Style','togglebutton',...
    'string',str_2,...
    'Value',initialize(2),...
    'Units','normalized',...
    'Position',[0.17+0.24,0.2,0.04,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Callback',{@onoff_2_callback});
% ======================================================================= %
% Amplification 2
% ======================================================================= %
V_2_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Amplification:',...
    'Units','normalized',...
    'Position',[0.025+0.24,0.12,0.16,0.06],...
    'FontUnits','normalized',...
    'FontSize',0.4);
voltup2 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmax),...
    'Units','normalized',...
    'Position',[0.175 + 0.24,0.12,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
voltlow2 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmin),...
    'Units','normalized',...
    'Position',[0.057 + 0.24,0.12,0.05,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
V_2_edit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(rcvapod(2)),...
    'Units','normalized',...
    'Position',[0.09+0.24,0.12,0.085,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@v2edit_callback});
% ======================================================================= %
% Channel 3 control
% ======================================================================= %
fram_3 = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.02+0.48,0.07,0.22,0.2]);
ch_3_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.03+0.48,0.19,0.195,0.07]);
ch_3_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Channel 3',...
    'Units','normalized',...
    'Position',[0.035+0.48,0.195,0.15,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.6);
% ======================================================================= %
% On/Off togglebutton 3
% ======================================================================= %
if initialize(3) == 0
    str_3 = 'OFF';
else 
    str_3 = 'ON';
end
onoff3_toggle = uicontrol('Parent',f,...
    'Style','togglebutton',...
    'string',str_3,...
    'Value',initialize(3),...
    'Units','normalized',...
    'Position',[0.17+0.48,0.2,0.04,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Callback',{@onoff_3_callback});
% ======================================================================= %
% Amplification 3
% ======================================================================= %
V_3_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Amplification:',...
    'Units','normalized',...
    'Position',[0.025+0.48,0.12,0.16,0.06],...
    'FontUnits','normalized',...
    'FontSize',0.4);
voltup3 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmax),...
    'Units','normalized',...
    'Position',[0.175+0.48,0.12,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
voltlow3 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmin),...
    'Units','normalized',...
    'Position',[0.057+0.48,0.12,0.05,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
V_3_edit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(rcvapod(3)),...
    'Units','normalized',...
    'Position',[0.09+0.48,0.12,0.085,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@v3edit_callback});
% ======================================================================= %
% Channel 4 control
% ======================================================================= %
fram_4 = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.02+0.72,0.07,0.22,0.2]);
ch_4_fram = uicontrol('Parent',f,...
    'Style','frame',...
    'Units','normalized',...
    'Position',[0.03+0.72,0.19,0.195,0.07]);
ch_4_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Channel 4',...
    'Units','normalized',...
    'Position',[0.035+0.72,0.195,0.15,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.6);
% ======================================================================= %
% On/Off togglebutton 4
% ======================================================================= %
if initialize(4) == 0
    str_4 = 'OFF';
else 
    str_4 = 'ON';
end
onoff4_toggle = uicontrol('Parent',f,...
    'Style','togglebutton',...
    'string',str_4,...
    'Value',initialize(4),...
    'Units','normalized',...
    'Position',[0.17+0.72,0.2,0.04,0.05],...
    'FontUnits','normalized',...
    'FontSize',0.4,...
    'Callback',{@onoff_4_callback});
% ======================================================================= %
% Amplification 4
% ======================================================================= %
V_4_text = uicontrol('Parent',f,...
    'Style','text',...
    'string','Amplification:',...
    'Units','normalized',...
    'Position',[0.025+0.72,0.12,0.16,0.06],...
    'FontUnits','normalized',...
    'FontSize',0.4);
voltup4 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmax),...
    'Units','normalized',...
    'Position',[0.175+0.72,0.12,0.03,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
voltlow4 = uicontrol('Parent',f,...
    'Style','text',...
    'String',num2str(ampmin),...
    'Units','normalized',...
    'Position',[0.057+0.72,0.12,0.05,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.8);
V_4_edit = uicontrol('Parent',f,...
    'Style','edit',...
    'string',num2str(rcvapod(4)),...
    'Units','normalized',...
    'Position',[0.09+0.72,0.12,0.085,0.03],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'Callback',{@v4edit_callback});
% ======================================================================= %
% Callback
% ======================================================================= %
% General Settings
% ======================================================================= %
% ----------------
% Save
% ----------------
function save_gui_callback(hObject,eventdata)
    [filename, pathname] = uiputfile(...
        {'*.png'; '*.jpg'; '*.fig'}, ...
        'Save GUI picture as',...
        'GUI');
    saveas(gcf, fullfile(pathname, filename));
end
function save_data_callback(hObject,eventdata)
    RData = evalin('base','RData');
    time_axis = evalin('base','time_axis');
    [filename, pathname] = uiputfile(...
        {'*.mat'}, ...
        'Save Receive Data as',...
        'Data');
    save(fullfile(pathname,filename),'RData','time_axis');
end
% ------------------
% Time Window Length
% ------------------
function spdedit_callback(hObject,eventdata)
    % ================================================================ %
    % UI input
    % ================================================================ %
    spd = str2double(get(spdedit,'String'));
    % ================================================================ %
    % Condition
    % ================================================================ %
    [spd] = limit(spdedit,spd,spdmax0/1e-6,spdmin0/1e-6);
    % ================================================================ %
    % System connection
    % ================================================================ %
    assignin('base','Time_Window_Length',spd * 1e-6);
end
% ------------------
% Number of Averages
% ------------------
function noaedit_callback(hObject,eventdata)
    % ================================================================ %
    % UI input
    % ================================================================ %
    noa = str2double(get(noaedit,'String'));
    % ================================================================ %
    % Condition
    % ================================================================ %
    [noa] = limit(noaedit, noa,noamax,noamin);
    % ================================================================ %
    % System connection
    % ================================================================ %
    assignin('base','number_of_averages',noa);

    SeqControl = evalin('base','SeqControl');
    SeqControl(1).argument = noa - 1;
    assignin('base','SeqControl',SeqControl);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'SeqControl'};
    assignin('base','Control',Control);
end
% ------------------
% Time Window Offset
% ------------------
function dtedit_callback(hObject,eventdata)
    dt = str2double(get(dtedit,'String'));
    [dt] = limit(dtedit, dt,delmax/1e-6,delmin/1e-6);
    % Put back to number to workspace
    assignin('base','Time_Window_Offset',dt*1e-6);
end
% ------------------
% Voltage axis limit
% ------------------
function ytoggle_callback(hObject,eventdata)
    yaxisstate = get(ytogglebutton,'Value');
    assignin('base','yaxisstate',yaxisstate);
    if yaxisstate == 0
        str_y = 'AUTO';
    else
        str_y = 'FIXED';
    end
    set(ytogglebutton,'string',str_y);
end
function yedit_callback(hObject,eventdata)
    amplitude = str2double(get(yedit, 'string'));
    [amplitude] = limit(yedit,amplitude,amplitudemax,amplitudemin);
    assignin('base','amplitude',amplitude);
end
% ====================================================================== %
% Individual settings
% ====================================================================== %
% ------------------
% on/off
% ------------------
    function onoff_1_callback(hObject,eventdata)
        onoff(onoff1_toggle,1);
    end
    function onoff_2_callback(hObject,eventdata)
        onoff(onoff2_toggle,2);
    end
    function onoff_3_callback(hObject,eventdata)
        onoff(onoff3_toggle,3);
    end
    function onoff_4_callback(~,eventdata)
        onoff(onoff4_toggle,4);
    end
% ------------------
% Amplification
% ------------------
    function v1edit_callback(hObject,eventdata,V_1_slid)
        vedit(V_1_edit,1);
    end
    function v2edit_callback(hObject,eventdata,V_2_slid)
        vedit(V_2_edit,2);
    end
    function v3edit_callback(hObject,eventdata,V_3_slid)
        vedit(V_3_edit,3);
    end
    function v4edit_callback(hObject,eventdata,V_4_slid)
        vedit(V_4_edit,4);
    end
% ====================================================================== %
% Functions for streamling the codes
% ====================================================================== %
% On Off callback function for each channel
function onoff(onoff_handle,element_number)
    initialize = evalin('base','initialize');
    rcvapod = evalin('base','rcvapod');
    number_of_averages = evalin('base','number_of_averages');
    Receive = evalin('base','Receive');
    Control = evalin('base','Control');
    % ================================================================ %
    % UI input
    % ================================================================ %
    onoff_val = get(onoff_handle,'Value');
    % ================================================================ %
    % Condition
    % ================================================================ %
    if onoff_val > -1
        set(onoff_handle,'string','ON');
        initialize(element_number) = 1;
        assignin('base','initialize',initialize);
    else
        set(onoff_handle,'string','OFF');
        initialize(element_number) = 0;
        assignin('base','initialize',initialize);
    end
    % ================================================================ %
    % System connection
    % ================================================================ %
    for i = 1:2
        Receive(i).Apod = initialize.*rcvapod;
    end
    Receive(3) = Receive(2);
    Receive(3).mode = 0;
    assignin('base','Receive',Receive);
    Control.Command = 'update&Run';
    Control.Parameters = {'Receive'};
    assignin('base','Control',Control);
end

% Amplification callback function for each channel
function vedit(vedit_handle,element_number)
    number_of_averages = evalin('base','number_of_averages');
    rcvapod = evalin('base','rcvapod');
    % ================================================================ %
    % UI input
    % ================================================================ %
    volt = str2double(get(vedit_handle,'string'));
    % ================================================================ %
    % Condition
    % ================================================================ %
    [volt] = limit(vedit_handle, volt,ampmax,ampmin);
    % ================================================================ %
    % System connection
    % ================================================================ %
    rcvapod(element_number) = volt;
    assignin('base','rcvapod',rcvapod);
    initialize = evalin('base','initialize');
    Receive = evalin('base','Receive');
    Control = evalin('base','Control');
    intermediate_rcvapod = initialize.*rcvapod;
    for i = 1:2
        Receive(i).Apod = intermediate_rcvapod;
    end
    assignin('base','Receive',Receive);
    Control.Command = 'update&Run';
    Control.Parameters = {'Receive'};
    assignin('base','Control',Control);
end

% Upper limit & lower limit function
function [edit_get] = limit(edit_handle, edit_get,max,min)
    if edit_get > max
        edit_get = max;
        set(edit_handle,'string',num2str(max));
    elseif edit_get < min
        edit_get = min;
        set(edit_handle,'string',num2str(min));
    end
end
% ====================================================================== %
% From default GUI
% ====================================================================== %
% Read in TPC structure that has been initialized by VSX.
TPC = evalin('base','TPC');

% - For more than one active profile, determine the number to use for the 2nd slider. Profile 5 
%   has priority over profiles 2-4.
hv2 = 0;     % hv2 will get the profile to use for 2nd slider.
assignin('base', 'hv2GUIprofile', hv2);


% Make the GUI visible, unless the call has requested that it be hidden.
visibility = 'on';
if((true == evalin('base', 'exist(''Mcr_GuiHide'', ''var'')')) && (1 == evalin('base', 'Mcr_GuiHide')))
    % Caller has requested that we do NOT show the GUI window.
    visibility = 'off';
end
set(f,'Visible', visibility);

    function closefunc(source,eventdata)
        assignin('base', 'exit', 1);
        if exist('hvtmr','var'), stop(hvtmr); delete(hvtmr); end
        delete(f);
    end

% if saveRF is in the temp directory, delete it
if exist([tempdir,'saveRF.m'],'file')
    delete([tempdir,'saveRF.m']);
end

end

