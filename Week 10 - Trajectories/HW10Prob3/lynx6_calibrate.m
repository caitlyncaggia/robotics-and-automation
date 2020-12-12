%============================ lynx6_calibrate ============================
%
%  A script used to calibrate the lynx6 manipulator.  Saves the calibration 
%  information to the specified file.
%
%============================ lynx6_calibrate ============================

%
%  Name:		lynx6_calibrate.m
%
%  Author:		Patricio A. Vela,	pvela@gatech.edu
%
%  Created:		2013/06/10
%
%============================ piktul_calibrate ===========================
function lynx6_calibrate(parms)

limits = zeros([3 6]);
limits(1,:) = 500;
limits(2,:) = 1500;
limits(3,:) = 2500;

jointLabels = {'Shoulder 1','Shoulder 2','Elbow','Wrist 1','Wrist 2','Gripper'};
minStrings  = {'-90','-90','-90','-90','-90','1.125'};
zeroStrings = {'0','0','0','0','0','0.5'};
maxStrings  = {'90','90','90','90','90','0.125'};
manDim = 6;

if (nargin > 0)
  fnames = fieldnames(parms);
  for ii = 1:length(fnames)
    switch fnames{ii}
      case 'limits',
        if (size(parms.limits,2) == manDim)
          if (size(parms.limits,1) == 2)
  	        limits(1,:) = parms.limits(1,:);
  	        limits(3,:) = parms.limits(2,:);
	        limits(2,:) = mean(limits,1);
          end
	      if (size(parms.limits,1) == 3)
	        limits(1,:) = parms.limits(1,:);
	        limits(2,:) = parms.limits(2,:);
	        limits(3,:) = parms.limits(3,:);
          end
        end
      case 'musecLims',
          limits = parms.musecLims;
      case 'alphaLims',
          if (isfield(parms, 'alphaOrient'))
            orientSign = parms.alphaOrient;
          else
            orientSign = ones([1, manDim]);
          end
          for ii = 1:(manDim)
            if (orientSign(ii) > 0)
              minStrings{ii}  = num2str(parms.alphaLims(1,ii));
              zeroStrings{ii} = num2str(parms.alphaLims(2,ii));
              maxStrings{ii}  = num2str(parms.alphaLims(3,ii));
            else
              minStrings{ii}  = num2str(parms.alphaLims(3,ii));
              zeroStrings{ii} = num2str(parms.alphaLims(2,ii));
              maxStrings{ii}  = num2str(parms.alphaLims(1,ii));
            end
          end
    end
  end
end

%== GUI Interface setup. ==
gM = 600;
gN = 700;
doCts = false;					% Continuous update?

fh = figure('Visible','Off','Position',[100 100 gN gM]);

hdone = uicontrol(fh,'Style','pushbutton','String','Done',...
                  'Position',[gN-100,10,60,20],...
		  'Callback',{@doneCallback});

hgo   = uicontrol(fh,'Style','pushbutton','String','Goto',...
                  'Position',[30,10,60,20],...
		  'Callback',{@gotoCallback});

hcts   = uicontrol(fh,'Style','togglebutton','String','Continuous',...
                  'Position',[100,10,80,20],...
		  'Callback',{@ctsCallback});

hfile   = uicontrol(fh,'Style','edit','String','lynxparms.mat',...
                  'Position',[250,10,200,20]);

hsave = uicontrol(fh,'Style','pushbutton','String','Save',...
                  'Position',[gN-220,10,60,20],...
		  'Callback',{@saveCallback});


ebox = [40 20];
tbox = [40 15];

lablPos = [10, gM-30, 100, 20];
slidPos = [150, gM-30, 500, 20];
eminPos = [150 , gM-60, ebox(1), ebox(2)];
emaxPos = [565, gM-60, ebox(1), ebox(2)];
zeroPos = [350, gM-60, ebox(1), ebox(2)];
eminVal = [150 , gM-60, ebox(1), ebox(2)];
emaxVal = [565, gM-60, ebox(1), ebox(2)];
zeroVal = [350, gM-60, ebox(1), ebox(2)];
%textPos = [1+2*tbox(1)+20, gM-50+30, tbox(1), tbox(2)];
currPos = [slidPos(1)+slidPos(3)+10, slidPos(2), tbox(1), tbox(2)];
dSlider = -90;


for ii = 1:manDim
  %--[A] Slider bar and slider bar value (editable)
  cPos = lablPos + [0, dSlider*(ii-1), 0, 0];
  uicontrol(fh,'Style','text','String',jointLabels{ii},'Position',cPos);

  cPos = slidPos + [0, dSlider*(ii-1), 0, 0];
  hslider{ii} = uicontrol(fh,'Style','slider','UserData',ii, ...
                             'Max', limits(3,ii), 'Min', limits(1,ii), ...
                             'Value', limits(2,ii), ...
                             'SliderStep', [0.005 0.05], ...
                             'Position', cPos, 'Callback', {@sliderCallback});
  set(hslider{ii},'Units','Normalized');

  %--[B] Minimum allowed value for signal (plus it's actual joint value)
  cPos = eminPos + [0, dSlider*(ii-1), 0, 0];
  tminh{ii} = uicontrol(fh,'Style','pushbutton','String','Min:', ...
                   'UserData', ii, 'Callback', {@tminCallback}, ...
                   'Position',cPos);
  cPos = eminPos + [ebox(1), dSlider*(ii-1),0, 0];
  eminh{ii} = uicontrol(fh,'Style','edit','String',limits(1,ii),...
                 'Value', ii, 'Position', cPos, 'Callback', {@minCallback});
  cPos = cPos + [0, -eminPos(4)-5, 0, 0];
  minv{ii} = uicontrol(fh,'Style','edit','String',minStrings{ii}, ...
                 'Value', ii, 'Position', cPos);

  %--[C] Zero or middle value for signal (plus it's actual joint value)
  cPos = zeroPos + [0, dSlider*(ii-1), 0, 0];
  tzeroh{ii} = uicontrol(fh,'Style','pushbutton','String','Zero:', ...
                   'UserData', ii, 'Callback', {@tzeroCallback}, ...
                   'Position',cPos);
  cPos = zeroPos + [ebox(1), dSlider*(ii-1),0, 0];
  zeroh{ii} = uicontrol(fh,'Style','edit','String',limits(2,ii),...
                 'Value', ii, 'Position', cPos, 'Callback', {@zeroCallback});
  cPos = cPos + [0, -zeroPos(4)-5, 0, 0];
  zerov{ii} = uicontrol(fh,'Style','edit','String',zeroStrings{ii}, ...
                 'Value', ii, 'Position', cPos);

  %--[D] Maximum allowed value for signal (plus it's actual joint value)
  cPos = emaxPos + [0, dSlider*(ii-1), 0, 0];
  tmaxh{ii} = uicontrol(fh,'Style','pushbutton','String','Max:', ...
                   'UserData', ii, 'Callback', {@tmaxCallback}, ...
                   'Position',cPos);
  cPos = emaxPos + [ebox(1), dSlider*(ii-1),0, 0];
  emaxh{ii} = uicontrol(fh,'Style','edit','String',num2str(limits(3,ii)), ...
                   'Value', ii, ...
                   'Position', cPos, 'Callback', {@maxCallback});
  cPos = cPos + [0, -emaxPos(4)-5, 0, 0];
  maxv{ii} = uicontrol(fh,'Style','edit','String',maxStrings{ii}, ...
                 'Value', ii, 'Position', cPos);

  cPos = currPos + [0, dSlider*(ii-1), 0, 5];
  tcurh{ii} = uicontrol(fh,'Style','edit','String',num2str(limits(2,ii)), ...
                   'Position',cPos,'UserData',ii,'Callback',{@currCallback});
end


set(fh,'Visible','on');
set(fh,'Name','Manipulator Interface Setup');


arm = lynx6;
uiwait;
arm.shutdown();

%
%
%============================ Helper Functions ===========================
%
%

  %========================== Slider Functions =========================
  %

  %--------------------------- sliderCallback --------------------------
  function sliderCallback(source, eventdata)
    pwmsig = round(get(source,'Value'));
    index = get(source, 'UserData');
    set(tcurh{index},'String',num2str(pwmsig,'%4d'));
  
    if (doCts)
      gotoSliderValues(3);
    end
  end

  %--------------------------- setSliderValue --------------------------
  function setSliderValue(index, value)
    set(hslider{index},'Value',value);
    set(tcurh{index},'String',num2str(value,'%4d'));
  end

  %-------------------------- gotoSliderValues -------------------------
  function gotoSliderValues(time)
    pwmsig = zeros([manDim,1]);
    for ii = 1:manDim
      pwmsig(ii) = round(get(hslider{ii},'Value'));
    end
    arm.setServos(pwmsig, time, false);
  end

  %========================== Musec Callbacks ==========================
  %

  %---------------------------- minCallback ----------------------------
  function minCallback(source, eventdata, handles)
    text = get(source,'String');
    index = get(source, 'Value');
    value = str2double(text);
    if isnan(value)
      errordlg('You must enter a numeric value','Bad Input','modal');
      return;
    end
    set(hslider{index},'Min',value);
    [index, value]
  end
  
  %--------------------------- currCallback ----------------------------
  function currCallback(source, eventdata)
    pwmsig = str2num(get(source,'String'));
    index = get(source, 'UserData')
    set(hslider{index}, 'Value', pwmsig);
  
    if (doCts)
      gotoSliderValues(1);
    end
  end

  %---------------------------- maxCallback ----------------------------
  function maxCallback(source, eventdata, handles)
    text  = get(source,'String');
    index = get(source, 'Value');
    maxval = str2double(get(source,'String'));
    hslider{index}
    get(hslider{index});
    curv = round(get(hslider{index},'Value'));
    if isnan(maxval)
      errordlg('You must enter a numeric value','Bad Input','modal');
      return;
    end
    if (curv > maxval)
      setSliderValue(index, maxval);
    end
    set(hslider{index},'Max',maxval);
    [index, maxval]
  end

  %---------------------------- zeroCallback ---------------------------
  function zeroCallback(source, eventdata, handles)
    text = get(source,'String');
    index = get(source, 'Value');
    value = str2double(text);
    if isnan(value)
      errordlg('You must enter a numeric value','Bad Input','modal');
      return;
    end
  end

  %---------------------------- tminCallback ---------------------------
  function tminCallback(source, eventdata, handles)
    jointNum = get(source,'UserData');
    sigVal = round(get(hslider{jointNum},'Value'));
    set(eminh{jointNum}, 'String', num2str(sigVal));
    set(hslider{jointNum}, 'Min', sigVal, 'Value', sigVal);
  end

  %--------------------------- tzeroCallback ---------------------------
  function tzeroCallback(source, eventdata, handles)
    jointNum = get(source,'UserData');
    sigVal = round(get(hslider{jointNum},'Value'));
    set(zeroh{jointNum}, 'String', num2str(sigVal));
  end

  %---------------------------- tmaxCallback ---------------------------
  function tmaxCallback(source, eventdata, handles)
    jointNum = get(source,'UserData');
    sigVal = round(get(hslider{jointNum},'Value'));
    set(emaxh{jointNum}, 'String', num2str(sigVal));
    set(hslider{jointNum},'Max', sigVal, 'Value', sigVal);
  end


  %===================== Interface Button Callbacks ====================
  %

  %---------------------------- ctsCallback ----------------------------
  function ctsCallback(source, eventdata)
    doCts = (get(source, 'Value') == 1);
  end

  %---------------------------- gotoCallback ---------------------------
  function gotoCallback(source, eventdata)
    gotoSliderValues(3);
  end

  %---------------------------- saveCallback ---------------------------
  function saveCallback(source, eventdata)

  alphaLims = zeros(3, manDim);
  musecLims = zeros(3, manDim);
  alphaOrient = ones(1, manDim);
  for ii=1:manDim
    alphaLims(1,ii) = str2num(get(minv{ii},'String'));
    alphaLims(2,ii) = str2num(get(zerov{ii},'String'));
    alphaLims(3,ii) = str2num(get(maxv{ii},'String'));
    if (alphaLims(1,ii) > alphaLims(3, ii))
      alphaOrient(ii) = -1;
      alphaLims([1 3],ii) = alphaLims([3 1], ii);
    end

    musecLims(1,ii) = str2num(get(eminh{ii},'String'));
    musecLims(2,ii) = str2num(get(zeroh{ii},'String'));
    musecLims(3,ii) = str2num(get(emaxh{ii},'String'));
  end

  filename = get(hfile,'String')

  save(filename, 'musecLims', 'alphaLims', 'alphaOrient');

  end

  %---------------------------- doneCallback ---------------------------
  function doneCallback(source, eventdata)
    clf(fh);
    close(fh);
  end

end

%
%
%============================== lynx6_setup ==============================
