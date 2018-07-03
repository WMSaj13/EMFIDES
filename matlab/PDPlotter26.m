function varargout = PDPlotter26(varargin)
% PDPLOTTER26 Application M-file for PDPlotter.fig v 2.4 W.M.Saj April 2006
%    FIG = PDPLOTTER26 launch PDPlotter GUI.
%    PDPLOTTER26('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 24-Jan-2007 14:59:59
warning off;
global pole;

if nargin == 0  % LAUNCH GUI
 
	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    handles=guihandles(gcbo);
	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        set(handles.text18,'String','');
	catch
        last_err=lasterr;
        set(handles.text18,'ForegroundColor',[1.0,0.0,0.0]);set(handles.text18,'String',last_err(:)');
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'PDPPlotter_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.PDPPlotter, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

[filename, pathname] = uigetfile('*.pdp', 'Pick an pdp-file');

if isequal(filename,0)|isequal(pathname,0)
    return
end
nazwapliku=[pathname filename];
load_plik(nazwapliku);

% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);
nazwapliku=get(handles.text1,'String');
if isequal(nazwapliku,'...')
    return
end
load_plik(nazwapliku)

% --------------------------------------------------------------------
function load_plik(nazwapliku)

global pole;

handles=guihandles(gcbo);
set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','loading file ...');
[pole,nazwapola,pola_typ]=load_pdp2(nazwapliku);

set(handles.text1,'String',nazwapliku);
set(handles.text13,'String',[ nazwapola '(' pola_typ ')']);

if nazwapliku(length(nazwapliku)-9:length(nazwapliku))=='_Freal.pdp'
    pole_temp=load_pdp2([nazwapliku(1:length(nazwapliku)-10) '_Fimag.pdp']);    
    pole=pole+sqrt(-1)*pole_temp;
    clear pole_temp;
    set(handles.text31,'String',[nazwapliku(1:length(nazwapliku)-10) '_Fimag.pdp']);
else
    if nazwapliku(length(nazwapliku)-9:length(nazwapliku))=='_Fimag.pdp'
        pole_temp=load_pdp2([nazwapliku(1:length(nazwapliku)-10) '_Freal.pdp']);    
        pole=sqrt(-1)*pole+pole_temp;
        clear pole_temp;
        set(handles.text31,'String',[nazwapliku(1:length(nazwapliku)-10) '_Freal.pdp']);
    else
        set(handles.text31,'String','...');
    end
end

[xs,ys,zs]=size(pole);

set(handles.text5,'String',num2str(xs));
set(handles.text6,'String',num2str(ys));
set(handles.text7,'String',num2str(zs));

maxzas;
plane=checkplane3(get(handles.edit2,'String'),1);
set(handles.edit1,'String',plane);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)

global pole;

if isempty(pole)
    return
end

handles=guihandles(gcbo);

[filename, pathname] = uiputfile('*.mat', 'Save as');

if isequal(filename,0)|isequal(pathname,0)
    return
end
nazwapola=get(handles.text13,'String');

set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','saving file ...');

[tablica,zm1,zm2,zm1_nazwa,zm2_nazwa]=getpole;

zm1=zm1(1,:);
zm2=zm2(:,1);

save([pathname filename],'tablica','nazwapola','zm1','zm2','zm1_nazwa','zm2_nazwa'); 

% --------------------------------------------------------------------
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',1.0);
set(handles.radiobutton2,'Value',0.0);
set(handles.radiobutton3,'Value',0.0);

set(handles.text8,'String','Z =');
set(handles.text9,'String','X =');
set(handles.text10,'String','Y =');

set(handles.edit1,'String',check_bound_z(get(handles.edit1,'String'),1));

set(handles.edit2,'String',check_bound_x(get(handles.edit2,'String'),1));
set(handles.edit3,'String',check_bound_x(get(handles.edit3,'String'),0));
set(handles.edit4,'String',check_bound_x(get(handles.edit4,'String'),1));

set(handles.edit5,'String',check_bound_y(get(handles.edit5,'String'),1));
set(handles.edit6,'String',check_bound_y(get(handles.edit6,'String'),0));
set(handles.edit7,'String',check_bound_y(get(handles.edit7,'String'),1));

set_sliders_value;

xs=str2num(get(handles.text5,'String'));
ys=str2num(get(handles.text6,'String'));
zs=str2num(get(handles.text7,'String'));

set_sliders_max(xs,ys,zs);

% --------------------------------------------------------------------
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',0.0);
set(handles.radiobutton2,'Value',1.0);
set(handles.radiobutton3,'Value',0.0);

set(handles.text8,'String','X =');
set(handles.text9,'String','Y =');
set(handles.text10,'String','Z =');

set(handles.edit1,'String',check_bound_x(get(handles.edit1,'String'),1));

set(handles.edit2,'String',check_bound_y(get(handles.edit2,'String'),1));
set(handles.edit3,'String',check_bound_y(get(handles.edit3,'String'),0));
set(handles.edit4,'String',check_bound_y(get(handles.edit4,'String'),1));

set(handles.edit5,'String',check_bound_z(get(handles.edit5,'String'),1));
set(handles.edit6,'String',check_bound_z(get(handles.edit6,'String'),0));
set(handles.edit7,'String',check_bound_z(get(handles.edit7,'String'),1));

set_sliders_value;

xs=str2num(get(handles.text5,'String'));
ys=str2num(get(handles.text6,'String'));
zs=str2num(get(handles.text7,'String'));

set_sliders_max(ys,zs,xs);

% --------------------------------------------------------------------
function varargout = radiobutton3_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',0.0);
set(handles.radiobutton2,'Value',0.0);
set(handles.radiobutton3,'Value',1.0);

set(handles.text8,'String','Y =');
set(handles.text9,'String','Z =');
set(handles.text10,'String','X =');

set(handles.edit1,'String',check_bound_y(get(handles.edit1,'String'),1));

set(handles.edit2,'String',check_bound_z(get(handles.edit2,'String'),1));
set(handles.edit3,'String',check_bound_z(get(handles.edit3,'String'),0));
set(handles.edit4,'String',check_bound_z(get(handles.edit4,'String'),1));

set(handles.edit5,'String',check_bound_x(get(handles.edit5,'String'),1));
set(handles.edit6,'String',check_bound_x(get(handles.edit6,'String'),0));
set(handles.edit7,'String',check_bound_x(get(handles.edit7,'String'),1));

set_sliders_value;

xs=str2num(get(handles.text5,'String'));
ys=str2num(get(handles.text6,'String'));
zs=str2num(get(handles.text7,'String'));

set_sliders_max(zs,xs,ys);

%%%%
function set_sliders_value()

handles=guihandles(gcbo);
set(handles.slider7,'Value',str2num(get(handles.edit1,'String')));

set(handles.slider1,'Value',str2num(get(handles.edit2,'String')));
set(handles.slider2,'Value',str2num(get(handles.edit3,'String'))-1);
set(handles.slider3,'Value',str2num(get(handles.edit4,'String')));

set(handles.slider4,'Value',str2num(get(handles.edit5,'String')));
set(handles.slider5,'Value',str2num(get(handles.edit6,'String'))-1);
set(handles.slider6,'Value',str2num(get(handles.edit7,'String')));

function set_sliders_max(is,js,ks)

handles=guihandles(gcbo);

if ks>1
    set(handles.slider7,'Max',ks-1);
    set(handles.slider7,'Visible','on');  
else
    set(handles.slider7,'Visible','off');      
end

if is>1
    set(handles.slider1,'Max',is-1);
    set(handles.slider2,'Max',is-1);
    set(handles.slider3,'Max',is-1);
    
    set(handles.slider1,'Visible','on');
    set(handles.slider2,'Visible','on');
    set(handles.slider3,'Visible','on');
    
else 
    set(handles.slider1,'Visible','off');
    set(handles.slider2,'Visible','off');
    set(handles.slider3,'Visible','off');
end

if js>1
    
    set(handles.slider4,'Max',js-1);
    set(handles.slider5,'Max',js-1);
    set(handles.slider6,'Max',js-1);
    
    set(handles.slider4,'Visible','on');
    set(handles.slider5,'Visible','on');
    set(handles.slider6,'Visible','on');
    
else 
    set(handles.slider4,'Visible','off');
    set(handles.slider5,'Visible','off');
    set(handles.slider6,'Visible','off');
end

% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);

if isequal(get(handles.text1,'String'),'...')
    return
end
set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','plotting...');

[tablica,zm1,zm2,string1,string2]=getpole;
figure;
if check1D(tablica,zm1,zm2,string1,string2)==1;
    return
end

tablica=forma_tab(tablica);
minv=min(min(tablica));
maxv=max(max(tablica));

if (get(handles.checkbox7,'Value')==0)
    v1=str2num(get(handles.edit22,'String'));
    v2=str2num(get(handles.edit23,'String'));
    v1=v1*(maxv-minv)+minv;
    v2=v2*(maxv-minv)+minv;
else
    v1=minv;
    v2=maxv;
end

surf(zm1,zm2,tablica);
afterplot3D(string1,string2,0,v1,v2);

% --------------------------------------------------------------------
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

if isequal(get(handles.text1,'String'),'...')
    return
end
set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','plotting...');

[tablica,zm1,zm2,string1,string2]=getpole;
figure;
if check1D(tablica,zm1,zm2,string1,string2)==1;
    return
end

tablica=forma_tab(tablica);
minv=min(min(tablica));
maxv=max(max(tablica));

if (get(handles.checkbox7,'Value')==0)
    v1=str2num(get(handles.edit22,'String'));
    v2=str2num(get(handles.edit23,'String'));
    v1=v1*(maxv-minv)+minv;
    v2=v2*(maxv-minv)+minv;
else
    v1=minv;
    v2=maxv;
end

pcolor(zm1,zm2,tablica);
afterplot3D(string1,string2,0,v1,v2);

% --------------------------------------------------------------------
function varargout = pushbutton6_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

if isequal(get(handles.text1,'String'),'...')
    return
end
set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','plotting...');

[tablica,zm1,zm2,string1,string2]=getpole;
figure;
if check1D(tablica,zm1,zm2,string1,string2)==1;
    return
end
tablica=forma_tab(tablica);

minv=min(min(tablica));
maxv=max(max(tablica));

if (get(handles.checkbox7,'Value')==0)
    v1=str2num(get(handles.edit22,'String'));
    v2=str2num(get(handles.edit23,'String'));
    v1=v1*(maxv-minv)+minv;
    v2=v2*(maxv-minv)+minv;
else
    v1=minv;
    v2=maxv;
end

v=[v1:((v2-v1)/str2num(get(handles.edit14,'String'))):v2];

[cs,h]=contour(zm1,zm2,tablica,v);
if isequal(get(handles.checkbox8,'Value'),1)
    set(h,'LineColor','w')
end

if isequal(get(handles.checkbox5,'Value'),1)
    h=clabel(cs,h);
    if isequal(get(handles.checkbox8,'Value'),1)
        set(h,'Color','w')
    end
end
afterplot3D(string1,string2,0,v1,v2);

% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane3(get(handles.edit1,'String'),0);
set(handles.edit1,'String',plane);
set(handles.slider7,'Value',str2num(get(handles.edit1,'String')));

% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane1(get(handles.edit2,'String'),1);
set(handles.edit2,'String',plane);
set(handles.slider1,'Value',str2num(get(handles.edit2,'String')));
% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane1(get(handles.edit3,'String'),0);
set(handles.edit3,'String',plane);
set(handles.slider2,'Value',str2num(get(handles.edit3,'String'))-1);
% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane1(get(handles.edit4,'String'),1);
set(handles.edit4,'String',plane);
set(handles.slider3,'Value',str2num(get(handles.edit4,'String')));
% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane2(get(handles.edit5,'String'),1);
set(handles.edit5,'String',plane);
set(handles.slider4,'Value',str2num(get(handles.edit5,'String')));
% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane2(get(handles.edit6,'String'),0);
set(handles.edit6,'String',plane);
set(handles.slider5,'Value',str2num(get(handles.edit6,'String'))-1);
% --------------------------------------------------------------------
function varargout = edit7_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
plane=checkplane2(get(handles.edit7,'String'),1);
set(handles.edit7,'String',plane);
set(handles.slider6,'Value',str2num(get(handles.edit7,'String')));
%% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)
%--------------------------------------------------------------------
function varargout = edit9_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit9,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit9,'String','1') 
else
    set(handles.edit9,'String',num2str(step(1))); 
end
% --------------------------------------------------------------------
function varargout = edit10_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit11_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit11,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit11,'String','1') 
else
  set(handles.edit11,'String',num2str(step(1))); 
end
% --------------------------------------------------------------------
function varargout = edit12_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit13_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit13,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit13,'String','1')  
else
  set(handles.edit13,'String',num2str(step(1))); 
end

% --------------------------------------------------------------------
function varargout = edit14_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit14,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit14,'String','20')  
else
  set(handles.edit14,'String',num2str(ceil(step(1)))); 
end
% --------------------------------------------------------------------
function varargout = edit15_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit16_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit16,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit16,'String','12')  
else
  set(handles.edit16,'String',num2str(ceil(step(1)))); 
end
% --------------------------------------------------------------------
function varargout = edit17_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit17,'String'));
if ((isempty(step)))
  set(handles.edit17,'String','0')  
else
  set(handles.edit17,'String',num2str(step(1))); 
end
% --------------------------------------------------------------------
function varargout = edit18_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit18,'String'));
if ((isempty(step)))
  set(handles.edit18,'String','0')  
else
  set(handles.edit18,'String',num2str(step(1))); 
end
% --------------------------------------------------------------------
function varargout = edit19_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit19,'String'));
if ((isempty(step)))
  set(handles.edit19,'String','0')  
else
  set(handles.edit19,'String',num2str(step(1))); 
end
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.radiobutton4,'Value',1.0);
set(handles.radiobutton5,'Value',0.0);
set(handles.radiobutton10,'Value',0.0);

% --------------------------------------------------------------------
function varargout = radiobutton5_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.radiobutton4,'Value',0.0);
set(handles.radiobutton5,'Value',1.0);
set(handles.radiobutton10,'Value',0.0);
% --------------------------------------------------------------------
function varargout = radiobutton10_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.radiobutton4,'Value',0.0);
set(handles.radiobutton5,'Value',0.0);
set(handles.radiobutton10,'Value',1.0);
% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function checkbox5_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function varargout = pushbutton8_Callback(h, eventdata, handles, varargin)
maxzas

% --------------------------------------------------------------------
function varargout = pushbutton10_Callback(h, eventdata, handles, varargin)

handles=guihandles(gcbo);
if isequal(get(handles.text1,'String'),'...')
    return
end
set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','plotting...');

[tablica,zm1,zm2,string1,string2]=getpole;
figure;
if check1D(tablica,zm1,zm2,string1,string2)==1;
    return
end
tablica=forma_tab(tablica);

minv=min(min(tablica));
maxv=max(max(tablica));

if (get(handles.checkbox7,'Value')==0)
    v1=str2num(get(handles.edit22,'String'));
    v2=str2num(get(handles.edit23,'String'));
    v1=v1*(maxv-minv)+minv;
    v2=v2*(maxv-minv)+minv;
else
    v1=minv;
    v2=maxv;
end

v=[v1:((v2-v1)/str2num(get(handles.edit14,'String'))):v2];
[cs,h]=contourf(zm1,zm2,tablica,v);

if isequal(get(handles.checkbox8,'Value'),1)
    set(h,'LineColor','w')
end
if isequal(get(handles.checkbox5,'Value'),1)
    h=clabel(cs,h);
    if isequal(get(handles.checkbox8,'Value'),1)
        set(h,'Color','w')
    end
end
afterplot3D(string1,string2,1,v1,v2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tablica,zm1,zm2,string1,string2]=getpole()

global pole;
handles=guihandles(gcbo);
plane_num=str2num(get(handles.edit1,'String'));

start1=str2num(get(handles.edit2,'String'));
krok1=str2num(get(handles.edit3,'String'));
end1=str2num(get(handles.edit4,'String'));

start2=str2num(get(handles.edit5,'String'));
krok2=str2num(get(handles.edit6,'String'));
end2=str2num(get(handles.edit7,'String'));

if get(handles.radiobutton1,'Value')==1.0
    tablica=squeeze(pole((start1+1):krok1:(end1+1),(start2+1):krok2:(end2+1),plane_num+1)).';
    start1=(start1+str2num(get(handles.edit17,'String')))*str2num(get(handles.edit9,'String'));
    krok1=krok1*str2num(get(handles.edit9,'String'));
    end1=(end1+str2num(get(handles.edit17,'String')))*str2num(get(handles.edit9,'String'));
    start2=(start2+str2num(get(handles.edit18,'String')))*str2num(get(handles.edit11,'String'));
    krok2=krok2*str2num(get(handles.edit11,'String'));
    end2=(end2+str2num(get(handles.edit18,'String')))*str2num(get(handles.edit11,'String'));   
    string1=get(handles.edit8,'String');
    string2=get(handles.edit10,'String');
end

if get(handles.radiobutton2,'Value')==1.0
    tablica=squeeze(pole(plane_num+1,(start1+1):krok1:(end1+1),(start2+1):krok2:(end2+1))).';
    start1=(start1+str2num(get(handles.edit18,'String')))*str2num(get(handles.edit11,'String'));
    krok1=krok1*str2num(get(handles.edit11,'String'));
    end1=(end1+str2num(get(handles.edit18,'String')))*str2num(get(handles.edit11,'String'));
    start2=(start2+str2num(get(handles.edit19,'String')))*str2num(get(handles.edit13,'String'));
    krok2=krok2*str2num(get(handles.edit13,'String'));
    end2=(end2+str2num(get(handles.edit19,'String')))*str2num(get(handles.edit13,'String'));       
    string1=get(handles.edit10,'String');
    string2=get(handles.edit12,'String');
end

if get(handles.radiobutton3,'Value')==1.0
    tablica=squeeze(pole((start2+1):krok2:(end2+1),plane_num+1,(start1+1):krok1:(end1+1)));
    start1=(start1+str2num(get(handles.edit19,'String')))*str2num(get(handles.edit13,'String'));
    krok1=krok1*str2num(get(handles.edit13,'String'));
    end1=(end1+str2num(get(handles.edit19,'String')))*str2num(get(handles.edit13,'String'));
    start2=(start2+str2num(get(handles.edit17,'String')))*str2num(get(handles.edit9,'String'));
    krok2=krok2*str2num(get(handles.edit9,'String'));
    end2=(end2+str2num(get(handles.edit17,'String')))*str2num(get(handles.edit9,'String'));   
    string1=get(handles.edit12,'String');
    string2=get(handles.edit8,'String');
end

[zm1,zm2]=meshgrid([start1:krok1:end1],[start2:krok2:end2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxzas()

global pole;

if isempty(pole)
    return
end

handles=guihandles(gcbo);

[xs,ys,zs]=size(pole);

set(handles.edit2,'String','0');
set(handles.edit5,'String','0');

if get(handles.radiobutton1,'Value')==1
    set(handles.edit4,'String',num2str(xs-1));
    set(handles.edit7,'String',num2str(ys-1));
    set_sliders_value;set_sliders_max(xs,ys,zs);
end

if get(handles.radiobutton2,'Value')==1
    set(handles.edit4,'String',num2str(ys-1));
    set(handles.edit7,'String',num2str(zs-1));
    set_sliders_value;set_sliders_max(ys,zs,xs);
end

if get(handles.radiobutton3,'Value')==1
    set(handles.edit4,'String',num2str(zs-1));
    set(handles.edit7,'String',num2str(xs-1));
    set_sliders_value;set_sliders_max(zs,xs,ys);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plane=check_bound_x(plane,shift)
handles=guihandles(gcbo);
plane=str2num(plane);
if isempty(plane)
    plane=1-shift;
end
plane=num2str(ceil(real(max(1-shift,min(plane(1),str2num(get(handles.text5,'String'))-shift)))));
%%%%%%%%%%%
function plane=check_bound_y(plane,shift)
handles=guihandles(gcbo);
plane=str2num(plane);
if isempty(plane)
    plane=1-shift;
end
plane=num2str(ceil(real(max(1-shift,min(plane(1),str2num(get(handles.text6,'String'))-shift)))));
%%%%%%%%%%%%
function plane=check_bound_z(plane,shift)
handles=guihandles(gcbo);
plane=str2num(plane);
if isempty(plane)
    plane=1-shift;
end
plane=num2str(ceil(real(max(1-shift,min(plane(1),str2num(get(handles.text7,'String'))-shift)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yesno=check1D(tablica,zm1,zm2,string1,string2)
handles=guihandles(gcbo);

[r1,r2]=size(zm1);
tablica=forma_tab(tablica);

if (r1==1 | r2==1)
    
    if (r1==1)
        if get(handles.checkbox3,'Value')==1.0 
            semilogy(squeeze(zm1),squeeze(tablica));
        else
            plot(squeeze(zm1),squeeze(tablica));
        end
    else
        if get(handles.checkbox3,'Value')==1.0 
            semilogy(squeeze(zm2),squeeze(tablica));
        else
            plot(squeeze(zm2),squeeze(tablica));
        end
    end
    
    axis tight;
    set(gcf,'Name','PDP Plotter Output Window');
    set(gca,'FontSize',str2num(get(handles.edit16,'String')));
    
    if (r1==1)
        xlabel(string1);
    else
        xlabel(string2);
    end
    
    ylabel(get(handles.edit15,'String'));
    
    if get(handles.checkbox4,'Value')==1.0
        grid on;
    else
        grid off;
    end
    yesno=1;
    return;
end

yesno=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function afterplot3D(string1,string2,isshad,v1,v2)
handles=guihandles(gcbo);

caxis([v1,v2]);
set(gcf,'Name','PDP Plotter Output Window');
if get(handles.checkbox1,'Value')==1.0 & isshad==0
    shading interp;    
end
if get(handles.radiobutton4,'Value')==1.0
    colormap jet;
else
    if get(handles.radiobutton5,'Value')==1.0
        colormap gray
    else
        temp=gray(64);temp=temp(64:-1:1,:);colormap(temp);
    end
end

axis tight;
set(gca,'FontSize',str2num(get(handles.edit16,'String')));
xlabel(string1);
ylabel(string2);
zlabel(get(handles.edit15,'String'));

if get(handles.checkbox2,'Value')==1.0
    h=colorbar;
    set(h,'FontSize',str2num(get(handles.edit16,'String')));
end

if get(handles.checkbox4,'Value')==1.0
    grid on;
else
    grid off;
end

title(get(handles.edit15,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plane=checkplane1(plane,shift)
handles=guihandles(gcbo);

if get(handles.radiobutton1,'Value')==1
    plane=check_bound_x(plane,shift);
end
if get(handles.radiobutton2,'Value')==1
    plane=check_bound_y(plane,shift);
end
if get(handles.radiobutton3,'Value')==1
    plane=check_bound_z(plane,shift);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plane=checkplane2(plane,shift)
handles=guihandles(gcbo);

if get(handles.radiobutton1,'Value')==1
    plane=check_bound_y(plane,shift);
end
if get(handles.radiobutton2,'Value')==1
    plane=check_bound_z(plane,shift);
end
if get(handles.radiobutton3,'Value')==1
    plane=check_bound_x(plane,shift);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plane=checkplane3(plane,shift)
handles=guihandles(gcbo);

if get(handles.radiobutton1,'Value')==1
    plane=check_bound_z(plane,shift);
end
if get(handles.radiobutton2,'Value')==1
    plane=check_bound_x(plane,shift);
end
if get(handles.radiobutton3,'Value')==1
    plane=check_bound_y(plane,shift);
end

% --------------------------------------------------------------------
function varargout = pushbutton11_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

if get(handles.radiobutton1,'Value')==1
    plane=get(handles.text7,'String');
end
if get(handles.radiobutton2,'Value')==1
    plane=get(handles.text5,'String');
end
if get(handles.radiobutton3,'Value')==1
    plane=get(handles.text6,'String');
end

plane=num2str(floor(str2num(plane)/2));
set(handles.edit1,'String',plane);
set(handles.slider7,'Value',str2num(get(handles.edit1,'String')));


% --------------------------------------------------------------------
function varargout = radiobutton6_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton6,'Value',1.0);
set(handles.radiobutton7,'Value',0.0);
set(handles.radiobutton8,'Value',0.0);
set(handles.radiobutton9,'Value',0.0);

% --------------------------------------------------------------------
function varargout = radiobutton7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton7.
handles=guihandles(gcbo);

set(handles.radiobutton6,'Value',0.0);
set(handles.radiobutton7,'Value',1.0);
set(handles.radiobutton8,'Value',0.0);
set(handles.radiobutton9,'Value',0.0);


% --------------------------------------------------------------------
function varargout = radiobutton8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton8.
handles=guihandles(gcbo);

set(handles.radiobutton6,'Value',0.0);
set(handles.radiobutton7,'Value',0.0);
set(handles.radiobutton8,'Value',1.0);
set(handles.radiobutton9,'Value',0.0);


% --------------------------------------------------------------------
function varargout = radiobutton9_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton9.
handles=guihandles(gcbo);

set(handles.radiobutton6,'Value',0.0);
set(handles.radiobutton7,'Value',0.0);
set(handles.radiobutton8,'Value',0.0);
set(handles.radiobutton9,'Value',1.0);

%------------------------------------------------------------------------
% load_pdp2.m ladowanie pliku pdp do tablicy matlaba
%% [tab,name,field_type]=load_pdp(filename)
function varargout = load_pdp2(name)
if nargout>3 
    disp('to many arguments');
    return;
end
    
% otwieramy plik
[fid,message]=fopen(name,'r');
if message
    error(message);
end

% nazwa pola
length_of_name=fread(fid,1,'int');
[nazwa,count]=fread(fid,length_of_name,'uchar');

nazwa=char(nazwa');

% wczytujemy rozmiary
xs=fread(fid,1,'int');
ys=fread(fid,1,'int');
zs=fread(fid,1,'int');

%wczytujemy typ pola
typ_pola=fread(fid,1,'int');
typ_pola=typ_pola+1;

pola_typ={'float','double','short','int','float','double'};
pola_typ=char(pola_typ(typ_pola));

% wczytujemy pole 
if typ_pola>4
    zs=2*zs;
end

pole=zeros(xs,ys,zs);
tmp=zeros(zs,ys);


for i=1:xs
    [tmp,count]=fread(fid,[zs,ys],pola_typ);
    pole(i,:,:)=tmp';
end

if typ_pola>4
    tmp=pole(:,:,1:2:zs)+sqrt(-1)*pole(:,:,2:2:zs);
    clear pole;
    pole=tmp;
    clear tmp;
    pola_typ=['complex ' pola_typ];
end

varargout(1)={pole};

if nargout==2 
    varargout(2)={nazwa};
end

if nargout==3 
    varargout(2)={nazwa};
    varargout(3)={pola_typ};
end

fclose(fid);

%%%------------------------------------------------
function outtablica = forma_tab(tablica)

handles=guihandles(gcbo);
if get(handles.radiobutton6,'Value')==1
    outtablica=real(tablica);
end
if get(handles.radiobutton7,'Value')==1
    outtablica=imag(tablica);
end
if get(handles.radiobutton8,'Value')==1
    outtablica=abs(tablica);
end
if get(handles.radiobutton9,'Value')==1
     outtablica=angle(tablica);
end


% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit2,'String', num2str(ceil(get(handles.slider1,'Value'))));
if get(handles.slider1,'Value')>get(handles.slider3,'Value')
    set(handles.slider3,'Value',get(handles.slider1,'Value'));
    set(handles.edit4,'String', num2str(ceil(get(handles.slider3,'Value'))));    
end


% --------------------------------------------------------------------
function varargout = slider2_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit3,'String', num2str(ceil(get(handles.slider2,'Value'))+1));

% --------------------------------------------------------------------
function varargout = slider3_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit4,'String', num2str(ceil(get(handles.slider3,'Value'))));
if get(handles.slider1,'Value')>get(handles.slider3,'Value')
    set(handles.slider1,'Value',get(handles.slider3,'Value'));
    set(handles.edit2,'String', num2str(ceil(get(handles.slider1,'Value'))));    
end

% --------------------------------------------------------------------
function varargout = slider4_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit5,'String', num2str(ceil(get(handles.slider4,'Value'))));
if get(handles.slider4,'Value')>get(handles.slider6,'Value')
    set(handles.slider6,'Value',get(handles.slider4,'Value'));
    set(handles.edit7,'String', num2str(ceil(get(handles.slider6,'Value'))));    
end
% --------------------------------------------------------------------
function varargout = slider5_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit6,'String', num2str(ceil(get(handles.slider5,'Value'))+1));
% --------------------------------------------------------------------
function varargout = slider6_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit7,'String', num2str(ceil(get(handles.slider6,'Value'))));
if get(handles.slider4,'Value')>get(handles.slider6,'Value')
    set(handles.slider4,'Value',get(handles.slider6,'Value'));
    set(handles.edit5,'String', num2str(ceil(get(handles.slider4,'Value'))));    
end
% --------------------------------------------------------------------
function varargout = slider7_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit1,'String', num2str(ceil(get(handles.slider7,'Value'))));


% --- Executes on slider movement.
function varargout = slider8_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit22,'String', num2str((get(handles.slider8,'Value')),3));
if get(handles.slider8,'Value')>get(handles.slider9,'Value')
    set(handles.slider9,'Value',get(handles.slider8,'Value'));
    set(handles.edit23,'String', num2str((get(handles.slider9,'Value')),3));    
end


% --- Executes on slider movement.
function varargout = slider9_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.edit23,'String', num2str((get(handles.slider9,'Value')),3));
if get(handles.slider8,'Value')>get(handles.slider9,'Value')
    set(handles.slider8,'Value',get(handles.slider9,'Value'));
    set(handles.edit22,'String', num2str((get(handles.slider8,'Value')),3));    
end

function varargout = edit22_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
v=str2num(get(handles.edit22,'String'));
if isempty(v)
    v=0.0
end
set(handles.edit22,'String',num2str(v,3));
set(handles.slider8,'Value',v);
slider9_Callback(hObject, eventdata, handles, varargin)

function varargout = edit23_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
v=str2num(get(handles.edit23,'String'));
if isempty(v)
    v=1.0
end
set(handles.edit23,'String',num2str(v,3));
set(handles.slider9,'Value',v);
slider8_Callback(hObject, eventdata, handles, varargin)

function varargout = checkbox7_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
if (get(handles.checkbox7,'Value')==0)
    set(handles.slider8,'Enable','On');
    set(handles.slider9,'Enable','On');
    set(handles.edit22,'Enable','On');
    set(handles.edit23,'Enable','On');
    edit22_Callback(hObject, eventdata, handles, varargin)
    edit23_Callback(hObject, eventdata, handles, varargin)
else
    set(handles.slider8,'Enable','Off');
    set(handles.slider9,'Enable','Off');
    set(handles.edit22,'Enable','Off');
    set(handles.edit23,'Enable','Off');    
end


% --- Executes on button press in checkbox8.
function  varargout = checkbox8_Callback(hObject, eventdata, handles, varargin)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)

global pole;

handles=guihandles(gcbo);

variables = evalin('base','whos');if isempty(variables) msgbox('There are no variables in workspace','Message','error');return; end
nam={variables.name};
[s,v] = listdlg('Name','Load variable from workspace','ListSize',[320,320],'PromptString','Select a variable:',...
                'SelectionMode','single',...
                'ListString',nam);
if (v==0) return; end
pole_tmp=evalin('base',nam{s});
if (isempty(pole_tmp) | iscell(pole_tmp)|~isnumeric(pole_tmp)| ndims(pole_tmp)>3) return; end
pole=pole_tmp;

set(handles.text18,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text18,'String','loading variable ...');

set(handles.text1,'String','N/A (loaded from workspace)');
if isreal(pole)
    set(handles.text13,'String',[nam{s} ' (real array)']);
else
    set(handles.text13,'String',[nam{s} ' (complex array)']);
end
set(handles.text31,'String','...');

[xs,ys,zs]=size(pole);

set(handles.text5,'String',num2str(xs));
set(handles.text6,'String',num2str(ys));
set(handles.text7,'String',num2str(zs));

maxzas;
plane=checkplane3(get(handles.edit2,'String'),1);
set(handles.edit1,'String',plane);



