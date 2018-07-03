function varargout = PDPCaughtinAct14(varargin)
% PDPCAUGHTINACT14 v 1.2 Application M-file for PDPCaughtinAct.fig 
% W.M.Saj 2007
%    FIG = PDPCAUGHTINACT14 launch PDPCaughtinAct GUI.
%    PDPCAUGHTINACT14('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 25-Jan-2007 02:00:37

warning off;

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
        set(handles.text7,'String','');
	catch
        last_err=lasterr;
        set(handles.text7,'ForegroundColor',[1.0,0.0,0.0]);set(handles.text7,'String',last_err(:)');
        if ~isempty(findobj('Tag','Output for PDP Caught in Act '))
            close(findobj('Tag','Output for PDP Caught in Act '));
            set(handles.pushbutton2,'Enable','on');
            set(handles.pushbutton2,'String','Start');
        end
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
%| 'PDPCaughtInAct_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.PDPCaughtInAct, handles.slider2. This
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
[filename, pathname] = uigetfile('*.pdp', 'Pick an pdp-file');

if isequal(filename,0)|isequal(pathname,0)
    return
end
nazwapliku=[pathname filename];

if (get(handles.checkbox2,'Value')==1.0)
    pos=findstr(nazwapliku,'_t='); if isempty(pos) return; end;
    prefix=nazwapliku(1:pos-1); set(handles.edit4,'String',prefix);
    set(handles.edit5,'String',nazwapliku(pos+3:(findstr(nazwapliku,'.pdp')-1)));
    set(handles.edit6,'String','0');
    set(handles.edit7,'String',num2str(2*str2num(get(handles.edit5,'String'))));
end

load_plik(nazwapliku);
set_slider;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pole,nazwapola]=load_plik(nazwapliku)

handles=guihandles(gcbo);
set(handles.text7,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text7,'String','loading file ...');
[pole,nazwapola]=load_pdp2mod(nazwapliku);
set(handles.text3,'String',nazwapliku);
[xs,ys,zs]=size(pole);

set(handles.text11,'String',num2str(xs));
set(handles.text12,'String',num2str(ys));
set(handles.text13,'String',num2str(zs));
%-------------------------------------
function set_slider()
handles=guihandles(gcbo);
val_x=str2num(get(handles.text11,'String'));
val_y=str2num(get(handles.text12,'String'));
val_z=str2num(get(handles.text13,'String'));

if get(handles.radiobutton1,'Value')==1
    set(handles.slider1,'Max',val_x-1);
    set(handles.slider1,'Value',ceil((val_x-1)/2));
    step=min(ceil(max(val_y,val_z))/10,10);
    if val_x==1
        set(handles.slider1,'Visible','off')
    else
        set(handles.slider1,'Visible','on')
    end
else 
   if get(handles.radiobutton2,'Value')==1
        set(handles.slider1,'Max',val_y-1);
        set(handles.slider1,'Value',ceil((val_y-1)/2));
        step=min(ceil(max(val_x,val_z))/10,10);
        if val_y==1
            set(handles.slider1,'Visible','off')
        else
            set(handles.slider1,'Visible','on')
        end
   else
        set(handles.slider1,'Max',val_z-1);
        set(handles.slider1,'Value',ceil((val_z-1)/2));
        step=min(ceil(max(val_x,val_y))/10,10);
        if val_z==1
            set(handles.slider1,'Visible','off')
        else
            set(handles.slider1,'Visible','on')
        end
   end
end

set(handles.slider2,'Max',step);
set(handles.slider2,'Value',min(get(handles.slider2,'Value'),step));
if step==1
    set(handles.slider2,'Visible','off')
else
    set(handles.slider2,'Visible','on')
end
        
set(handles.text17,'String',num2str(ceil(get(handles.slider1,'Value'))));    
set(handles.text20,'String',num2str(ceil(get(handles.slider2,'Value'))));  
%------------------------------------------------------------------------
% load_pdp2.m ladowanie pliku pdp do tablicy matlaba
%% [tab,name,field_type]=load_pdp(filename)
function [pole,nazwa] = load_pdp2mod(name)

% otwieramy plik
[fid,message]=fopen(name,'r');
if message
    pole=[];nazwa=[];
    return;
end

% nazwa pola
try
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

    fclose(fid);

    if typ_pola>4
        tmp=pole(:,:,1:2:zs)+sqrt(-1)*pole(:,:,2:2:zs);
        clear pole;
        pole=tmp;
        clear tmp;
        pole=abs(pole);
    end
    
catch
    pole=[];nazwa=[];
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

nazwa_pliku=get(handles.text3,'String');
if (isequal(nazwa_pliku,'...') & (get(handles.checkbox2,'Value')==0.0))
    return
end

set(handles.pushbutton1,'Enable','off');
set(handles.pushbutton1,'String','Working...');

set(handles.pushbutton2,'Enable','off');
set(handles.pushbutton2,'String','Working...');

set(handles.checkbox2,'Enable','off');

xs=str2num(get(handles.text11,'String'));
ys=str2num(get(handles.text12,'String'));
zs=str2num(get(handles.text13,'String'));

if (get(handles.checkbox2,'Value')==1.0)
    prefix=get(handles.edit4,'String');
    t_min=str2num(get(handles.edit5,'String'));
    t_stp=str2num(get(handles.edit6,'String'));
    t_max=str2num(get(handles.edit7,'String'));
    t=t_min;
end

if ~isempty(findobj('Tag','Output for PDP Caught in Act '))
    h=findobj('Tag','Output for PDP Caught in Act ');
else
    h=figure;
end

set(h,'Name',['PDP Caught in Act view for file: ' nazwa_pliku]);
set(h,'Color',[0.1 0.3 1])
set(h,'Tag','Output for PDP Caught in Act ')
set(h,'CloseRequestFcn',';delete(gcf);')
set(h,'DoubleBuffer','on')
is_start=1;


while(1)
    
    dt=str2num(get(handles.edit1,'String'));
    if isempty(dt)
        dt=3;
    else
        dt=ceil(abs(dt(1)));
    end
    pause(0.01+dt);
    
    if (get(handles.checkbox2,'Value')==0.0)
        [pole,nazwapola]=load_plik(nazwa_pliku);
        if isempty(pole) 
            pause(dt);
            continue; 
        end
    else
        nazwa_pliku=[prefix '_t=' num2str(t) '.pdp'];
        if exist(nazwa_pliku,'file') 
            [pole,nazwapola]=load_plik(nazwa_pliku);
        else
            pause(dt);
            continue;
        end
        t=t+t_stp;
        if t>t_max
            h=findobj('Tag','Output for PDP Caught in Act ');
            if ~isempty(h)
                close(h);
            end
        end
    end
    

    
    cut_plane=floor(get(handles.slider1,'Value'))+1;
    step=floor(get(handles.slider2,'Value'));
    
    if get(handles.radiobutton1,'Value')==1    
        pole=squeeze(pole(cut_plane,1:step:ys,1:step:zs));
    else 
        if get(handles.radiobutton2,'Value')==1
            pole=squeeze(pole(1:step:xs,cut_plane,1:step:zs));
        else
            pole=squeeze(pole(1:step:xs,1:step:ys,cut_plane));             
        end
    end
    
    if isempty(findobj('Tag','Output for PDP Caught in Act '))
        set(handles.pushbutton1,'Enable','on');
        set(handles.pushbutton1,'String','Choose file');
        set(handles.pushbutton2,'Enable','on');
        set(handles.pushbutton2,'String','Start');
        set(handles.checkbox2,'Enable','On');
        break
    end
    
    roz=size(pole);
    set(handles.text7,'ForegroundColor',[0.0,1.0,0.0]);set(handles.text7,'String','Plotting...');
    
    figure(findobj('Tag','Output for PDP Caught in Act '));
    view_matrix=view;
    
    if (roz(1)==1 | roz(2)==1)
        
        plot(pole,'LineWidth',2.0);
        set(gca,'XTick',[]);
        set(gca,'YColor','w');

        if get(handles.checkbox1,'Value')==1
            
            y1=str2num(get(handles.edit2,'String'));
            y2=str2num(get(handles.edit3,'String'));
            if (isempty(y1) | isempty(y2) | (y1(1)>=y2(1)))
                y1=-1;y2=1;
            end
            ylim(gca,[y1(1) y2(1)]);              
        else
            ylim(gca,'auto');
            axis tight;
        end

        
    else 
        
        if get(handles.radiobutton4,'Value')==1
            view_matrix=view;
            surf(pole);shading interp;
            if is_start
                view([45 75]);
            else
                view(view_matrix);               
            end
        else
            if get(handles.radiobutton5,'Value')==1
                contourf(pole,20);
            else
                pcolor(pole);shading interp;
            end
        end
        
        figure(findobj('Tag','Output for PDP Caught in Act ')); %% jezeli sie zgubilo
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        set(gca,'ZColor','w');
       
        if get(handles.checkbox1,'Value')==1
          
            z1=str2num(get(handles.edit2,'String'));
            z2=str2num(get(handles.edit3,'String'));
            
            if (isempty(z1) | isempty(z2) | (z1(1)>=z2(1)))
                z1=-1;z2=1;    
            end
            zlim(gca,[z1(1) z2(1)]); caxis(gca,[z1(1) z2(1)]);             
        else
            zlim(gca,'auto');caxis('auto');
            axis tight;
        end   

    end
    
    title(nazwapola,'Color','w');

    
    if is_start==1
        is_start=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step=str2num(get(handles.edit1,'String'));
if ((isempty(step)) | (step(1)<=0))
  set(handles.edit1,'String','3')  
else
  set(handles.edit1,'String',num2str(ceil(step))); 
end

function varargout = edit2_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
step1=str2num(get(handles.edit2,'String'));
step2=str2num(get(handles.edit3,'String'));
if ((isempty(step1)) | (step1(1) >= step2(1)))
  set(handles.edit2,'String','-1')  
else
  set(handles.edit2,'String',num2str(step1(1))); 
end

function varargout = edit3_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

step1=str2num(get(handles.edit2,'String'));
step2=str2num(get(handles.edit3,'String'));

if ((isempty(step2)) | (step2(1) <= step1(1)))
  set(handles.edit3,'String','1')  
else
  set(handles.edit3,'String',num2str(step2(1))); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',1.0);
set(handles.radiobutton2,'Value',0.0);
set(handles.radiobutton3,'Value',0.0);

set(handles.text18,'String','x =');
set_slider;

function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',0.0);
set(handles.radiobutton2,'Value',1.0);
set(handles.radiobutton3,'Value',0.0);

set(handles.text18,'String','y =');
set_slider;

function varargout = radiobutton3_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton1,'Value',0.0);
set(handles.radiobutton2,'Value',0.0);
set(handles.radiobutton3,'Value',1.0);

set(handles.text18,'String','z =');
set_slider;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton4,'Value',1.0);
set(handles.radiobutton5,'Value',0.0);
set(handles.radiobutton6,'Value',0.0);

h=findobj('Tag','Output for PDP Caught in Act ');
if ~isempty(h)
        view([45 75]);        
end

function varargout = radiobutton5_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton4,'Value',0.0);
set(handles.radiobutton5,'Value',1.0);
set(handles.radiobutton6,'Value',0.0);

function varargout = radiobutton6_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);

set(handles.radiobutton4,'Value',0.0);
set(handles.radiobutton5,'Value',0.0);
set(handles.radiobutton6,'Value',1.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.text17,'String',num2str(ceil(get(handles.slider1,'Value'))));


% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)
handles=guihandles(gcbo);
if get(handles.checkbox1,'Value')==0.0
    set(handles.edit2,'Enable','off');
    set(handles.edit3,'Enable','off');
else
    set(handles.edit2,'Enable','on');
    set(handles.edit3,'Enable','on');
end


% --- Executes on slider movement.
function varargout = slider2_Callback(hObject, eventdata, handles, varargin)
handles=guihandles(gcbo);
set(handles.text20,'String',num2str(ceil(get(handles.slider2,'Value')))); 



function edit4_Callback(hObject, eventdata, handles)
handles=guihandles(gcbo);
set(handles.text3,'String','...');

nazwapliku=[get(handles.edit4,'String') '_t=' get(handles.edit5,'String') '.pdp'];
if exist(nazwapliku,'file') 
    load_plik(nazwapliku);
    set_slider;
end

function edit5_Callback(hObject, eventdata, handles)
handles=guihandles(gcbo);
w=str2num(get(handles.edit5,'String'));
if isempty(w)
    w=1000;
else
    w=ceil(abs(w(1)));
    if (w>str2num(get(handles.edit7,'String')))
        w=str2num(get(handles.edit7,'String'));
    end
end
set(handles.edit5,'String',num2str(w));
set(handles.text3,'String','...');

nazwapliku=[get(handles.edit4,'String') '_t=' get(handles.edit5,'String') '.pdp'];
if exist(nazwapliku,'file') 
    load_plik(nazwapliku);
    set_slider;
end

function edit6_Callback(hObject, eventdata, handles)
handles=guihandles(gcbo);
w=str2num(get(handles.edit6,'String'));
if isempty(w)
    w=1000;
else
    w=ceil(abs(w(1)));
end
set(handles.edit6,'String',num2str(w));

function edit7_Callback(hObject, eventdata, handles)
handles=guihandles(gcbo);
w=str2num(get(handles.edit7,'String'));
if isempty(w)
    w=10000;
else
    w=ceil(abs(w(1)));
    if (w<str2num(get(handles.edit5,'String')))
        w=str2num(get(handles.edit5,'String'));
    end
end
set(handles.edit7,'String',num2str(w));

function checkbox2_Callback(hObject, eventdata, handles)
handles=guihandles(gcbo);
if get(handles.checkbox2,'Value')==0.0
    set(handles.edit4,'Visible','Off');
    set(handles.edit5,'Visible','Off');
    set(handles.edit6,'Visible','Off');
    set(handles.edit7,'Visible','Off');
    
    set(handles.text22,'Visible','Off');
    set(handles.text23,'Visible','Off');
    set(handles.text24,'Visible','Off');
else
    nazwapliku=get(handles.text3,'String');  
    pos=findstr(nazwapliku,'_t=');
    if ~isempty(pos)
        prefix=nazwapliku(1:pos-1); set(handles.edit4,'String',prefix);
        set(handles.edit5,'String',nazwapliku(pos+3:(findstr(nazwapliku,'.pdp')-1)));
        set(handles.edit6,'String','0');
        set(handles.edit7,'String',num2str(2*str2num(get(handles.edit5,'String'))));
    end
    
    set(handles.edit4,'Visible','On');
    set(handles.edit5,'Visible','On');
    set(handles.edit6,'Visible','On');
    set(handles.edit7,'Visible','On');
    
    set(handles.text22,'Visible','On');
    set(handles.text23,'Visible','On');
    set(handles.text24,'Visible','On');
end


