% LOAD_PDP2.m ladowanie pliku pdp do tablicy Matlaba v 2.0
% [TAB,CIAG,TYP_POLA]=LOAD_PDP2(NAZWA_PLIKU)
% NAZWA_PLIKU - nazwa wczytywanego pliku
% TAB - nazwa tablicy do ktorej wczytujemy dane
% CIAG - nazwa ciagu do ktorego wczywtujemy ciag znakow 
% przechowywyanych w pliku pdp
% TYP_POLA - format wczytywanego pola , wartosci 'float',
% 'double','short','int','complex float','complex double'
%
% TAB=LOAD_PDP2(NAZWA_PLIKU) i [TAB,CIAG]=LOAD_PDP2(NAZWA_PLIKU)
% to rowiez poprawne wywolania funkcji
%
% zobacz takze SAVE2PDP2

function varargout = load_pdp2(name)

if nargout>3 
    error('to many arguments');
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

pole_typ={'float','double','short','int','float','double'};
% wczytujemy pole            

if typ_pola>4
    zs=2*zs;
end

pole=zeros(xs,ys,zs);
tmp=zeros(zs,ys);

try
    for i=1:xs
        [tmp,count]=fread(fid,[zs,ys],char(pole_typ(typ_pola)));
        pole(i,:,:)=tmp';
    end
    fclose(fid);
catch
    fclose(fid);
    error('file not readed');
    varargout=[];
    return;
end

if typ_pola>4
    tmp=pole(:,:,1:2:zs)+sqrt(-1)*pole(:,:,2:2:zs);
    clear pole;
    pole=tmp;
    clear tmp;
    pole_typ=['complex ' char(pole_typ(typ_pola))];
else
    pole_typ=pole_typ(typ_pola);
end

varargout(1)={pole};

if nargout==2 
    varargout(2)={nazwa};
end

if nargout==3 
    varargout(2)={nazwa};
    varargout(3)={pole_typ};
end

