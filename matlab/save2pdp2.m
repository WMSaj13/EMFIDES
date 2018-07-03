% SAVE2PDP2(TAB,NAZWA_PLIKU,CIAG,PRECYZJA)
% zapisuje tablice Matlaba do pliku pdp
% TAB - nazwa tablicy (NDIMS(TAB)<=3)
% NAZWA_PLIKU - nazwa pliku w formacie pdp do ktorego zapisujemy tablice
% CIAG - ciag zankow zapisywany do pliku pdp
% PRECYZJA - precyzja zapisu , akceptowane wartosci : 'float',
% 'double','short','int','complex float','complex double'
%
% zobacz takze LOAD_PDP2
function save2pdp(tablica,nazwa_pliku,nazwa_pola,precyzja)

if ndims(tablica)>3
    disp('too many dimensions');
    return;
end

%typ pola
pole_typ={'float','double','short','int'};

is_real_type=0;is_complex_type=0; typ=-1;

for i=1:4
    if isequal(pole_typ(i),precyzja)    
        typ=i-1;
        break
    end
end

if isequal('complex float',precyzja)
    typ=4;precyzja='float';
end

if isequal('complex double',precyzja)
    typ=5;precyzja='double';
end

if (typ==-1)
    disp('unknown field type');
end

[fid,message]=fopen(nazwa_pliku,'w');
if message
    error(message);
end

xs=size(tablica,1);
ys=size(tablica,2);
zs=size(tablica,3);

% nazwa pola
fwrite(fid,length(nazwa_pola),'int');
fwrite(fid,nazwa_pola,'uchar');

% rozmiary
fwrite(fid,xs,'int');
fwrite(fid,ys,'int');
fwrite(fid,zs,'int');

%% typ pola
fwrite(fid,typ,'int');

try
    if (typ>=4)
        tab=zeros(xs,ys,2*zs);
        tab(:,:,1:2:2*zs)=real(tablica(:,:,:));
        tab(:,:,2:2:2*zs)=imag(tablica(:,:,:));
        tablica=tab;
        clear tab;
    end

    % zapis pola            
    for i=1:xs
        tmp=squeeze(tablica(i,:,:))';
        fwrite(fid,tmp,precyzja);
    end
    fclose(fid);
catch
    fclose(fid);
    delete(nazwa_pliku);
    error('file not saved')
end