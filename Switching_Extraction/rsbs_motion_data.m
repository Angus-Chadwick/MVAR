if 0,
    figure;plot(AUXPOS.tax,AUXPOS.x,'r')
    hold on;plot(AUXPOS.tax,AUXPOS.y,'b')
end

if 0,
    dum1=diff(AUXPOS.x);
    dum1=[0;dum1];    
    dum1(abs(dum1)>5)=NaN;
    
    dum2=AUXPOS.xraw;
    [dum1,dum2];
    A=dum1;
    B=dum2;
    tax=AUXPOS.tax;
    
    figure;plot(tax,dum1,'r');
    hold on;plot(tax,dum2,'b--','LineWidth',2);    
end

if ismember(idinf{NID}.type,{'SEMIPLAY'}),
    xpos=diff(AUXPB.x);
    xpos=[0;xpos];
    xpos(abs(xpos)>5)=NaN;
else
    xpos=diff(AUXPOS.x);
    xpos=[0;xpos];
    xpos(abs(xpos)>5)=NaN;
end

if strcmp(idinf{NID}.type,'rec4playback')|strcmp(idinf{NID}.type,'playback_clamped')| strcmp(idinf{NID}.type,'playback_notclamped'),
    %    ok=not(AUXPB.x==999);
    %    pos = [AUXPB.tax(ok),AUXPB.x(ok),AUXPB.y(ok),NaN(size(AUXPB.x(ok)))];
    %else
    %    pos = [AUXPOS.tax,AUXPOS.x,AUXPOS.y,AUXPOS.trltyp];
end

timstrt = sel.start;
timstop = sel.lngth+sel.start;
clear MOT;
MOT.dat=cell([length(sel.trl),1]);
MOT.xpos=cell([length(sel.trl),1]);
MOT.tim=cell([length(sel.trl),1]);
if ismember(idinf{NID}.type,{'SEMIPLAY'}),
    MOT.tim2=cell([length(sel.trl),1]);
end
for T=1:length(sel.trl),
    rngix=find(AUXPOS.tax>sel.trl(T)+timstrt&AUXPOS.tax<sel.trl(T)+timstop);
    MOT.dat{T}=AUXPOS.xraw( rngix );    
    MOT.tim{T}=AUXPOS.tax( rngix )-sel.trl(T);
    if ismember(idinf{NID}.type,{'SEMIPLAY'}),
        rngix=find(AUXPB.tax>sel.trl(T)+timstrt&AUXPB.tax<sel.trl(T)+timstop);
        MOT.xpos{T}=xpos( rngix );
        MOT.tim2{T}=AUXPB.tax( rngix )-sel.trl(T);
    else
        MOT.xpos{T}=xpos( rngix );
    end
end

% interpolate data on common axis
%MOT.timaxi = timstrt:0.02:timstop;
MOT.timaxi = timstrt:0.005:timstop;
MOT.dati = NaN([length(sel.trl),length(MOT.timaxi)]); % initialize
MOT.rawdati  = NaN([length(sel.trl),length(MOT.timaxi)]); % initialize
MOT.rawxposi = NaN([length(sel.trl),length(MOT.timaxi)]); % initialize

% speed just before grating onset
%rawdati    = NaN([length(sel.trl),length(MOT.timaxi)]); % raw velocity
velstrt    = -round(0.5/0.02); % start of averaging window
velnsmp    = round(0.5/0.02); % length of averaging window
for T=1:length(sel.trl),
    tmp1=MOT.tim{T};
    tmp2=MOT.dat{T};
    if ~isempty(tmp2),
        tmpix=find(diff(tmp1)==0)+1; % otherwise problems with interpolation
        %tmp1(tmpix(1)-1:tmpix(1)+1);
        tmp1(tmpix)=[];tmp2(tmpix)=[];
        tmp3=interp1(tmp1,tmp2,MOT.timaxi);
        %rawdati(T,:)=tmp3;
        MOT.rawdati(T,:)=tmp3;
        
        % no smoothing
        if 1,
            ix1=find(~isnan(tmp3),1,'first'); ix2=find(~isnan(tmp3),1,'last'); % do not smooth with NaN's in vector, that is very slow!
            %tmp3(ix1:ix2)=smooth(tmp3(ix1:ix2),'rlowess',10);
            tmp3(ix1:ix2)=smooth(tmp3(ix1:ix2),'moving',10); % 20140530 faster
        end
        MOT.dati(T,:)=tmp3;
        
        if isfield(MOT,'tim2'),
            tmp1=MOT.tim2{T};
            tmp2=MOT.xpos{T};
            tmpix=find(diff(tmp1)==0)+1; % otherwise problems with interpolation
            tmp1(tmpix)=[];tmp2(tmpix)=[];
            tmp3=interp1(tmp1,tmp2,MOT.timaxi);
            MOT.rawxposi(T,:)=tmp3;
        else
            tmp2b=MOT.xpos{T};
            tmp2b(tmpix)=[];
            MOT.rawxposi(T,:)=interp1(tmp1,tmp2b,MOT.timaxi);
        end
    end % if ~isempty(tmp2),
end
% exclude first and last datapoints to avoid NaNs after interpolation
MOT.timaxi=MOT.timaxi(2:end-1);
MOT.dati=MOT.dati(:,2:end-1);
MOT.rawdati=MOT.rawdati(:,2:end-1);
MOT.rawxposi=MOT.rawxposi(:,2:end-1);

% figure;imagesc(MOT.timaxi,[1:size(rawdati,1)],rawdati); 
% figure;plot([1:size(rawdati,1)],MOT.prevel); 
%MOT.prevel = nanmean(rawdati(:,find(MOT.timaxi==0)+[velstrt:velstrt+velnsmp]),2);
MOT.prevel = nanmean(MOT.dati(:,find(MOT.timaxi==0)+[velstrt:velstrt+velnsmp]),2);

if 0,
    for T=1:length(sel.trl),
        figure;plot(MOT.tim{T},MOT.dat{T},'b-o');
        hold on;plot(MOT.timaxi,MOT.dati(T,:),'r');
        pause;close;
    end
end

if 0,
    
    % plot reward signal together with position
    figure;
    ptim=AUXPOS.tax(AUXPOS.x~=999); px=AUXPOS.x(AUXPOS.x~=999);
    px=(px-min(px))./(max(px)-min(px));
    plot(ptim,px,'bx'); % position info
    dat=(dat-min(dat))./(max(dat)-min(dat));
    hold on;plot(tim,dat,'r');    % lick on PC1, should always be shorter!
    set(gca,'XLim',[DIG.tax(1),DIG.tax(end)]);
    
    figure;hist(diff(AUXPOS.tax))
    
    % compute velocity
    ptim=AUXPOS.tax(AUXPOS.x~=999); px=AUXPOS.x(AUXPOS.x~=999);
    
    vlc = diff(px)./diff(ptim);
    figure;plot(vlc)
    
end

sacctims=NaN;

% find triggers within each trial
MOT.trltrglabel = {'G1on','G1off','G2on','G2off','Lck','Ron','Roff','sacc','eyecloseL','eyecloseR','eyeopenL','eyeopenR'};
clear tmp;
tmp{1}=gratingonnormal;
tmp{2}=gratingoffnormal;
tmp{3}=gratingon2;
tmp{4}=gratingoff2;
tmp{5}=licktimes;
tmp{6}=valveonind;
tmp{7}=valveoffind;
tmp{8}=sacctims*1000;
tmp{9}=blinktimes{1}*1000;
tmp{10}=blinktimes{2}*1000;
tmp{11}=blinktimesopen{1}*1000;
tmp{12}=blinktimesopen{2}*1000;
for i=1:length(tmp), tmp{i}=tmp{i}./1000; end % convert ms to seconds
for T=1:length(sel.trl),
    for i=1:length(tmp),
        rngix=find(tmp{i}>sel.trl(T)+timstrt&tmp{i}<sel.trl(T)+timstop);
        MOT.trg{T}{i}=tmp{i}(rngix)-sel.trl(T);
    end
end

if not(isnan(PT)),
    % extract eye position data
    SFeye=1./median(diff(PT));
    frmstrteye = round(sel.start*SFeye);
    frmstopeye = round((sel.lngth+sel.start)*SFeye);
    % init
    clear EYE;
    EYE.t=NaN([length(sel.trl),length([frmstrteye:frmstopeye])]);
    nam1={'x','y','r','avg','eh','ev','er','d'};
    for C=1:2,
        for i=1:length(nam1),
            EYE.(nam1{i}){C}=NaN([length(sel.trl),length([frmstrteye:frmstopeye])]);
        end
    end    
    taxeye = [frmstrteye:frmstopeye]./SFeye;
    for T=1:length(sel.trl), % T=2 T=3
        rngix=find(PT>sel.trl(T)+timstrt&PT<sel.trl(T)+timstop); % any times in this range?
        if ~isempty(rngix),
            [Y,I]=min(abs(taxeye-(PT(rngix(1))-sel.trl(T))));            
            %taxeye(I:I+length(taxeye)-I)            
            rngix1=I:I+length(taxeye)-I;
            rngix2=rngix(1):rngix(1)+length(taxeye)-I;
            
            rngix1=rngix1( rngix2<=length(PT) );
            rngix2=rngix2( rngix2<=length(PT) );
            
            EYE.t(T,rngix1)=PT( rngix2 )-sel.trl(T);
            
            for C=1:2,
                if ~isempty(PX{C}),
                    EYE.x{C}(T,rngix1)=PX{C}( rngix2 );
                    EYE.y{C}(T,rngix1)=PY{C}( rngix2 );
                    EYE.r{C}(T,rngix1)=PR{C}( rngix2 );
                    EYE.avg{C}(T,rngix1)=PA{C}( rngix2 );
                    EYE.eh{C}(T,rngix1)=Eh{C}( rngix2 );
                    EYE.ev{C}(T,rngix1)=Ev{C}( rngix2 );
                    EYE.er{C}(T,rngix1)=Er{C}( rngix2 );
                    if exist(sprintf('d%d',C)),
                        eval(sprintf('EYE.d{C}(T,rngix1)=d%d( rngix2 );',C)); % outliers (see rsbs_eye_data) but does not seem to be used
                    end
                end % if ~isempty(PX{C}),
            end % for C=1:2,       
        end % if ~isempty(rngix)   
    end % for T=1:length(sel.trl),
else
    EYE=[];
end % if not(isnan(PT)),