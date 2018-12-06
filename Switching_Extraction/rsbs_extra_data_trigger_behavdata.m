function [sel,BTRGLABEL,MOT,ODOURDELAY,gratingonnormal,gratingon2]=rsbs_extra_data_trigger_behavdata(out,ADDINF,idinf,NID,EYE,EYEINF)

AUXPOS = ADDINF.AUXPOS;

ODOURDELAY = 500; % ms 

%--------------------------------------------------------------------------
% edit rsbs_VR_main_behavdata.m

RUNONSET = mouseVR2p_runonset(out,ADDINF.AUXPOS,idinf{NID});

% get runonset after reward onset
RUNONSETafter = mouseVR2p_runonset_after(out,ADDINF,idinf{NID},RUNONSET);


PARAM.gratingonnormal = ADDINF.gratingonnormal;
PARAM.gratingon2 = ADDINF.gratingon2;
STMONSET = mouseVR2p_greyonset(out,ADDINF.AUXPOS,idinf{NID},PARAM);

% for onsetgrey, find and remove timepoints where grating has appeared

if 0,%strcmp(idinf{NID}.type,'SD'),
    
    strmatch('onsetgrey',SIG.BTRGLABEL)
    A=BTRG{ strmatch('onsetgrey',BTRGLABEL) };
    B=[ BTRG{ strmatch('Von',BTRGLABEL) }; BTRG{ strmatch('Aon',BTRGLABEL) }];
    size(A),size(B)
    C=NaN(size(A)); % find next stim onset
    for i=1:length(A),
        C(i)=find(B>A(i),1,'first');
    end
    
    ( gratingonnormal - STMONSET.onsetgrey{1} )/1000
    ( STMONSET.onsetgrey{1} - STMONSET.onsetcircle{1} )/1000
    figure;plot(ADDINF.AUXPOS.tax,ADDINF.AUXPOS.y,'r');
    hold on;plot(ADDINF.AUXPOS.tax,ADDINF.AUXPOS.x,'b');
    dum=get(gca,'YLim');
    dum2=gratingonnormal./1000;    
    L=strcat('g1on',cellstr(num2str([1:length(gratingonnormal)]')))
    hold on;text(dum2,repmat(mean(dum),size(dum2)),L,'Rotation',90);
    dum2=STMONSET.onsetgrey{1}./1000;    
    L=strcat('grey1on',cellstr(num2str([1:length(gratingonnormal)]')))
    hold on;text(dum2,repmat(mean(dum),size(dum2)),L,'Rotation',90);
    dum2=STMONSET.onsetcircle{1}./1000;    
    L=strcat('circle1on',cellstr(num2str([1:length(gratingonnormal)]')))
    hold on;text(dum2,repmat(mean(dum),size(dum2)),L,'Rotation',90);
    
end
    

% make BTRG structure
if strcmp(idinf{NID}.type,'SD'),
    gratingonnormal  = ADDINF.gratingonnormal;
    gratingoffnormal = ADDINF.gratingoffnormal;
    gratingon2       = ADDINF.gratingon2;
    gratingoff2      = ADDINF.gratingoff2;
    licktimes        = ADDINF.licktimes;
    valveonind       = ADDINF.valveonind;
    valveoffind      = ADDINF.valveoffind;
  
    ons             = RUNONSET.offs; % offset of still = onset of running
    offs            = RUNONSET.ons;  % onset of still  = offset of running
    
    % only on and offsets in runup corridor, and stay at least 1s in runup corridor
    dum=[RUNONSET.offsrunupexit-RUNONSET.offs]/1000;
    ons             = RUNONSET.offs(dum>1); % offset of still = onset of running
    offs            = RUNONSET.ons(dum>1);  % onset of still  = offset of running
    
    onsetcircle     = STMONSET.onsetcircle;
    onsetgrey       = STMONSET.onsetgrey;
    
    % make triggers of following times
    cix=0; clear BTRG BTRGLABEL;
    cix=cix+1;BTRG{cix}=gratingonnormal;       BTRGLABEL{cix}='Von';
    cix=cix+1;BTRG{cix}=gratingon2;            BTRGLABEL{cix}='Aon';
    
    valvefirst=NaN(size(gratingonnormal));
    %clear valvefirst;vix=0;
    for i=1:length(gratingonnormal),
        tmpix=find(valveonind>gratingonnormal(i),1,'first');
        if ~isempty(tmpix),
            %vix=vix+1;
            valvefirst(i)=valveonind(tmpix);
        end
    end
    valvefirst = valvefirst(~isnan(valvefirst));
    cix=cix+1;BTRG{cix}=valvefirst;       BTRGLABEL{cix}='Ron';
    
    cix=cix+1;BTRG{cix}=ons;             BTRGLABEL{cix}='ons';
    cix=cix+1;BTRG{cix}=offs;            BTRGLABEL{cix}='offs';
    
    %     for cnd=1:2,
    %         cix=cix+1;BTRG{cix}=onsetcircle{cnd}(~isnan(onsetcircle{cnd})); BTRGLABEL{cix}=sprintf('onsetcircle%d',cnd);
    %         cix=cix+1;BTRG{cix}=onsetgrey{cnd}(~isnan(onsetgrey{cnd}));   BTRGLABEL{cix}=sprintf('onsetgrey%d',cnd);
    %     end
    
    % pool conditions
    cix=cix+1; BTRGLABEL{cix}='onsetcircle'; BTRG{cix}=[]; 
    for cnd=1:2,
        okix=find(~isnan(onsetcircle{cnd}));
        BTRG{cix}=[BTRG{cix};onsetcircle{cnd}(okix)];
    end
    cix=cix+1; BTRGLABEL{cix}='onsetgrey'; BTRG{cix}=[]; %BTRGgrating{cix}=[];
    for cnd=1:2,
        okix=find(~isnan(onsetgrey{cnd}));
        BTRG{cix}=[BTRG{cix};onsetgrey{cnd}(okix)];
        %BTRGgrating{cix}=[BTRGgrating{cix};[okix,repmat(cnd,size(okix))]];
    end
    
    cix=cix+1; BTRGLABEL{cix}='onsRwd'; BTRG{cix}=RUNONSETafter.onsRwd;
    
elseif strcmp(idinf{NID}.type,'SWITCH'),
    
    gratingonnormal = ADDINF.gratingonnormal;
    gratingon2      = ADDINF.gratingon2;
%     if ismember(idinf{NID}.id,{'M70_20141106','M71_20141101','M80_20141110','M87_20141110'}), % Oct 2015 removed because don't know why here
%         odour1ontimes   = [];
%         odour2ontimes   = [];
%         gratingolf1     = [];
%         gratingolf2     = [];
%     else
        odour1ontimes   = ADDINF.odour1ontimes;
        odour2ontimes   = ADDINF.odour2ontimes;
        gratingolf1     = ADDINF.gratingolf1;
        gratingolf2     = ADDINF.gratingolf2;
%     end
    
    gratingoffnormal = ADDINF.gratingoffnormal;
    gratingoff2      = ADDINF.gratingoff2;
    licktimes        = ADDINF.licktimes;    
    valveonind      = ADDINF.valveonind;
    valveoffind      = ADDINF.valveoffind;  
    
    ons             = RUNONSET.offs; % offset of still = onset of running
    offs            = RUNONSET.ons;  % onset of still  = offset of running    
    
    % only on and offsets in runup corridor, and stay at least 1s in runup corridor
    dum=[RUNONSET.offsrunupexit-RUNONSET.offs]/1000;
    ons             = RUNONSET.offs(dum>1); % offset of still = onset of running
    offs            = RUNONSET.ons(dum>1);  % onset of still  = offset of running
        
    onsetgrey       = STMONSET.onsetgrey;
    
    offsrunupexit = RUNONSET.offsrunupexit; % time it exits the runup corridor
    offsthreshold = RUNONSET.offsthreshold; % time it falls below the speed threshold again
        
    %[offs,ons,offsrunupexit,offsthreshold]
    %[offsthreshold-ons]./1000
        
     % make triggers of following times
        cix=0; clear BTRG BTRGLABEL;
        cix=cix+1;BTRG{cix}=gratingonnormal; BTRGLABEL{cix}='Von';
        cix=cix+1;BTRG{cix}=gratingon2;      BTRGLABEL{cix}='Aon';
        cix=cix+1;BTRG{cix}=odour1ontimes+ODOURDELAY;   BTRGLABEL{cix}='o1on';
        cix=cix+1;BTRG{cix}=odour2ontimes+ODOURDELAY;   BTRGLABEL{cix}='o2on';
        cix=cix+1;BTRG{cix}=gratingolf1;     BTRGLABEL{cix}='VonIrr';
        cix=cix+1;BTRG{cix}=gratingolf2;     BTRGLABEL{cix}='AonIrr';
            
    valvefirst=NaN(size(gratingonnormal));
    %clear valvefirst;vix=0;
    for i=1:length(gratingonnormal),
        tmpix=find(valveonind>gratingonnormal(i),1,'first');
        if ~isempty(tmpix),
            %vix=vix+1;
            valvefirst(i)=valveonind(tmpix);
        end
    end
    valvefirst = valvefirst(~isnan(valvefirst));
    cix=cix+1;BTRG{cix}=valvefirst;       BTRGLABEL{cix}='RonV';
    
    valvefirst=NaN(size(odour1ontimes));
    for i=1:length(odour1ontimes),
        tmpix=find(valveonind>odour1ontimes(i),1,'first');
        if ~isempty(tmpix),
            valvefirst(i)=valveonind(tmpix);
        end
    end
    valvefirst = valvefirst(~isnan(valvefirst));
    cix=cix+1;BTRG{cix}=valvefirst;       BTRGLABEL{cix}='RonO';
   
    cix=cix+1;BTRG{cix}=ons;             BTRGLABEL{cix}='ons';
    cix=cix+1;BTRG{cix}=offs;            BTRGLABEL{cix}='offs';
    
    cix=cix+1;BTRG{cix}=onsetgrey(~isnan(onsetgrey)); BTRGLABEL{cix}='onsetgrey';
   
    % start and stop of visual and olfactory block
    BLK      = ADDINF.BLK;
    STRTSTOP = ADDINF.STRTSTOP;
    ix = strmatch('onsetgrey',BTRGLABEL)
    ix1= BTRG{cix};
    ix1BLK = NaN(size(ix1));    
    for i=1:size(STRTSTOP), % block start and stop
        TM= ADDINF.AUXPOS.tax( STRTSTOP(i,:));
        inrng = find(ix1./1000>=TM(1)&ix1./1000<TM(2));
        if BLK(i)==1, % visual block
            ix1BLK(inrng)=1;
        else % olfactory block
            ix1BLK(inrng)=2;
        end
    end
    if any(isnan(ix1BLK)), error('unassigned grey onsets'); end
    
    cix=cix+1;BTRG{cix}=ix1(ix1BLK==1); BTRGLABEL{cix}='onsetgreyV';
    cix=cix+1;BTRG{cix}=ix1(ix1BLK==2); BTRGLABEL{cix}='onsetgreyO';
    
    if 0, % define whether grey onset is part of visual or olfactory block
        figure;hold on;plot(ADDINF.AUXPOS.tax, ADDINF.AUXPOS.y./100,'r') % NOTE time in .tax in seconds
        axl = get (gca,'YLim');
        
        ix1=ADDINF.gratingonnormal; % onsets of rewarded grating in visual block
        hold on;plot(ix1./1000,mean(axl),'bo','LineWidth',3); % trigger time in ms, convert to seconds
        ix1= ADDINF.gratingon2;% nonrewarded grating in visual block
        hold on;plot(ix1./1000,mean(axl),'ro','LineWidth',3);
        ix1= ADDINF.gratingolf1; % rewarded grating but now olfactory block
        hold on;plot(ix1./1000,mean(axl),'go','LineWidth',3);
        ix1= ADDINF.gratingolf2; % nonrewarded grating, olfactory block
        hold on;plot(ix1./1000,mean(axl),'ko','LineWidth',3);
%         
        cix = strmatch('onsetgrey',BTRGLABEL,'exact')
        ix1= BTRG{cix}; 
        hold on;plot(ix1./1000,mean(axl),'mx','LineWidth',3); % grey onsets
        
        cix = strmatch('onsetgreyV',BTRGLABEL)
        ix1= BTRG{cix}; 
        hold on;plot(ix1./1000,mean(axl)-.5,'b+','LineWidth',3); % grey onsets
        cix = strmatch('onsetgreyO',BTRGLABEL)
        ix1= BTRG{cix}; 
        hold on;plot(ix1./1000,mean(axl)-.5,'r+','LineWidth',3); % grey onsets
        
        for i=1:size(STRTSTOP), % block start and stop
            TM= ADDINF.AUXPOS.tax( STRTSTOP(i,:));
            if BLK(i)==1, % visual block
                hold on;h=patch(TM([1,2,2,1]),axl([1,1,2,2]),[0 0 1],'FaceAlpha',0.2,'EdgeColor','none');
            else % olfactory block
                hold on;patch(TM([1,2,2,1]),axl([1,1,2,2]),[1 0 0],'FaceAlpha',0.2,'EdgeColor','none');
            end
        end        
    end
    
    
elseif strcmp(idinf{NID}.type,'running_dark'),
    
    gratingonnormal=NaN; gratingoffnormal=NaN;
    gratingon2=NaN; gratingoff2=NaN;    
    ons             = RUNONSET.offs; % offset of still = onset of running
    offs            = RUNONSET.ons;  % onset of still  = offset of running
    valveonind    = ADDINF.valveonind;
    valveoffind   = ADDINF.valveoffind;
    licktimes     = ADDINF.licktimes;
    
    % make triggers of following times
    cix=0; clear BTRG BTRGLABEL;
    cix=cix+1;BTRG{cix}=ons;             BTRGLABEL{cix}='ons';
    cix=cix+1;BTRG{cix}=offs;            BTRGLABEL{cix}='offs';
    
    beep;
    cix=cix+1;BTRG{cix}=licktimes';   BTRGLABEL{cix}='Lck';
    cix=cix+1;BTRG{cix}=EYEINF.STIMl.*1000; BTRGLABEL{cix}='saccL';
    cix=cix+1;BTRG{cix}=EYEINF.STIMr.*1000; BTRGLABEL{cix}='saccR';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesclose{1}.*1000; BTRGLABEL{cix}='blinkcloseL';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesclose{2}.*1000; BTRGLABEL{cix}='blinkcloseR';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesopen{1}.*1000; BTRGLABEL{cix}='blinkopenL';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesopen{2}.*1000; BTRGLABEL{cix}='blinkopenR';
    
elseif strcmp(idinf{NID}.type,'odor_dark'),
        
    odour1ontimes = ADDINF.odour1ontimes;
    odour2ontimes = ADDINF.odour2ontimes;
    ons           = RUNONSET.offs; % offset of still = onset of running
    offs          = RUNONSET.ons;  % onset of still  = offset of running
    valveonind    = ADDINF.valveonind;    
    valveoffind   = ADDINF.valveoffind;
    licktimes     = ADDINF.licktimes;
    
    % make triggers of following times
    cix=0; clear BTRG BTRGLABEL;
    cix=cix+1;BTRG{cix}=odour1ontimes+ODOURDELAY;   BTRGLABEL{cix}='o1on';
    cix=cix+1;BTRG{cix}=odour2ontimes+ODOURDELAY;   BTRGLABEL{cix}='o2on';
    cix=cix+1;BTRG{cix}=ons;             BTRGLABEL{cix}='ons';
    cix=cix+1;BTRG{cix}=offs;            BTRGLABEL{cix}='offs';
    cix=cix+1;BTRG{cix}=valveonind';           BTRGLABEL{cix}='Rwd';
    
    cix=cix+1;BTRG{cix}=licktimes';   BTRGLABEL{cix}='Lck';
    cix=cix+1;BTRG{cix}=EYEINF.STIMl.*1000; BTRGLABEL{cix}='saccL';
    cix=cix+1;BTRG{cix}=EYEINF.STIMr.*1000; BTRGLABEL{cix}='saccR';    
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesclose{1}.*1000; BTRGLABEL{cix}='blinkcloseL';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesclose{2}.*1000; BTRGLABEL{cix}='blinkcloseR';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesopen{1}.*1000; BTRGLABEL{cix}='blinkopenL';
    cix=cix+1;BTRG{cix}=EYE.TOT.EYEDAT.blinktimesopen{2}.*1000; BTRGLABEL{cix}='blinkopenR';
    
    % min(EYEINF.STIMl),max(EYEINF.STIMl)
    % min(EYE.TOT.EYEDAT.blinktimesclose{1}),max(EYE.TOT.EYEDAT.blinktimesclose{1})
    % min(odour1ontimes),max(odour1ontimes)
    %edit rsbs_ROI_CRSP_rsp_switch_pop_eyedat.m
    
    lickfirst=[];
    gratingonnormal=NaN; gratingoffnormal=NaN;
    gratingon2=NaN; gratingoff2=NaN;
    
end

%-----------------------------------------------------------
trglabel{1} = 'word';
trglabel{2} = 'stim_on';
trglabel{3} = 'targ_on';
tWRD = strmatch('word',trglabel);
tSTM = strmatch('stim_on',trglabel);
tTRG = strmatch('targ_on',trglabel);

NTRG=0;
for i=1:length(BTRG),
    NTRG=NTRG+length(BTRG{i});
end

tol=1e-9;
clear ADD;
ADD.id = (tSTM-1)*ones(NTRG,1); % stim onset bit
ADD.val = zeros(NTRG,1);
ADD.tax = [];
for i=1:length(BTRG),
    ADD.tax = cat(1,ADD.tax,BTRG{i}./1000);
end
% word bit
for i=1:length(BTRG),
    ADD.id  = cat(1,ADD.id,(tWRD-1)*ones(size(BTRG{i})));
    ADD.val = cat(1,ADD.val,i*ones(size(BTRG{i})));
    ADD.tax = cat(1,ADD.tax,(BTRG{i}./1000)-tol); % word comes before stim on bit
end

DIG.id  = [];
DIG.val = [];
DIG.tax = [];
% append
DIG.id  = cat(1,DIG.id,ADD.id);
DIG.val  = cat(1,DIG.val,ADD.val);
DIG.tax = cat(2,DIG.tax,ADD.tax');

% sort bit times
[DY,DI] = sort(DIG.tax);
DIG.id  = DIG.id(DI);
DIG.val = DIG.val(DI);
DIG.tax = DIG.tax(DI);

trgix = find(DIG.id==tSTM-1); % should be equal to NTRG

clear TRL;
TRL.tim = NaN([length(trgix),3]);
TRL.wrd = NaN([length(trgix),1]);
for i=1:length(trgix)%i=length(trgix)
    TRL.tim(i,1) = DIG.tax(trgix(i));
    ix = intersect(find(DIG.id==tTRG-1),find(DIG.tax>DIG.tax(trgix(i))));
    if ~isempty(ix), TRL.tim(i,2) = DIG.tax(ix(1)); end
    
    ix = intersect(find(DIG.id==tWRD-1),find(DIG.tax<DIG.tax(trgix(i)))); % word comes before stim on bit
    if ~isempty(ix),
        TRL.tim(i,3) = DIG.tax(ix(end));
        TRL.wrd(i,1) = DIG.val(ix(end));
    end
end
% cut out data around trigger

% cut out data around trigger
sel.lngth =  10;
sel.start = -5;

sel.trl   = TRL.tim(:,1); % stim on in this case
if isempty(TRL.tim)
    sel.cond=[];
else
    sel.cond  = TRL.wrd;
end

% rsbs_ROI_CRSP_bhv_eye
%   rsbs_motion_eyes
do_eyes=0;
clear blinktimes*;
clear PT PX PY PR PA Eh Ev Er d1 d2;
if do_eyes
    if ~isempty(EYE), % already made by rsbs_ROI_CRSP_bhv_eye
        EYEDAT = EYE.TOT.EYEDAT;
        
        % this ensures the data is extract in rsbs_motion_data
        PT=EYEDAT.PT;
        PX=EYEDAT.PX;
        PY=EYEDAT.PY;
        PR=EYEDAT.PR;
        PA=EYEDAT.PA;
        Eh=EYEDAT.Eh;
        Ev=EYEDAT.Ev;
        Er=EYEDAT.Er;
        blinktimes     = EYEDAT.blinktimesclose;
        blinktimesopen = EYEDAT.blinktimesopen;        
    else
        % JP 20160102: should not need to be called 
        % rsbs_eye_data; % resscn_eye_data20130410;
        
    end
else
    PT=NaN;
    sacctims=NaN;
    blinktimes{1}=NaN;
    blinktimes{2}=NaN;
    blinktimesopen{1}=NaN;
    blinktimesopen{2}=NaN;
end

rsbs_motion_data;
