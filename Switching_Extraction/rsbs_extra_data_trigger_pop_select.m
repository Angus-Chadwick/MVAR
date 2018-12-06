function [DPOP]=rsbs_extra_data_trigger_pop_select(SIG,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);

%--------------------------------------------------------------
% if switching identify the valid trials

if strcmp(idinf{NID}.type,'SWITCH'),
    SIG.sel.valid=zeros(size(SIG.sel.cond));
    SIG.sel.valid( ~ismember(SIG.sel.cond,[1:6]) ) = 1; % only check for gratings and odour onsets
    TRLSEL = LN.TOTALLredux{RIXSES(1)}.TRL{RIXSES(2)};
    for i=1:6%length(TRLSEL.CN),
        if 0,
            %unique(TRLSEL.sel.cond( TRLSEL.CN{i} ))
            length(find(SIG.sel.cond==1))
            length(find(TRLSEL.sel.cond==1))
            length(find(SIG.sel.cond==3))
            length(find(TRLSEL.sel.cond==3))
            % '\M73_20141112_B2'
            
            [SIG.sel.trl( find(SIG.sel.cond==3) )-DELAY/1000,TRLSEL.sel.trl( find(TRLSEL.sel.cond==3) )]
            [SIG.sel.trl( find(SIG.sel.cond==3) ) - TRLSEL.sel.trl( find(TRLSEL.sel.cond==3) )]
        end
        
        if ismember(i,[3,4])
            DELAY = ODOURDELAY
        else
            DELAY = 0;
        end
        
        for j=1:length(TRLSEL.CN{i}), % for each trial
            %mtch = find( (SIG.sel.trl-DELAY/1000 ) == TRLSEL.sel.trl( TRLSEL.CN{i}(j) ));
            mtch = find( abs((SIG.sel.trl-DELAY/1000)-(TRLSEL.sel.trl( TRLSEL.CN{i}(j) )))<0.001 );
            if isempty(mtch),
                if ismember(idinf{NID}.id,{'M70_20141106','M71_20141101','M80_20141110','M87_20141110'}),
                    
                else
                    error('here');
                end
            else
                SIG.sel.valid( mtch ) = 1;
            end
        end
    end
else
    SIG.sel.valid=ones(size(SIG.sel.cond));
end

%--------------------------------------------------------------
% preprocess

if 1,
    
    ses = 1;
    
    for layer=1:length(SIG.DAT),
        badchn=find(isnan(nanmean(nanmean(SIG.DAT{layer},3),1))); % throw out chn with only NaNs
        okchn=find(~isnan(nanmean(nanmean(SIG.DAT{layer},3),1))); % throw out chn with only NaNs
        if not(isempty(badchn)),
            SIG.DAT{layer}          = SIG.DAT{layer}(:,okchn,:);
            SIG.CRSPR{layer}.conversion=SIG.CRSPR{layer}.conversion(okchn,:);
        end % if not(isempty(badchn)),
    end %  for layer=1:length(TOT{ses}.SIG.DAT),
    
    % restrict data points
    rok=find(SIG.tax>-2&SIG.tax<5);
    SIG.tax=SIG.tax(rok);
    for layer=1:length(SIG.DAT),
        SIG.DAT{layer}=SIG.DAT{layer}(:,:,rok);
    end
    
    if length(SIG.DAT)==1
        POOLLAYER=0;
        warning('not pooling layers');
    else
        POOLLAYER=1;
    end
    if POOLLAYER,
        POOLDAT=[]; POOLCV=[];
        for layer=1:length(SIG.DAT),
            POOLDAT = cat(2,POOLDAT,SIG.DAT{layer});
            POOLCV  = cat(1,POOLCV,[repmat(layer,[size(SIG.CRSPR{layer}.conversion,1),1]),SIG.CRSPR{layer}.conversion(:,2)]);
        end % for layer=1:length(SIG.DAT),
        TOTINF.POOLCV{ses} = POOLCV;
        SIG.DAT = {POOLDAT};
    else
        layer = 1;
        POOLCV = [repmat(layer,[size(SIG.CRSPR{layer}.conversion,1),1]),SIG.CRSPR{layer}.conversion(:,2)];
        TOTINF.POOLCV{ses} = POOLCV;
    end
    NL=length(SIG.DAT); % =1 after pooling, otherwise 4
    
    trgs=[1:length(SIG.BTRGLABEL)]';
    
    DPOP=cell([1,1]); % contains session info
    
    %------------------------------------------------------
    PSTHMOT=single(NaN([length(trgs),size(SIG.MOT.dati,2)]));
    PSTHMOTERR=single(NaN([length(trgs),size(SIG.MOT.dati,2)]));
    PSTHMOTn  =single(NaN([length(trgs),size(SIG.MOT.dati,2)]));
    for C=1:length(trgs),
        trls = find(SIG.sel.cond==trgs(C)&SIG.sel.valid); % selected trials
        PSTHMOT(C,:)=nanmean(SIG.MOT.dati(trls,:),1);
        PSTHMOTn(C,:)=sum(~isnan(SIG.MOT.dati(trls,:)),1);
        PSTHMOTERR(C,:,:) = squeeze( nanstd(SIG.MOT.dati(trls,:),[],1)./sqrt(sum( not(isnan(SIG.MOT.dati(trls,:))) ,1)) );
    end
    
    for layer=1:NL,
        PSTH   =single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
        PSTHERR=single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
        PSTHn  =single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
        for C=1:length(trgs),
            trls       = find(SIG.sel.cond==trgs(C)&SIG.sel.valid); % selected trials (exclude trials with only NaNs)
            PSTHn(C,:,:)   = sum(~isnan(SIG.DAT{layer}(trls,:,:)),1);
            PSTH(C,:,:)    = nanmean(SIG.DAT{layer}(trls,:,:),1);
            PSTHERR(C,:,:) = squeeze( nanstd(SIG.DAT{layer}(trls,:,:),[],1)./sqrt(sum( not(isnan(SIG.DAT{layer}(trls,:,:))) ,1)) );
        end % for C=1:length(trgs),
        
        DPOP{ses}.PL{layer}.PSTH    = PSTH;
        DPOP{ses}.PL{layer}.PSTHERR = PSTHERR;
        DPOP{ses}.PL{layer}.PSTHn   = PSTHn;
        
    end %  for layer=1:length(TOT{ses}.SIG.DAT),
    
    if strcmp(idinf{NID}.type,'SD'),
        C=strmatch('onsetgrey',SIG.BTRGLABEL);
        
        % for grey onset condition, set timepoint after grating onset to NaN
        
        %------------------------------------------------------
        G=[gratingonnormal;gratingon2]./1000;
        trls   = find(SIG.sel.cond==trgs(C)&SIG.sel.valid); % selected trials
        clear CURDATMOT CURDAT;
        CURDATMOT = SIG.MOT.dati(trls,:);
        for layer=1:NL,
            CURDAT{layer} = SIG.DAT{layer}(trls,:,:);
        end
        DIFF   = NaN(length(trls),1);
        if 0,
            figure;imagesc(SIG.MOT.timaxi,[1:size(CURDAT,1)],CURDAT);
        end
        for i=1:length(trls),
            tim_onset_grey    = SIG.sel.trl(trls(i));
            tim_onset_grating = G(find(G>tim_onset_grey,1,'first'));
            
            DIFF(i)=tim_onset_grating-tim_onset_grey;
        end % for i=1:length(trls),
        
        for i=1:length(DIFF),
            if DIFF(i)<SIG.MOT.timaxi(end)
                [~,Ix]=min(abs( SIG.MOT.timaxi-DIFF(i) ));
                CURDATMOT(i,Ix:end)=NaN;
            end
            for layer=1:NL,
                if DIFF(i)<SIG.tax(end)
                    [~,Ix]=min(abs( SIG.tax-DIFF(i) ));
                    CURDAT{layer}(i,:,Ix:end)=NaN;
                end
            end
        end
        if 0,
            figure;imagesc(SIG.MOT.timaxi,[1:size(CURDAT,1)],CURDAT);
            hold on;plot(DIFF,[1:size(CURDAT,1)],'w*');
        end
        
        C=size(PSTHMOT,1)+1;
        SIG.BTRGLABEL{C} = 'onsetgrey_untilgrating';
        PSTHMOT(C,:)=nanmean(CURDATMOT,1);
        PSTHMOTn(C,:)=sum(~isnan(CURDATMOT),1);
        PSTHMOTERR(C,:,:) = squeeze( nanstd(CURDATMOT,[],1)./sqrt(sum( not(isnan(CURDATMOT)) ,1)) );
        for layer=1:NL,
            DPOP{ses}.PL{layer}.PSTHn(C,:,:)   = sum(~isnan(CURDAT{layer}),1);
            DPOP{ses}.PL{layer}.PSTH(C,:,:)    = nanmean(CURDAT{layer},1);
            DPOP{ses}.PL{layer}.PSTHERR(C,:,:) = squeeze( nanstd(CURDAT{layer},[],1)./sqrt(sum( not(isnan(CURDAT{layer})) ,1)) );
        end %  for layer=1:length(TOT{ses}.SIG.DAT),
        DPOP{ses}.P.DIFF = DIFF;
    end
    
    DPOP{ses}.P.t          = SIG.tax;
    DPOP{ses}.P.tmot       = SIG.MOT.timaxi;
    DPOP{ses}.P.SFi        = SIG.SFi;
    DPOP{ses}.P.trgs       = trgs;
    DPOP{ses}.P.BTRGLABEL  = SIG.BTRGLABEL;
    DPOP{ses}.P.PSTHMOT    = PSTHMOT;
    DPOP{ses}.P.PSTHMOTERR = PSTHMOTERR;
    DPOP{ses}.P.PSTHMOTn   = PSTHMOTn;
    DPOP{ses}.P.sel        = SIG.sel;
    
    clear CN;
    for c=1:length(SIG.BTRGLABEL),
        CN{c}=find(SIG.sel.cond==c & SIG.sel.valid==1);
    end
    t=SIG.tax;
    DPOP{ses}.P.CN        = CN;
    %--------------------------------------------------------------
    if 1, % alternating trials for control
        
        clear CNrng;
        for ccnd=1:length(CN),
            CNrng{ccnd}{1}=CN{ccnd}( 1:2:length(CN{ccnd}) );
            CNrng{ccnd}{2}=CN{ccnd}( 2:2:length(CN{ccnd}) );
        end
        DPOP{ses}.P.CNrng        = CNrng;
        for CNrngix=1:2,
            for layer=1:NL,
                PSTH   =single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
                PSTHERR=single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
                PSTHn  =single(NaN([length(trgs),size(SIG.DAT{layer},2),size(SIG.DAT{layer},3)]));
                for C=1:length(trgs),
                    trls       =  CNrng{C}{CNrngix}; % selected trials (exclude trials with only NaNs)
                    PSTHn(C,:,:)   = sum(~isnan(SIG.DAT{layer}(trls,:,:)),1);
                    PSTH(C,:,:)    = nanmean(SIG.DAT{layer}(trls,:,:),1);
                    PSTHERR(C,:,:) = squeeze( nanstd(SIG.DAT{layer}(trls,:,:),[],1)./sqrt(sum( not(isnan(SIG.DAT{layer}(trls,:,:))) ,1)) );
                end % for C=1:length(trgs),
                
                DPOP{ses}.PL{layer}.CNrng = CNrng;
                DPOP{ses}.PL{layer}.ALT{CNrngix}.PSTH    = PSTH;
                DPOP{ses}.PL{layer}.ALT{CNrngix}.PSTHERR = PSTHERR;
                DPOP{ses}.PL{layer}.ALT{CNrngix}.PSTHn   = PSTHn;
                if 0,
                    figure;imagesc(squeeze(DPOP{ses}.PL{layer}.ALT{1}.PSTH(1,:,:)));
                    figure;imagesc(squeeze(DPOP{ses}.PL{layer}.ALT{2}.PSTH(1,:,:)));
                end
            end %  for layer=1:length(TOT{ses}.SIG.DAT),
        end % for CNrngix=1:2,
        
    end
    
    
end % if 1,

if 1, % slow
    
    %--------------------------------------------------------------
    % noise corr
    
    for layer=1:NL,
        
        t=SIG.tax;
        tMOT = SIG.MOT.timaxi;
        ixSTRT=find(t==0);
        [~,ixBSL] =min(abs( (t--1) )); 
        [~,ixSTOP]=min(abs( (t-1) )); 
        
        rngSTM = ixSTRT:ixSTOP; % t(rngNC)
        rngBSL = ixBSL:ixSTRT-1; % t(rngBSL)
        
        nr=size(SIG.DAT{layer},2);
        
        for ccnd=1:length(SIG.BTRGLABEL), % if only interested in V and A only do first 2
            
            % CAN PREPROCESS DATA HERE FURTHER IF NEEDED
            
            group1=SIG.DAT{layer}(CN{ccnd},:,:); % trl*ROI*time
            group2=SIG.MOT.dati(CN{ccnd},:); 
                       
            clear DCOR;
            DCOR.group1=group1;
            DCOR.group2=group2;
	
            DCOR.t      = t;
            DCOR.tMOT   = tMOT;
            DCOR.rngSTM  = rngSTM;
            DCOR.rngBSL = rngBSL;
            DPOP{ses}.PL{layer}.DCOR{ccnd}=DCOR;
        end % for ccnd=1:length(SIG.BTRGLABEL),
    end % for layer=1:length(TOT{ses}.SIG.DAT),
    
end % if 0,

DPOP{ses}.POOLCV = POOLCV;
