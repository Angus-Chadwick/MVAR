function [DPOP]=rsbs_extra_data_trigger_pop(SIG,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2)

%-----------------------------------------------------------
% edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m


% edit resscn_dFoverF
try
    fprintf('Trying to open matlabpool...\n')
    matlabpool open
catch
    n_workers =  matlabpool('size');
    fprintf('Already open, you are connected to %d workers \n', n_workers)
end
n_workers =  matlabpool('size');

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
    
    % % restrict data points
    % rok=find(SIG.tax>-2&SIG.tax<2);
    % SIG.tax=SIG.tax(rok);
    % for layer=1:4,
    %     SIG.DAT{layer}=SIG.DAT{layer}(:,:,rok);
    % end
    
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
    
    
    clear CN;
    for c=1:length(SIG.BTRGLABEL),
        CN{c}=find(SIG.sel.cond==c & SIG.sel.valid==1);
    end
    t=SIG.tax;
    
end % if 1,

%--------------------------------------------------------------
% TwinPP, TwinNC, DCORrng

if 1, % slow
    
    % define instead in nsamples
    %Twin.timwinext=2:length(t)-1;
    %Twin.twextsmp=[Twin.timwinext'-1,Twin.timwinext'+1];
    Twin.twextsmp=[];
    % add some more windows
    ADDWIN=[-1 0; 0 1; -0.25 0; 0 0.25; 0.25 0.5; 0 0.4; 0 0.5; 0.1 0.25; 0.1 0.35; 0.15 0.25; 0.15 0.3; 0.15 0.35; 0.2 0.3; 0.2 0.35; 0.05 0.3; 0.3 1.0; 0+0.6 1+0.6];
    ADDWIN=[0 1];
    for i=1:size(ADDWIN),
        Twin.twextsmp=[Twin.twextsmp; find(t>=ADDWIN(i,1),1,'first'),find(t<=ADDWIN(i,2),1,'last')];
    end
    Twin.twext=t(Twin.twextsmp);
    
    for layer=1:NL,
        
        clear TwinPP;
        TwinPP.twext = Twin.twext;
        TwinPP.twextsmp = Twin.twextsmp;
        
        TwinPP.twref=[-0.5 0];
        
        nr=size(SIG.DAT{layer},2); nw=size(TwinPP.twext,1); ncnd=length(SIG.BTRGLABEL);
        TwinPP.m1        = single(NaN([nw,nr,ncnd]));
        TwinPP.m2        = single(NaN([nw,nr,ncnd]));
        TwinPP.v1        = single(NaN([nw,nr,ncnd]));
        TwinPP.v2        = single(NaN([nw,nr,ncnd]));
        TwinPP.n1        = single(NaN([nw,nr,ncnd]));
        TwinPP.n2        = single(NaN([nw,nr,ncnd]));
        TwinPP.S1c       = single(NaN([nw,nr,nr,ncnd]));
        TwinPP.S2c       = single(NaN([nw,nr,nr,ncnd]));
        TwinPP.DPi       = single(NaN([nw,nr,ncnd]));
        TwinPP.p         = single(NaN([nw,nr,ncnd]));
        TwinPP.md        = single(NaN([nw,nr,ncnd]));
        for w=1:nw,%w=43,TwinPP.twext(w,:)
            for ccnd=1:length(SIG.BTRGLABEL),
                rng1=TwinPP.twextsmp(w,1):TwinPP.twextsmp(w,2);
                rng2=find(t>=TwinPP.twref(1)&t<=TwinPP.twref(2)); % pre-stimulus baseline
                
                if strcmp(SIG.BTRGLABEL{ccnd},'onsetgrey_untilgrating'),
                    group1=nanmean(CURDAT{layer}(:,:,rng1),3)';
                    group2=nanmean(CURDAT{layer}(:,:,rng2),3)';
                else
                    group1=nanmean(SIG.DAT{layer}(CN{ccnd},:,rng1),3)';
                    group2=nanmean(SIG.DAT{layer}(CN{ccnd},:,rng2),3)';
                end
                if 0,
                    % edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
                    % edit rsbsI_sd_pop_trigger.m
                    % edit rsbs_sd_pop_trigger_prepro.m
                    % rsbsI_sd_pop_trigger_makeDPOP
                    rng=find(DPOP{ses}.P.t>=0&DPOP{ses}.P.t<=1);
                    A=nanmean(DPOP{ses}.PL{layer}.PSTH(ccnd,:,rng),3)
                    B=nanmean(group1,2)';
                    if 0,
                        figure;plot(A,B);
                        A-B
                    end
                    if ~isequalwithequalnans_tolerance(A,B,0.001),
                        error('here');
                    end
                end
                if size(group1,2)>10&size(group2,2)>10,
                    % compute covariance matrix, in matlab each row is an observation, and each column a variable
                    S1=nancov(group1');   S2=nancov(group2');
                    try,
                        S1c=corrcov(S1);      S2c=corrcov(S2);
                    catch
                        S1(isnan(S1))=0;
                        S2(isnan(S2))=0;
                        S1c=corrcov(S1);
                        S2c=corrcov(S2);
                    end
                    
                    m1=nanmean(group1,2); m2=nanmean(group2,2);
                    n1=sum(~isnan(group1),2); n2=sum(~isnan(group2),2);
                    
                    TwinPP.m1(w,:,ccnd)      = m1;
                    TwinPP.m2(w,:,ccnd)      = m2;
                    TwinPP.n1(w,:,ccnd)      = n1;
                    TwinPP.n2(w,:,ccnd)      = n2;
                    TwinPP.S1c(w,:,:,ccnd)   = S1c;
                    TwinPP.S2c(w,:,:,ccnd)   = S2c;
                    
                    % compute individual d-primes
                    v1=nanvar(group1,[],2); v2=nanvar(group2,[],2);
                    poolstd=sqrt( (v1.*(n1-1)+v2.*(n2-1))./(n1+n2-2) ); % pooledstd
                    DPi=(m1-m2)./poolstd;
                    TwinPP.DPi(w,:,ccnd) = DPi; % sum(DPi.^2) % is equal to d_sq_shuf
                    
                    TwinPP.v1(w,:,ccnd)      = v1;
                    TwinPP.v2(w,:,ccnd)      = v2;
                    
                    % auroc
                    for i=1:size(group1,1),
                        A=group1(i,:); B=group2(i,:); A=A(~isnan(A)); B=B(~isnan(B));
                        %TwinPP.auc(w,i,ccnd)=calc_ROC(A,B); % Twin.auc(w,:)
                        TwinPP.p(w,i,ccnd)  =ranksum(A,B);
                        TwinPP.md(w,i,ccnd) =nanmedian(group1(i,:))-nanmedian(group2(i,:));
                    end
                end % if size(group1,2)>10&size(group1,2)>10,
            end % for ccnd=1:2,
        end % for w=1:nw,
        DPOP{ses}.PL{layer}.TwinPP=TwinPP;
    end % for layer=1:length(TOT{ses}.SIG.DAT),
    
    %--------------------------------------------------------------
    % tot corr, signal corr and noise corr
    
    for layer=1:NL,
        
        t=SIG.tax;
        ixSTRT=find(t==0);
        [~,ixBSL] =min(abs( (t--0.5) ));
        [~,ixBSL2] =min(abs( (t--1) )); % t(ixBSL2)
        [~,ixSTOP]=min(abs( (t-1.1) )); % t(ixSTOP)
        [~,ixSTOP2]=min(abs( (t-1) )); % t(ixSTOP2)
        
        rngSC=ixBSL:ixSTOP; % t(rngSC) include some baseline
        rngNC =ixSTRT:ixSTOP; % t(rngNC)
        rngNC2=ixSTRT:ixSTOP2; % t(rngNC)
        nr=size(SIG.DAT{layer},2);
        
        rngBSL=ixBSL2:ixSTRT-1 % t(rngBSL)
        
        for ccnd=1:length(SIG.BTRGLABEL),
            
            group1=SIG.DAT{layer}(CN{ccnd},:,:); % trl*ROI*time
            group1AVG=nanmean(group1,1); % 1*ROI*time
            group1S=group1-repmat(group1AVG,[size(group1,1),1,1]); % trl*ROI*time
            group1AVG=squeeze(group1AVG); % ROI*time
            
            group1AVGoverTIME=nanmean(group1(:,:,rngNC2),3); % 1*ROI*time
            group1AVGoverTIMEbsl=nanmean(group1(:,:,[rngBSL,rngBSL(end)+1]),3); % 1*ROI*time t([rngBSL,rngBSL(end)+1])
            
            CORTOT=single(NaN([nr,nr]));
            CORNSE1=single(NaN([nr,nr])); % NN method
            CORNSE2=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
            CORNSE3=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
            
            % same but for baseline
            CORNSE4=single(NaN([nr,nr])); % NN method
            CORNSE5=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
            CORNSE6=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
            
            if isfield(SIG,'PARTIAL'),
                
                PACORTYP=2;
                if PACORTYP==1,
                    % correcting for all other channels
                    PACORNSE1=single(NaN([nr,nr])); % NN method
                    PACORNSE2=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
                    PACORNSE3=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
                else
                    dum = 1:length(SIG.PARTIAL.rng{1})-2; % 4 channels, correlation between 1 and 2, 2 other remaining channels
                    clear PAC;
                    for i=1:length(dum),
                        PAC{i} = nchoosek(1:length(dum),i);
                    end
                    for k=1:length(PAC),
                        for m=1:size(PAC{k},1),
                            % correcting for all other channels
                            PACORNSE1{k}{m}=single(NaN([nr,nr])); % NN method
                            PACORNSE2{k}{m}=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
                            PACORNSE3{k}{m}=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
                        end % for m=1:size(PAC{k},1),
                    end % for k=1:length(PAC),
                end
            end
            
            CORSIG=single(NaN([nr,nr]));
            
            % METHOD2: noise correlation: % compute corr in window, then average
            % define instead in nsamples
            % edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
            TwinNC.timwinext=rngNC(2:9);
            TwinNC.twextsmp=[TwinNC.timwinext'-1,TwinNC.timwinext'+1];
            TwinNC.twext=t(TwinNC.twextsmp); % TwinNC.twext
            
            TwinNCbsl.timwinext=rngBSL(1:end);
            TwinNCbsl.twextsmp=[TwinNCbsl.timwinext'-1,TwinNCbsl.timwinext'+1];
            TwinNCbsl.twext=t(TwinNCbsl.twextsmp); % TwinNCbsl.twext
            if size(group1,1)>10,
                
                dum2=single(NaN([nr,nr,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                for tp=1:size(TwinNC.twext,1),
                    dum=nanmean(group1(:,:,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3);
                    dum2(:,:,tp)=corr(dum,'rows','complete');
                end
                CORNSE2=nanmean(dum2,3);
                
                dum2=single(NaN([nr,nr,size(TwinNCbsl.twext,1)])); % average within long time window (0-1s), then correlation
                for tp=1:size(TwinNCbsl.twext,1),
                    dum=nanmean(group1(:,:,TwinNCbsl.twextsmp(tp,1):TwinNCbsl.twextsmp(tp,2)),3);
                    dum2(:,:,tp)=corr(dum,'rows','complete');
                end
                CORNSE5=nanmean(dum2,3);
                
                %tic;
                for i=1:nr,
                    A1=squeeze(group1(:,i,rngNC))'; % V: tim x trl
                    for j=1:nr,%j=2
                        if j>i,
                            A2=squeeze(group1(:,j,rngNC))'; % V: tim x trl
                            dum=corrcov(nancov([A1(:)],[A2(:)]));
                            CORTOT(i,j)=dum(1,2);
                        end
                    end
                end
                %t=toc
                %------------------------------------------------------
                
                % METHOD1: noise correlation
                for i=1:nr,
                    A1=squeeze(group1S(:,i,rngNC))';
                    for j=1:nr,%j=2
                        if j>i,
                            A2=squeeze(group1S(:,j,rngNC))';
                            dum=corrcov(nancov(A1(:),A2(:)));
                            CORNSE1(i,j)=dum(1,2);
                        end
                    end
                end
                
                % METHOD1: noise correlation
                for i=1:nr,
                    A1=squeeze(group1S(:,i,rngBSL))';
                    for j=1:nr,%j=2
                        if j>i,
                            A2=squeeze(group1S(:,j,rngBSL))';
                            dum=corrcov(nancov(A1(:),A2(:)));
                            CORNSE4(i,j)=dum(1,2);
                        end
                    end
                end
                
                CORNSE3=corr(group1AVGoverTIME,'rows','complete');
                %isequalwithequalnans_tolerance(CORNSE3(triu(CORNSE3,1)~=0),dum(triu(CORNSE3,1)~=0),0.001)
                
                CORNSE6=corr(group1AVGoverTIMEbsl,'rows','complete');
                
                if isfield(SIG,'PARTIAL'),
                    for i=1:nr,%i=5
                        A1=squeeze(group1S(:,i,rngNC))';
                        if ~isempty(find(SIG.PARTIAL.rng{1}==i)),
                            rngix=1;
                        else
                            rngix=2;
                        end
                        for j=1:nr,%j=6
                            if ~(i==j),
                                if ismember(j,SIG.PARTIAL.rng{rngix}),
                                    if PACORTYP==2,
                                        chnotherALL=SIG.PARTIAL.rng{rngix}(~ismember(SIG.PARTIAL.rng{rngix},[i,j]))
                                        for k=1:length(PAC),
                                            for m=1:size(PAC{k},1),
                                                chnother = chnotherALL( PAC{k}(m,:) )
                                                dum2=single(NaN([1,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                                                for tp=1:size(TwinNC.twext,1),
                                                    dum=nanmean(group1(:,:,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3);
                                                    dum2(1,tp)=partialcorr(dum(:,i),dum(:,j),dum(:,chnother),'rows','complete');
                                                end
                                                PACORNSE2{k}{m}(i,j)=nanmean(dum2,2);
                                                A2=squeeze(group1S(:,j,rngNC))';
                                                %dum=corrcov(nancov(A1(:),A2(:)));
                                                %corr(A1(:),A2(:),'rows','complete');
                                                A3=[]; for k=1:length(chnother), A3tmp=squeeze(group1S(:,chnother(k),rngNC))'; A3=[A3,A3tmp(:)]; end
                                                PACORNSE1{k}{m}(i,j)=partialcorr(A1(:),A2(:),A3,'rows','complete');
                                                PACORNSE3{k}{m}(i,j)=partialcorr(group1AVGoverTIME(:,i),group1AVGoverTIME(:,j),group1AVGoverTIME(:,chnother),'rows','complete');
                                            end
                                        end
                                    else
                                        dum2=single(NaN([1,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                                        chnother=SIG.PARTIAL.rng{rngix}(~ismember(SIG.PARTIAL.rng{rngix},[i,j]));
                                        for tp=1:size(TwinNC.twext,1),
                                            dum=nanmean(group1(:,:,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3);
                                            dum2(1,tp)=partialcorr(dum(:,i),dum(:,j),dum(:,chnother),'rows','complete');
                                        end
                                        PACORNSE2(i,j)=nanmean(dum2,2);
                                        A2=squeeze(group1S(:,j,rngNC))';
                                        %dum=corrcov(nancov(A1(:),A2(:)));
                                        %corr(A1(:),A2(:),'rows','complete');
                                        A3=[]; for k=1:length(chnother), A3tmp=squeeze(group1S(:,chnother(k),rngNC))'; A3=[A3,A3tmp(:)]; end
                                        PACORNSE1(i,j)=partialcorr(A1(:),A2(:),A3,'rows','complete');
                                        PACORNSE3(i,j)=partialcorr(group1AVGoverTIME(:,i),group1AVGoverTIME(:,j),group1AVGoverTIME(:,chnother),'rows','complete');
                                    end
                                end
                            end % if ~(i==j),
                        end % for j=1:nr,
                    end % for i=1:nr,
                end % if isfield(SIG,'PARTIAL'),
                
                if 0,
                    A=CORNSE3;
                    B=squeeze( TwinPP.S1c(:,:,:,1) );
                    figure;imagesc(A);figure;imagesc(B)
                    A(1:5,1:5),B(1:5,1:5)
                    isequalwithequalnans_tolerance(A,B,0.001)
                end
                % signal corelation
                for i=1:nr, % i=33,j=70
                    A=[group1AVG(i,rngSC)];
                    for j=1:nr,
                        if j>i,
                            B=[group1AVG(j,rngSC)];
                            dum=corrcov(nancov(A',B'));
                            CORSIG(i,j)=dum(1,2);
                        end
                    end
                end
            end % if size(group1,1)>10,
            clear DCOR;
            DCOR.CORTOT=CORTOT;
            DCOR.CORSIG=CORSIG;
            DCOR.CORNSE1=CORNSE1;
            DCOR.CORNSE2=CORNSE2;
            DCOR.CORNSE3=CORNSE3;
            
            DCOR.CORNSE4=CORNSE4;
            DCOR.CORNSE5=CORNSE5;
            DCOR.CORNSE6=CORNSE6;
            if isfield(SIG,'PARTIAL'),
                DCOR.PACORNSE1=PACORNSE1;
                DCOR.PACORNSE2=PACORNSE2;
                DCOR.PACORNSE3=PACORNSE3;
            end
            DCOR.TwinNSE2 = TwinNC;
            if 0,
                figure;imagesc(CORNSE1);colorbar;dum=caxis;
                figure;imagesc(CORNSE2);caxis(dum);
                figure;imagesc(CORNSE3);caxis(dum);
            end
            DPOP{ses}.PL{layer}.DCOR{ccnd}=DCOR;
        end % for ccnd=1:length(SIG.BTRGLABEL),
    end % for layer=1:length(TOT{ses}.SIG.DAT),
    
    %--------------------------------------------------------------------------
    % compute first half second half correlation coefficients
    
    for layer=1:NL,
        
        t=SIG.tax;
        ixSTRT=find(t==0);
        [~,ixBSL] =min(abs( (t--0.5) )); % t(ixBSL)
        [~,ixSTOP]=min(abs( (t-1.1) )); % t(ixSTOP)
        [~,ixSTOP2]=min(abs( (t-1) )); % t(ixSTOP2)
        
        rngSC=ixBSL:ixSTOP; % t(rngSC) include some baseline
        rngNC =ixSTRT:ixSTOP; % t(rngNC)
        rngNC2=ixSTRT:ixSTOP2; % t(rngNC2)
        nr=size(SIG.DAT{layer},2);
        
        % divide trials into first and second half
        CNnumA=cellfun(@length,CN);
        CNnum=CNnumA;
        CNnum(mod(CNnum,2)==1)=CNnum(mod(CNnum,2)==1)-1;
        CNnum=CNnum/2
        clear CNrng;
        for ccnd=1:length(CN),
            CNrng{ccnd}{1}=1:CNnum(ccnd);
            CNrng{ccnd}{2}=CNnumA(ccnd)-CNnum(ccnd)+1:CNnumA(ccnd);
        end
        
        for CNrngix=1:2,
            
            for ccnd=1:length(SIG.BTRGLABEL),
                
                group1=SIG.DAT{layer}(CNrng{ccnd}{CNrngix},:,:); % trl*ROI*time
                group1AVG=nanmean(group1,1); % 1*ROI*time
                group1S=group1-repmat(group1AVG,[size(group1,1),1,1]); % trl*ROI*time
                group1AVG=squeeze(group1AVG); % ROI*time
                
                group1AVGoverTIME=nanmean(group1(:,:,rngNC2),3); % 1*ROI*time t(rngNC2)
                
                CORTOT=single(NaN([nr,nr]));
                CORNSE1=single(NaN([nr,nr])); % NN method
                CORNSE2=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
                CORNSE3=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
                
                if 0,
                    if isfield(SIG,'PARTIAL'),
                        % correcting for all other channels
                        PACORNSE1=single(NaN([nr,nr])); % NN method
                        PACORNSE2=single(NaN([nr,nr])); % compute corr in window, then average 0-1s
                        PACORNSE3=single(NaN([nr,nr])); % average within long time window (0-1s), then correlation
                    end
                end
                
                CORSIG=single(NaN([nr,nr]));
                
                % METHOD2: noise correlation: % compute corr in window, then average
                % define instead in nsamples
                % edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
                TwinNC.timwinext=rngNC(2:9);
                TwinNC.twextsmp=[TwinNC.timwinext'-1,TwinNC.timwinext'+1];
                TwinNC.twext=t(TwinNC.twextsmp); % TwinNC.twext
                if size(group1,1)>10,
                    
                    dum2=single(NaN([nr,nr,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                    for tp=1:size(TwinNC.twext,1),
                        dum=nanmean(group1(:,:,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3);
                        dum2(:,:,tp)=corr(dum,'rows','complete');
                    end
                    CORNSE2=nanmean(dum2,3);
                    
                    for i=1:nr,
                        A1=squeeze(group1(:,i,rngNC))'; % V: tim x trl
                        for j=1:nr,%j=2
                            if j>i,
                                A2=squeeze(group1(:,j,rngNC))'; % V: tim x trl
                                dum=corrcov(nancov([A1(:)],[A2(:)]));
                                CORTOT(i,j)=dum(1,2);
                            end
                        end
                    end
                    %------------------------------------------------------
                    
                    % METHOD1: noise correlation
                    for i=1:nr,
                        A1=squeeze(group1S(:,i,rngNC))';
                        for j=1:nr,%j=2
                            if j>i,
                                A2=squeeze(group1S(:,j,rngNC))';
                                dum=corrcov(nancov(A1(:),A2(:)));
                                CORNSE1(i,j)=dum(1,2);
                            end
                        end
                    end
                    
                    CORNSE3=corr(group1AVGoverTIME,'rows','complete');
                    
                    if 0,
                        if isfield(SIG,'PARTIAL'),
                            for i=1:nr,%i=5
                                A1=squeeze(group1S(:,i,rngNC))';
                                if ~isempty(find(SIG.PARTIAL.rng{1}==i)),
                                    rngix=1;
                                else
                                    rngix=2;
                                end
                                for j=1:nr,%j=6
                                    if ismember(j,SIG.PARTIAL.rng{rngix}),
                                        dum2=single(NaN([1,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                                        chnother=SIG.PARTIAL.rng{rngix}(~ismember(SIG.PARTIAL.rng{rngix},[i,j]));
                                        for tp=1:size(TwinNC.twext,1),
                                            dum=nanmean(group1(:,:,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3);
                                            dum2(1,tp)=partialcorr(dum(:,i),dum(:,j),dum(:,chnother),'rows','complete');
                                        end
                                        PACORNSE2(i,j)=nanmean(dum2,2);
                                        A2=squeeze(group1S(:,j,rngNC))';
                                        %dum=corrcov(nancov(A1(:),A2(:)));
                                        %corr(A1(:),A2(:),'rows','complete');
                                        A3=[]; for k=1:length(chnother), A3tmp=squeeze(group1S(:,chnother(k),rngNC))'; A3=[A3,A3tmp(:)]; end
                                        PACORNSE1(i,j)=partialcorr(A1(:),A2(:),A3,'rows','complete');
                                        PACORNSE3(i,j)=partialcorr(group1AVGoverTIME(:,i),group1AVGoverTIME(:,j),group1AVGoverTIME(:,chnother),'rows','complete');
                                    end
                                end % for j=1:nr,
                            end % for i=1:nr,
                        end % if isfield(SIG,'PARTIAL'),
                    end
                    
                    % signal corelation
                    for i=1:nr, % i=33,j=70
                        A=[group1AVG(i,rngSC)];
                        for j=1:nr,
                            if j>i,
                                B=[group1AVG(j,rngSC)];
                                dum=corrcov(nancov(A',B'));
                                CORSIG(i,j)=dum(1,2);
                            end
                        end
                    end
                end % if size(group1,1)>10,
                clear DCOR;
                DCOR.CORTOT=CORTOT;
                DCOR.CORSIG=CORSIG;
                DCOR.CORNSE1=CORNSE1;
                DCOR.CORNSE2=CORNSE2;
                DCOR.CORNSE3=CORNSE3;
                
                if 0,
                    if isfield(SIG,'PARTIAL'),
                        DCOR.PACORNSE1=PACORNSE1;
                        DCOR.PACORNSE2=PACORNSE2;
                        DCOR.PACORNSE3=PACORNSE3;
                    end
                end
                DCOR.TwinNSE2 = TwinNC;
                
                DPOP{ses}.PL{layer}.DCORrng{ccnd}{CNrngix}=DCOR;
            end % for ccnd=1:length(SIG.BTRGLABEL),
            
        end % for CNrng=1:2,
    end % for layer=1:length(TOT{ses}.SIG.DAT),
    
end % if 0,


%--------------------------------------------------------------------------

stop=0;
if isfield(SIG,'label'),
    if (~isempty(strmatch('pc',SIG.label)))|(~isempty(strmatch('sm',SIG.label))),
        stop=1;
    end
end;
if ~stop, % population activity etc
    
    %
    % see edit rsbs_popactivity
    
    % find cell type of each cell
    if 1,
        SIG.POOLCV = POOLCV;
        CELLLAB=NaN(size(SIG.POOLCV,1),1); % index into TL
        for chn=1:size(SIG.POOLCV,1),
            ix=find( (LN.TL.rix==RIXSES(1))&(LN.TL.ses==RIXSES(2))...
                &(LN.TL.lay==SIG.POOLCV(chn,1))&(LN.TL.B==SIG.POOLCV(chn,2)) );
            if not(isempty(ix)),
                CELLLAB(chn) = LN.TL.labJ(ix);
            else
                warning(sprintf('cannot find chn %d',chn));
            end
        end
        
        CLASSLAB={'','PV','SOM','VIP'};
        CLASSLAB2={'PYR','PVB','SOM','VIP'};
        CLASSNUM=zeros(size(CELLLAB));
        for CLASS=1:length(CLASSLAB),
            CLASSNUM(CELLLAB==strmatch(CLASSLAB{CLASS},LN.TL.labB,'exact'))=CLASS;
        end
    end
    
    for layer=1:NL,
        
        t=SIG.tax;
        ixSTRT=find(t==0);
        [~,ixBSL] =min(abs( (t--0.5) ));
        [~,ixBSL2] =min(abs( (t--1) )); % t(ixBSL2)
        [~,ixSTOP]=min(abs( (t-1.1) )); % t(ixSTOP)
        [~,ixSTOP2]=min(abs( (t-1) )); % t(ixSTOP2)
        
        rngSC=ixBSL:ixSTOP; % t(rngSC) include some baseline
        rngNC =ixSTRT:ixSTOP; % t(rngNC)
        rngNC2=ixSTRT:ixSTOP2; % t(rngNC)
        nr=size(SIG.DAT{layer},2);
        
        rngBSL=ixBSL2:ixSTRT-1; % t(rngBSL)
        
        for ccnd=1:length(SIG.BTRGLABEL),
            
            group1=SIG.DAT{layer}(CN{ccnd},:,:); % trl*ROI*time
            group1AVG=nanmean(group1,1); % 1*ROI*time
            group1S=group1-repmat(group1AVG,[size(group1,1),1,1]); % trl*ROI*time
            group1AVG=squeeze(group1AVG); % ROI*time
            
            group1AVGoverTIME=nanmean(group1(:,:,rngNC2),3); % 1*ROI*time
            group1AVGoverTIMEbsl=nanmean(group1(:,:,[rngBSL,rngBSL(end)+1]),3); % 1*ROI*time t([rngBSL,rngBSL(end)+1])
            
            PA1=single(NaN([nr,2])); % NN method: 2 refers to regular PopCor and CorOfSummedActivity
            PA2=single(NaN([nr,2])); % compute corr in window, then average 0-1s
            
            PA1ct=single(NaN([nr,2,4])); % NN method: 2 refers to regular PopCor and CorOfSummedActivity
            PA2ct=single(NaN([nr,2,4])); % compute corr in window, then average 0-1s
            
            % METHOD2: noise correlation: % compute corr in window, then average
            % define instead in nsamples
            % edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
            TwinNC.timwinext=rngNC(2:9);
            TwinNC.twextsmp=[TwinNC.timwinext'-1,TwinNC.timwinext'+1];
            TwinNC.twext=t(TwinNC.twextsmp); % TwinNC.twext
            
            if size(group1,1)>10,
                
                nchn=size(group1,2);
                pcA=single(NaN([nchn,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                pcB=single(NaN([nchn,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                for chn=1:nchn, % for each channel
                    if mod(chn,100)==0, fprintf('%d(of %d)\n',chn,nchn); end
                    for tp=1:size(TwinNC.twext,1),
                        A=nanmean(group1S(:,chn,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3)';
                        chnothr=setdiff([1:nchn],chn);
                        B=nanmean(group1S(:,chnothr,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3)';
                        
                        SD=nanstd(A);
                        SF=repmat(A,[length(chnothr),1]).*B;
                        pcA(chn,tp)=1/SD * nansum(nansum(SF,1));
                        %Note that ci is proportional to the Pearson correlation of
                        %cell i with the summed activity of all other neurons
                        [pcB(chn,tp)]=corr(A', nansum(B,1)','rows','complete');
                    end
                end % for chn=1:size(out.POOLDAT,1),
                
                PA1(:,1) = nanmean(pcA,2);
                PA1(:,2) = nanmean(pcB,2);
                
                %------------------------------------------------------
                % METHOD1: noise correlation
                
                for chn=1:nchn, % for each channel
                    if mod(chn,100)==0, fprintf('%d(of %d)\n',chn,nchn); end
                    A=group1S(:,chn,rngNC);
                    chnothr=setdiff([1:nchn],chn);
                    B=squeeze(group1S(:,chnothr,rngNC)); % mean across channels
                    
                    SD=nanstd(A(:)');
                    SF=repmat(A,[1,length(chnothr),1]).*B;
                    PA2(chn,1)=1/SD * nansum(nansum(nansum(SF)));
                    %Note that ci is proportional to the Pearson correlation of
                    %cell i with the summed activity of all other neurons
                    An=squeeze(A); Bn=squeeze(nansum(B,2));
                    [PA2(chn,2)]=corr(An(:), Bn(:),'rows','complete');
                end % for chn=1:size(out.POOLDAT,1),
                
                if 0,
                    figure;plot(PA1(:,1),PA1(:,2));
                    figure;plot(PA2(:,1),PA2(:,2));
                    figure;plot(PA1(:,2),PA2(:,2));
                end
                
                for lab=1:4,
                    
                    nchn=size(group1,2);
                    pcA=single(NaN([nchn,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                    pcB=single(NaN([nchn,size(TwinNC.twext,1)])); % average within long time window (0-1s), then correlation
                    for chn=1:nchn, % for each channel
                        if mod(chn,100)==0, fprintf('lab%d:%d(of %d)\n',lab,chn,nchn); end
                        for tp=1:size(TwinNC.twext,1),
                            A=nanmean(group1S(:,chn,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3)';
                            %chnothr=setdiff([1:nchn],chn);
                            chnothr=setdiff(find(CLASSNUM==lab),chn);
                            B=nanmean(group1S(:,chnothr,TwinNC.twextsmp(tp,1):TwinNC.twextsmp(tp,2)),3)';
                            
                            SD=nanstd(A);
                            SF=repmat(A,[length(chnothr),1]).*B;
                            pcA(chn,tp)=1/SD * nansum(nansum(SF,1));
                            %Note that ci is proportional to the Pearson correlation of
                            %cell i with the summed activity of all other neurons
                            [pcB(chn,tp)]=corr(A', nansum(B,1)','rows','complete');
                        end
                    end % for chn=1:size(out.POOLDAT,1),
                    
                    PA1ct(:,1,lab) = nanmean(pcA,2);
                    PA1ct(:,2,lab) = nanmean(pcB,2);
                    
                    for chn=1:nchn, % for each channel
                        if mod(chn,100)==0, fprintf('lab%d:%d(of %d)\n',lab,chn,nchn); end
                        A=group1S(:,chn,rngNC);
                        %chnothr=setdiff([1:nchn],chn);
                        chnothr=setdiff(find(CLASSNUM==lab),chn);
                        %    B=squeeze(group1S(:,chnothr,rngNC)); % mean across channels
                        %if length(chnothr)==1,
                        %   B=squeeze(group1S(:,chnothr,rngNC)); % mean across channels
                        %else
                        B=group1S(:,chnothr,rngNC); % mean across channels
                        %end
                        SD=nanstd(A(:)');
                        SF=repmat(A,[1,length(chnothr),1]).*B;
                        PA2ct(chn,1,lab)=1/SD * nansum(nansum(nansum(SF)));
                        %Note that ci is proportional to the Pearson correlation of
                        %cell i with the summed activity of all other neurons
                        An=squeeze(A); Bn=squeeze(nansum(B,2));
                        %[PA2ct(chn,2,lab)]=corr(An(:), Bn(:),'rows','complete');
                        PA2ct(chn,2,lab)=corr(An(:), Bn(:),'rows','complete');
                    end % for chn=1:size(out.POOLDAT,1),
                end % for lab=1:4,
                
            end % if size(group1,1)>10,
            clear DCOR;
            DCOR.PA1=PA1;
            DCOR.PA2=PA2;
            DCOR.PA1ct=PA1ct;
            DCOR.PA2ct=PA2ct;
            DCOR.TwinNSE2 = TwinNC;
            DCOR.t     = t;
            DCOR.rngNC = rngNC;
            DPOP{ses}.PL{layer}.PA{ccnd}=DCOR;
        end % for ccnd=1:length(SIG.BTRGLABEL),
    end % for layer=1:length(TOT{ses}.SIG.DAT),
    
end % if stop,

%test = 1;
if 0,
    close all;
    figure;plot(DPOP{ses}.PL{layer}.PA{1}.PA1(:,2),DPOP{ses}.PL{layer}.PA{1}.PA1ct(:,2,1),'bo')
    figure;plot(DPOP{ses}.PL{layer}.PA{1}.PA1(:,2),DPOP{ses}.PL{layer}.PA{1}.PA1ct(:,2,2),'bo')
    figure;plot(DPOP{ses}.PL{layer}.PA{1}.PA1(:,2),DPOP{ses}.PL{layer}.PA{1}.PA1ct(:,2,3),'bo')
    figure;plot(DPOP{ses}.PL{layer}.PA{1}.PA1(:,2),DPOP{ses}.PL{layer}.PA{1}.PA1ct(:,2,4),'bo')
end
%--------------------------------------------------------------------------

% ramping stuff
if 0,
    for layer=1:NL,
        
        Twin.twextsmp=[];
        ix0=find(t==0);
        Twin.twextsmp=[ix0-5 ix0-2; ix0-1 ix0+2; ix0+4 ix0+7]; % Twin.twextsmp(:,1):Twin.twextsmp(:,2)
        Twin.twext=t(Twin.twextsmp);
        
        clear TwinPP;
        TwinPP.twext = Twin.twext;
        TwinPP.twextsmp = Twin.twextsmp;
        
        TwinPP.twref=[-0.5,-0.1; -0.1 0.3; 0.3 0.8];
        
        TwinPP.C=nchoosek(1:size(TwinPP.twref,1),2);
        nr=size(SIG.DAT{layer},2);
        nw=size(TwinPP.C,1);
        ncnd=length(SIG.BTRGLABEL);
        TwinPP.m1        = single(NaN([nw,nr,ncnd]));
        TwinPP.m2        = single(NaN([nw,nr,ncnd]));
        TwinPP.v1        = single(NaN([nw,nr,ncnd]));
        TwinPP.v2        = single(NaN([nw,nr,ncnd]));
        TwinPP.n1        = single(NaN([nw,nr,ncnd]));
        TwinPP.n2        = single(NaN([nw,nr,ncnd]));
        TwinPP.S1c       = single(NaN([nw,nr,nr,ncnd]));
        TwinPP.S2c       = single(NaN([nw,nr,nr,ncnd]));
        TwinPP.DPi       = single(NaN([nw,nr,ncnd]));
        TwinPP.p         = single(NaN([nw,nr,ncnd]));
        TwinPP.md        = single(NaN([nw,nr,ncnd]));
        for w=1:nw,%w=43,TwinPP.twext(w,:)
            for ccnd=1:length(SIG.BTRGLABEL),
                rng1=TwinPP.twextsmp(TwinPP.C(w,1),1):TwinPP.twextsmp(TwinPP.C(w,1),2);
                rng2=TwinPP.twextsmp(TwinPP.C(w,2),1):TwinPP.twextsmp(TwinPP.C(w,2),2);
                
                group1=nanmean(SIG.DAT{layer}(CN{ccnd},:,rng1),3)';
                group2=nanmean(SIG.DAT{layer}(CN{ccnd},:,rng2),3)';
                
                if size(group1,2)>10&size(group2,2)>10,
                    % compute covariance matrix, in matlab each row is an observation, and each column a variable
                    S1=nancov(group1');   S2=nancov(group2');
                    try,
                        S1c=corrcov(S1);      S2c=corrcov(S2);
                    catch
                        S1(isnan(S1))=0;
                        S2(isnan(S2))=0;
                        S1c=corrcov(S1);
                        S2c=corrcov(S2);
                    end
                    
                    m1=nanmean(group1,2); m2=nanmean(group2,2);
                    n1=sum(~isnan(group1),2); n2=sum(~isnan(group2),2);
                    
                    TwinPP.m1(w,:,ccnd)      = m1;
                    TwinPP.m2(w,:,ccnd)      = m2;
                    TwinPP.n1(w,:,ccnd)      = n1;
                    TwinPP.n2(w,:,ccnd)      = n2;
                    TwinPP.S1c(w,:,:,ccnd)   = S1c;
                    TwinPP.S2c(w,:,:,ccnd)   = S2c;
                    
                    % compute individual d-primes
                    v1=nanvar(group1,[],2); v2=nanvar(group2,[],2);
                    poolstd=sqrt( (v1.*(n1-1)+v2.*(n2-1))./(n1+n2-2) ); % pooledstd
                    DPi=(m1-m2)./poolstd;
                    TwinPP.DPi(w,:,ccnd) = DPi; % sum(DPi.^2) % is equal to d_sq_shuf
                    
                    TwinPP.v1(w,:,ccnd)      = v1;
                    TwinPP.v2(w,:,ccnd)      = v2;
                    
                    % auroc
                    for i=1:size(group1,1),
                        A=group1(i,:); B=group2(i,:); A=A(~isnan(A)); B=B(~isnan(B));
                        %TwinPP.auc(w,i,ccnd)=calc_ROC(A,B); % Twin.auc(w,:)
                        TwinPP.p(w,i,ccnd)  =ranksum(A,B);
                        TwinPP.md(w,i,ccnd) =nanmedian(group1(i,:))-nanmedian(group2(i,:));
                    end
                end % if size(group1,2)>10&size(group1,2)>10,
            end % for ccnd=1:2,
        end % for w=1:nw,
        DPOP{ses}.PL{layer}.TwinRAMP=TwinPP;
    end % for layer=1:length(TOT{ses}.SIG.DAT),
end % if 0, % ramping stuff

DPOP{ses}.POOLCV = POOLCV;