function [CORLAG]=rsbs_extra_data_trigger_popact_model_corrlag(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
% edit rsbs_predict_resp_corrlag.m
% edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
% - [LN.CORLAG]=rsbsI_sd_pop_signoise_makeCMB_CORLAG(LN.TOTALL,LN.TL,LN.RSG,LN.CMB);


%rsbs_placecorr
% - [LN.TCOR]=rsbsI_rd_pop_makeTCOR(LN.TOTALL,LN.TL,LN.RSG,LN.CMB,PARAM)
% - rsbsI_rd_pop_corlag(LN.TOTALL,LN.RSG,LN.TL,LN.DPOP,LN.CMB,VRSSEL,LN.TCOR);

%[CLAG]=rsbs_predict_resp_corrlag(out,ADDINF);
%[LN.CORLAG]=rsbsI_sd_pop_signoise_makeCMB_CORLAG(LN.TOTALL,LN.TL,LN.RSG,LN.CMB);

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

%--------------------------------------------------------------------------
% preprocess data

if 1,
    TEST = 1;
    
    t = SIG.tax;
    ix0   = find(t==0)
    ixbsl = find(t>-1,1,'first');
    ixstm = find(t>1,1,'first');
    
    ixbsl = find(t>-4,1,'first');
    ixstm = find(t>4,1,'first');
    
    trng = [fliplr([ix0-4:-4:ixbsl]), [ix0:4:ixstm]] ;
    t( trng )
    navg = 4;
    
    DATRAW = SIGPOP.DAT{1}; % trl x cnd x time
    DAT    = NaN([size(DATRAW,1),size(DATRAW,2),length(trng)]);
    for nsmp=1:length(trng),
        DAT(:,:,nsmp)=nanmean(DATRAW(:,:,trng(nsmp):trng(nsmp)+navg-1),3);
    end

    clear CN;
    for c=1:length(SIG.BTRGLABEL),
        CN{c}=find(SIG.sel.cond==c & SIG.sel.valid==1);
    end
    t=SIG.tax;
    
    if TEST
        if strcmp(idinf{NID}.type,'SWITCH'),
            sel = ismember(SIG.BTRGLABEL,{'Von','Aon','VonIrr','AonIrr'});
        else
            sel = ismember(SIG.BTRGLABEL,{'Von','Aon','ons'});
        end
        for i=1:length(sel),
            if sel(i)==0, CN{i}=[]; end
        end
    end
    tim=t( trng );
    
end % if 1, % preprocess data

%--------------------------------------------------------------------------
% pick training and testing trials

numfolds = 1; clear ALLinds;
for ccnd=1:length(SIG.BTRGLABEL),
    if length(CN{ccnd})>10,
        Ntrials = length(CN{ccnd});
        rperm = randperm(Ntrials);
        n15 = ceil(Ntrials/5);
        for i = 1:numfolds
            ALLinds{ccnd}(i).test = (i-1)*n15 + (1:n15);
            ALLinds{ccnd}(i).test(ALLinds{ccnd}(i).test>Ntrials) = [];
            ALLinds{ccnd}(i).test = rperm(ALLinds{ccnd}(i).test);
            ALLinds{ccnd}(i).train = find(~ismember(1:Ntrials , ALLinds{ccnd}(i).test ));
        end % for i = 1:numfolds
    end
end

%%--------------------------------------------------------------------------
if 1,
    
    clear DEP IND;
    DEP{1} = strmatch('smPYR',SIGPOP.label,'exact'); % predicted variable
    IND{1} = [strmatch('smPVB',SIGPOP.label,'exact'),strmatch('smSOM',SIGPOP.label,'exact'),strmatch('smVIP',SIGPOP.label,'exact')];
    
    DEP{2} = strmatch('smPVB',SIGPOP.label,'exact'); % predicted variable
    IND{2} = [strmatch('smPYR',SIGPOP.label,'exact'),strmatch('smSOM',SIGPOP.label,'exact'),strmatch('smVIP',SIGPOP.label,'exact')];
    
    DEP{3} = strmatch('smSOM',SIGPOP.label,'exact'); % predicted variable
    IND{3} = [strmatch('smPYR',SIGPOP.label,'exact'),strmatch('smPVB',SIGPOP.label,'exact'),strmatch('smVIP',SIGPOP.label,'exact')];
    
    DEP{4} = strmatch('smVIP',SIGPOP.label,'exact'); % predicted variable
    IND{4} = [strmatch('smPYR',SIGPOP.label,'exact'),strmatch('smPVB',SIGPOP.label,'exact'),strmatch('smSOM',SIGPOP.label,'exact')];
    
    
    alldts{1} = -5:-1; % preceding samples
    alldts{1} = -15:5:-1; % preceding samples
    
    clear CORLAG; % estimate model
    for LAGTYP = 1:length(alldts),
        
        dts = alldts{LAGTYP};
        
        BSLSUBTRACT=1;
        CORLAG{LAGTYP}.BSLSUBTRACT = BSLSUBTRACT;
        CORLAG{LAGTYP}.DEP       = DEP;
        CORLAG{LAGTYP}.IND       = IND;
        CORLAG{LAGTYP}.BTRGLABEL = SIG.BTRGLABEL;
        CORLAG{LAGTYP}.label     = SIGPOP.label;
        CORLAG{LAGTYP}.t         = tim;
        CORLAG{LAGTYP}.dts       = dts;
        
        CORLAG{LAGTYP}.CNnum = cellfun(@length,CN);
        CORLAG{LAGTYP}.sel   = sel;
        
        Nroi    = size(DAT,2);
        Nsmp    = size(DAT,3);
        rngix   = -min(dts(1))+1:Nsmp;
        CORLAG{LAGTYP}.rngix = rngix;
        
        for ccnd=1:length(SIG.BTRGLABEL),
            if length(CN{ccnd})>10,
                fprintf('%s ccnd=%d\n',datestr(now),ccnd);
                
                CURDAT = DAT(CN{ccnd},:,:); % trl x ROI x time
                                 
                if BSLSUBTRACT,
                    CURDAT = CURDAT-repmat(nanmean(CURDAT,1),[size(CURDAT,1),1,1]);                    
                end
                
                Ntrials = size(CURDAT,1);
                
                inds=ALLinds{ccnd};
                
                for MODEL=1:length(DEP),
                    CURDEP = DEP{MODEL};
                    
                    % take all combinations of predictors
                    clear PREDICTOR;
                    for MODELPREDICTORNUM=1:length(IND{MODEL}),  % for each number of predictors (1 or 2 indep variables)
                        % NB: correlation: can only have 1 predictor!! i==1
                        PREDICTOR{MODELPREDICTORNUM} = nchoosek(IND{MODEL},MODELPREDICTORNUM);
                    end
                    
                    for MODELPREDICTORNUM=1:length(PREDICTOR), % for each number of predictors
                        for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each version with given number of predictors
                            fprintf('DEP[%s] predsiz%d CURPRED%d[%s]\n',  sprintf('%d ',DEP{MODEL}), MODELPREDICTORNUM,MODELPREDICTOR, sprintf('%d ',PREDICTOR{MODELPREDICTORNUM}(MODELPREDICTOR,:)));
                        end % for each version with given number of predictors
                    end % for each number of predictors
                    
                    for MODELPREDICTORNUM=1:length(PREDICTOR), % for each model type with specific number of predictors
                        %MODELPREDICTORNUM = 2;
                        
                        for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each set of predictors
                            
                            CURPREDICTORSET = PREDICTOR{MODELPREDICTORNUM}(MODELPREDICTOR,:);
                            
                            for i = 1:numfolds,
                                testind  = inds(i).test;
                                trainind = inds(i).train;
                                
                                % VRS=1 seperately for each sample
                                % VRS=2 for all samples together
                                for VRS=2%:2, % VRS=2;
                                    clear RT;
                                    if VRS==1, % seperately for each sample
                                        datYtst  = NaN([Nsmp,length(testind)]);
                                        datYprd = NaN([Nsmp,length(testind)]);
                                        
                                        for S = -min(dts)+1:Nsmp-max(dts), % for each sample
                                            % other predictors at times preceding time S
                                            X0trn = CURDAT(trainind,CURPREDICTORSET,S+dts);
                                            X0tst = CURDAT(testind,CURPREDICTORSET,S+dts);
                                            Xtrn = reshape(X0trn,[size(X0trn,1),size(X0trn,2)*size(X0trn,3)]);
                                            % squeeze(X0trn(12,:,:)), Xtrn(12,:)
                                            Xtst = reshape(X0tst,[size(X0tst,1),size(X0tst,2)*size(X0tst,3)]);
                                            
                                            % activity at time S (to be predicted)
                                            Ytrn = CURDAT(trainind,CURDEP,S);
                                            Ytst = CURDAT(testind,CURDEP,S);
                                            
                                            datYtst(S,:) = Ytst;
                                            
                                            eoForest = TreeBagger(32, Xtrn, Ytrn, 'Method', 'regression', ... % X Each row represents an observation and each column represents a predictor or feature
                                                'OOBPred','On', 'MinLeaf', 5);
                                            Yprd = eoForest.predict(Xtst);
                                            datYprd(S,:) = Yprd;
                                        end % % for each sample
                                    else % combine data of different samples
                                        datYtst  = NaN([1,length(testind)]);
                                        datYprd = NaN([1,length(testind)]);
                                        
                                        Xtrn = NaN(length(trainind)*Nsmp,length(CURPREDICTORSET)*length(dts));
                                        Xtst = NaN(length(testind)*Nsmp,length(CURPREDICTORSET)*length(dts));
                                        Ytrn = NaN(length(trainind)*Nsmp,1);
                                        Ytst = NaN(length(testind)*Nsmp,1);
                                        for S = -min(dts(1))+1:Nsmp,
                                            % other predictors at times preceding time S
                                            X0trn = CURDAT(trainind,CURPREDICTORSET,S+dts);
                                            X0tst = CURDAT(testind,CURPREDICTORSET,S+dts);
                                            Xtrn((S-1)*length(trainind)+1:(S)*length(trainind),:) = reshape(X0trn,[size(X0trn,1),size(X0trn,2)*size(X0trn,3)]);
                                            %dum=reshape(X0trn,[size(X0trn,1),size(X0trn,2)*size(X0trn,3)]);
                                            % squeeze(X0trn(12,:,:)),dum(12,:)
                                            Xtst((S-1)*length(testind)+1:(S)*length(testind),:) = reshape(X0tst,[size(X0tst,1),size(X0tst,2)*size(X0tst,3)]);
                                            
                                            % activity at time S (to be predicted)
                                            Ytrn((S-1)*length(trainind)+1:(S)*length(trainind),:) = CURDAT(trainind,CURDEP,S);
                                            Ytst((S-1)*length(testind)+1:(S)*length(testind),:) = CURDAT(testind,CURDEP,S);
                                            % dum =CURDAT(trainind,CURDEP,S);
                                        end % for S = length(dts)+1:length(Nsmp),
                                        datYtst = Ytst;
                                        eoForest = TreeBagger(32,Xtrn,Ytrn,'Method', 'regression', ... % X Each row represents an observation and each column represents a predictor or feature
                                            'OOBPred','On', 'MinLeaf', 5);
                                        Yprd = eoForest.predict(Xtst);
                                        datYprd = Yprd;
                                        
                                    end % if VRS==1, % seperately for each sample
                                    
                                    CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ytst{i}{VRS} = datYtst;
                                    CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Yprd{i}{VRS} = datYprd;
                                    
                                end % for VRS=1:2,
                                CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ntrial          = [length(inds(i).train),length(inds(i).test)];
                                
                            end % end % for i = 1:numfolds,
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURDEP          = CURDEP;
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURPREDICTORSET = CURPREDICTORSET;
                            
                        end % for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each set of predictors
                    end % for MODELPREDICTORNUM=1:length(IND{MODEL}), % for each model type with specific number of predictors
                end % for lab=1:4,
            end % if length(CN{ccnd})>10,
        end % for ccnd=1:length(SIG.BTRGLABEL),
        
    end % for LAGTYP = 1:length(dts),
    
end

if 0,
    clc
    RSG.savedirplot = 'M:\BehaviourImaging\rsbs_workflow\rsbs_trigger_popact_model_predict\';
    if ~isdir(RSG.savedirplot), mkdir(RSG.savedirplot); end
    
    interactive=0;toplot=1;
    ext ='-dtiff';fext    ='.tiff';
    %ext='-dpsc';fext='.ps';
    
    if interactive, vis='on'; else vis='off'; end
    
    for LAGTYP = 1:length(CORLAG),
        for ccnd=1:length(CORLAG{LAGTYP}.dat),
            fprintf('%%-------------------\nccnd%d\n',ccnd);
            if ~isempty(CORLAG{LAGTYP}.dat{ccnd}),
                for MODEL=1:length(CORLAG{LAGTYP}.dat{ccnd}),
                    cix1 = 1;
                    cix2 = 1;
                    for MODELPREDICTORNUM=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}), % MODELPREDICTORNUM=2 MODELPREDICTORNUM=3
                        for MODELPREDICTOR=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}),
                            
                            %ccnd=1;MODEL=1;MODELPREDICTORNUM=2;MODELPREDICTOR=1;
                            
                            VRS=2;
                            i=1;
                            CURDEP          = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURDEP;
                            CURPREDICTORSET = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURPREDICTORSET;
                            datYtst = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ytst{i}{VRS};
                            datYprd = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Yprd{i}{VRS};
                            Ntrial       = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ntrial;
                            t            = CORLAG{LAGTYP}.t;
                            rngix        = CORLAG{LAGTYP}.rngix;
                            
                            Ytst = datYtst(:);
                            Yprd = datYprd(:);
                            SSE = nansum((Ytst-Yprd).^2); SST = nansum((Ytst-repmat(nanmean(Ytst),[size(Ytst,1),1])).^2); rsquare = 1-SSE./SST;
                            
                            str=[];
                            for i=1:length(CURPREDICTORSET),
                                str=[str,sprintf('%s,',SIGPOP.label{CURPREDICTORSET(i)})];
                            end
                            str(end)='';
                            fprintf('DEP[%s] predsiz%d CURPRED%d[%s], R=%1.2f\n',  sprintf('%s ',SIGPOP.label{CURDEP}), MODELPREDICTORNUM,MODELPREDICTOR, str,rsquare);
                            
                            
                        end % MODELPREDICTOR=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}),
                    end % MODELPREDICTORNUM=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}), % MODELPREDICTORNUM=2 MODELPREDICTORNUM=3
                    
                end %  for MODEL=1:length(DEP),
            end % if ~isempty(CORLAG{LAGTYP}.dat{ccnd}),
        end % for ccnd=1:length(CORLAG{LAGTYP}.dat),
    end % for LAGTYP = 1:length(dts),
    
end % if 1, % plotting


