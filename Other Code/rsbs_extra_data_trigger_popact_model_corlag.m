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
            sel = ismember(SIG.BTRGLABEL,{'Von','Aon','Ron','ons','onsetcircle','onsetgrey'});
        end
        for i=1:length(sel),
            if sel(i)==0, CN{i}=[]; end
        end
    end
    tim=t( trng );
    
    % only use sssions with data for each cell type
    if any(isnan(nanmean(nanmean( DAT(:,find(ismember(SIGPOP.label,{'smPYR','smPVB','smSOM','smVIP'})),:) ,1),3)))
        for c=1:length(SIG.BTRGLABEL),
            CN{c}=[];
        end
    end
    
end % if 1, % preprocess data

%%--------------------------------------------------------------------------
if any(~cellfun(@isempty,CN)),
    
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
    alldts{1} = -15:5:15; % preceding samples
    
    clear CORLAG; % estimate model
    for LAGTYP = 1:length(alldts),
        
        dts = alldts{LAGTYP};
        
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
                
                Ntrials = size(CURDAT,1);
                
                for MODEL=1:length(DEP),
                    CURDEP = DEP{MODEL};
                    
                    % take all combinations of predictors
                    clear PREDICTOR PAC;
                    for MODELPREDICTORNUM=1%:length(IND{MODEL}),  % for each number of predictors (1 or 2 indep variables)
                        % NB: correlation: can only have 1 predictor!! i==1
                        PREDICTOR{MODELPREDICTORNUM} = nchoosek(IND{MODEL},MODELPREDICTORNUM);
                        for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1),  % for each version with given number of predictors
                            dum=setdiff(IND{MODEL},PREDICTOR{MODELPREDICTORNUM}(MODELPREDICTOR,:)); % find remaining predictors for partial corr
                            clear tmp; tix=1;
                            tmp{tix} = [NaN]; tix=tix+1; % default: no to be partialized variables
                            for k=1:length(dum),
                                dum2=nchoosek(dum,k);
                                for l=1:size(dum2,1),
                                    tmp{tix} = dum2(l,:); tix=tix+1;
                                end
                            end
                            PAC{MODELPREDICTORNUM}{MODELPREDICTOR} = tmp;
                        end
                    end
                    
                    for MODELPREDICTORNUM=1:length(PREDICTOR), % for each number of predictors
                        for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each version with given number of predictors
                            if ~isempty(PAC{MODELPREDICTORNUM}{MODELPREDICTOR}), % for each combination of to be partialized variables
                                tmp=PAC{MODELPREDICTORNUM}{MODELPREDICTOR};
                                for k=1:length(tmp),
                                    fprintf('DEP[%s] predsiz%d CURPRED%d[%s] PART[%s] \n',  sprintf('%d ',DEP{MODEL}), MODELPREDICTORNUM,MODELPREDICTOR, sprintf('%d ',PREDICTOR{MODELPREDICTORNUM}(MODELPREDICTOR,:)), sprintf('%d ',tmp{k}));
                                end
                            end % for each combination of to be partialized variables
                        end % for each version with given number of predictors
                    end % for each number of predictors
                    
                    for MODELPREDICTORNUM=1:length(PREDICTOR), % for each model type with specific number of predictors
                        %MODELPREDICTORNUM = 2;
                        
                        for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each set of predictors
                            
                            CURPREDICTORSET = PREDICTOR{MODELPREDICTORNUM}(MODELPREDICTOR,:);
                            CURPAC = PAC{MODELPREDICTORNUM}{MODELPREDICTOR}; % versions with differnt to be partialized
                            
                            % VRS=1 seperately for each sample
                            % VRS=2 for all samples together
                            for VRS=1:2, % VRS=2;
                                clear RT;
                                if VRS==1, % seperately for each sample
                                    for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                        RT{VRS}{MODELTOPARTIALIZEOUT}=NaN([Nsmp,length(dts)]);
                                    end
                                    for S = -min(dts)+1:Nsmp-max(dts), % for each sample
                                        for lag = 1:length(dts),
                                            % autopredictors at times preceding time S
                                            X = squeeze(CURDAT(:,CURPREDICTORSET,S+dts(lag))); % predictors
                                            Y = CURDAT(:,CURDEP,S); % dependent variable
                                            for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                                if isnan(CURPAC{MODELTOPARTIALIZEOUT}),
                                                    RT{VRS}{MODELTOPARTIALIZEOUT}(S,lag) = corr(X,Y,'rows','complete');
                                                else
                                                    Z  = squeeze(CURDAT(:,CURPAC{MODELTOPARTIALIZEOUT},S+dts(lag))); % partializing variables
                                                    for i=1:length(dts),
                                                        RT{VRS}{MODELTOPARTIALIZEOUT}(S,lag) = partialcorr(X,Y,Z,'rows','complete');
                                                    end
                                                end
                                            end % for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                        end % for lag = 1:length(dts),
                                    end % % for each sample
                                else % combine data of different samples
                                    for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                        RT{VRS}{MODELTOPARTIALIZEOUT}=NaN([1,length(dts)]);
                                    end
                                    for lag = 1:length(dts),
                                        X=[];
                                        Y=[];
                                        clear Z;
                                        for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                            Z{MODELTOPARTIALIZEOUT}=[];
                                        end
                                        for S = -min(dts)+1:Nsmp-max(dts), % for each sample
                                            Y = cat(1,Y,CURDAT(:,CURDEP,S)); % dependent variable
                                            X = cat(1,X,CURDAT(:,CURPREDICTORSET,S+dts(lag))); % predictors
                                            for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                                if ~isnan(CURPAC{MODELTOPARTIALIZEOUT}),
                                                    Z{MODELTOPARTIALIZEOUT}=cat(1,Z{MODELTOPARTIALIZEOUT},squeeze(CURDAT(:,CURPAC{MODELTOPARTIALIZEOUT},S+dts(lag)))); % partializing variables
                                                end
                                            end
                                        end % for each sample
                                        for MODELTOPARTIALIZEOUT=1:length(CURPAC),
                                            if isnan(CURPAC{MODELTOPARTIALIZEOUT}),
                                                RT{VRS}{MODELTOPARTIALIZEOUT}(lag) = corr(X,Y,'rows','complete');
                                            else
                                                RT{VRS}{MODELTOPARTIALIZEOUT}(lag) = partialcorr(X,Y,Z{MODELTOPARTIALIZEOUT},'rows','complete');
                                            end
                                        end
                                    end % for lag = 1:length(dts),
                                end % if VRS==1, % seperately for each sample
                            end % for VRS=1:2,
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURDEP          = CURDEP;
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURPREDICTORSET = CURPREDICTORSET;
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURPAC          = CURPAC;
                            CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.RT              = RT;
                            
                        end % for MODELPREDICTOR=1:size(PREDICTOR{MODELPREDICTORNUM},1), % for each set of predictors
                    end % for MODELPREDICTORNUM=1:length(IND{MODEL}), % for each model type with specific number of predictors
                end % for lab=1:4,
            end % if length(CN{ccnd})>10,
        end % for ccnd=1:length(SIG.BTRGLABEL),
        
    end % for LAGTYP = 1:length(dts),
    
    
    %     RSG.savedirplot = 'M:\BehaviourImaging\rsbs_workflow\rsbs_trigger_model\';
    %     if ~isdir(RSG.savedirplot), mkdir(RSG.savedirplot); end
    %
    %     interactive=0;toplot=1;
    %     ext ='-dtiff';fext    ='.tiff';
    %     %ext='-dpsc';fext='.ps';
    %
    %     if interactive, vis='on'; else vis='off'; end
    %
    %     for LAGTYP = 1:length(CORLAG),
    %         for ccnd=1:length(CORLAG{LAGTYP}.dat),
    %             if ~isempty(CORLAG{LAGTYP}.dat{ccnd}),
    %                 filename1 = sprintf('%sER%s%d_%s_B%d%s',RSG.savedirplot,SIG.BTRGLABEL{ccnd},LAGTYP,idinf{NID}.id,idinf{NID}.block,fext);
    %                 filename2 = sprintf('%sST%s%d_%s_B%d%s',RSG.savedirplot,SIG.BTRGLABEL{ccnd},LAGTYP,idinf{NID}.id,idinf{NID}.block,fext);
    %                 if toplot, delete(filename1); end
    %                 if toplot, delete(filename2); end
    %
    %                 for MODEL=1:length(CORLAG{LAGTYP}.dat{ccnd}),
    %                     f1=figure('Units','centimeters','Position',[2,2,16,30]);
    %                     f2=figure('Units','centimeters','Position',[2,2,40,60]);
    %                     cix1 = 1;
    %                     cix2 = 1;
    %                     for MODELPREDICTORNUM=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}), % MODELPREDICTORNUM=2 MODELPREDICTORNUM=3
    %                         for MODELPREDICTOR=1:length(CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}),
    %
    %                             %ccnd=1;MODEL=1;MODELPREDICTORNUM=2;MODELPREDICTOR=1;
    %
    %                             i=1;
    %                             CURDEP          = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURDEP;
    %                             CURPREDICTORSET = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.CURPREDICTORSET;
    %                             datYtst      = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ytst{i};
    %                             datYpredrstr = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ypredrstr{i};
    %                             datYpredfull = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ypredfull{i};
    %                             Ntrial       = CORLAG{LAGTYP}.dat{ccnd}{MODEL}{MODELPREDICTORNUM}{MODELPREDICTOR}.Ntrial;
    %                             t            = CORLAG{LAGTYP}.t;
    %                             rngix        = CORLAG{LAGTYP}.rngix;
    %
    %                             Ytest = datYtst(:);
    %                             Ypred = datYpredrstr(:);
    %                             SSE = nansum((Ytest-Ypred).^2); SST = nansum((Ytest-repmat(nanmean(Ytest),[size(Ytest,1),1])).^2); rsquare = 1-SSE./SST;
    %                             RSrstr = rsquare;
    %                             Ytest = datYtst(:);
    %                             Ypred = datYpredfull(:);
    %                             SSE = nansum((Ytest-Ypred).^2); SST = nansum((Ytest-repmat(nanmean(Ytest),[size(Ytest,1),1])).^2); rsquare = 1-SSE./SST;
    %                             RSfull = rsquare;
    %
    %                             Ytest = datYtst';
    %                             Ypred = datYpredrstr';
    %                             SSE = nansum((Ytest-Ypred).^2); SST = nansum((Ytest-repmat(nanmean(Ytest),[size(Ytest,1),1])).^2); rsquare = 1-SSE./SST
    %                             RSrstr_tim = rsquare;
    %                             Ytest = datYtst';
    %                             Ypred = datYpredfull';
    %                             SSE = nansum((Ytest-Ypred).^2); SST = nansum((Ytest-repmat(nanmean(Ytest),[size(Ytest,1),1])).^2); rsquare = 1-SSE./SST
    %                             RSfull_tim = rsquare;
    %
    %                             figure(f1);
    %                             subplot(7,2,cix1);plot(tim,rsquare);cix1=cix1+1;
    %                             plot(tim,nanmean(datYtst,2),'k');
    %                             hold on;plot(tim,nanmean(datYpredrstr,2),'bx-','LineWidth',1);
    %                             hold on;plot(tim,nanmean(datYpredfull,2),'ro-','LineWidth',1);
    %                             legend({sprintf('rstr%1.2f',RSrstr),sprintf('full%1.2f',RSfull)},'Location','Best');
    %                             subplot(7,2,cix1);plot(tim,RSrstr_tim,'b');cix1=cix1+1;
    %                             hold on;plot(tim,RSfull_tim,'r');
    %
    %                             str=[];
    %                             for i=1:length(CURPREDICTORSET),
    %                                 str=[str,SIGPOP.label{CURPREDICTORSET(i)},','];
    %                             end
    %                             str(end)=[];
    %                             title(sprintf('%s: %s[%s]',SIG.BTRGLABEL{ccnd},SIGPOP.label{CURDEP},str));
    %
    %                             if 1, % single trials
    %                                 figure(f2);
    %                                 subplot(4,3,cix2);cix2=cix2+1; % figure;
    %                                 tax=t;ix0=find(tax==0);
    %                                 ix=[1:(length(t)):(length(t))*Ntrial(2)]+ix0-1;
    %                                 tax = [1:(length(t))*Ntrial(2)];
    %                                 dum=datYtst;dum=dum(:)'; %dum(isnan(dum))=0;
    %                                 hold on;plot(tax,dum,'k','LineWidth',1);
    %                                 dum=datYpredrstr;
    %                                 hold on;plot(tax,dum(:)','bx-','LineWidth',1);
    %                                 dum=datYpredfull;
    %                                 hold on;plot(tax,dum(:)','ro-','LineWidth',1);
    %                                 hold on;plot([ix',ix']',repmat(get(gca,'YLim'),[length(ix),1])');
    %                                 legend({sprintf('rstr%1.2f',RSrstr),sprintf('full%1.2f',RSfull)},'Location','Best');
    %                                 axis tight;
    %                                 str=[];
    %                                 for i=1:length(CURPREDICTORSET),
    %                                     str=[str,SIGPOP.label{CURPREDICTORSET(i)},','];
    %                                 end
    %                                 str(end)=[];
    %                                 title(sprintf('%s: %s[%s]',SIG.BTRGLABEL{ccnd},SIGPOP.label{CURDEP},str));
    %                             end
    %                         end % for MODELPREDICTOR
    %                     end % for MODELPREDICTORNUM=1
    %                     if interactive, pause; end
    %                     if toplot,
    %                         figure(f1);
    %                         set(gcf,'PaperpositionMode','auto');
    %                         if toplot, eval(sprintf('print -f%d -append %s %s',gcf,ext,filename1)); end
    %                         close;
    %
    %                         figure(f2);
    %                         set(gcf,'PaperpositionMode','auto');
    %                         if toplot, eval(sprintf('print -f%d -append %s %s',gcf,ext,filename2)); end
    %                         close;
    %                     end
    %                 end %  for MODEL=1:length(DEP),
    %             end % if ~isempty(CORLAG{LAGTYP}.dat{ccnd}),
    %         end % for ccnd=1:length(CORLAG{LAGTYP}.dat),
    %     end % for LAGTYP = 1:length(dts),
else
    CORLAG=[];
end % if 1, % plotting


