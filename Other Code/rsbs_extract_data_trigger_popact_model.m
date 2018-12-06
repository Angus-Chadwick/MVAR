function [TOT]=rsbs_extract_data_trigger_popact_model(out,ADDINF,idinf,NID,EYE,EYEINF,LN,RIXSES);

[sel,BTRGLABEL,MOT,ODOURDELAY,gratingonnormal,gratingon2]=rsbs_extra_data_trigger_behavdata(out,ADDINF,idinf,NID,EYE,EYEINF);

[SIG,SIGPOP]=rsbs_extra_data_trigger_VR_main_model(out,EYE, sel,BTRGLABEL,MOT,LN,RIXSES);

if 0,
    for chn=5:8
        sel=find(SIGPOP.sel.cond==1);
        figure;plot(SIGPOP.tax,squeeze(nanmean(SIGPOP.DAT{1}(sel,chn,:),1)),'k');
        title(sprintf('%s',SIGPOP.label{chn}));
        pause;close;
    end
end

%--------------------------------------------------------------------------
% 1. rsbs_placecorr (done at other place: rsbs_extract_data_trigger_popact)
% % if 0,
%     %rsbs_popactivity
%     % rsbs_extract_data_trigger_popact
%     
%       out.POOLCV  = DPOP{1}.POOLCV;
%     out.POOLDAT = DUM;
%     
%     PF=rsbs_placecorr(out,ADDINF,LN,RIXSES);
%      
%     % partial correlation
%     A=[1:4]; % 4 celltypes
%     out.PARTIAL.rng{1} = A;
%     out.PARTIAL.rng{2} = 4+A;
%     
%     TCOR=rsbs_placecorr_trigger(out,ADDINF,LN,RIXSES);
% end

%--------------------------------------------------------------------------
% 2. rsbs_predict_resp_corrlag
if 1,
    
    [CORLAG]=rsbs_extra_data_trigger_popact_model_corlag(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
    
    %[DPOP]=rsbs_extra_data_trigger_pop_predict_allsmp(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
    TOT.CORLAG=CORLAG;
    
end

%--------------------------------------------------------------------------
% 3. model
if 0,
    
    % edit rsbs_extra_data_trigger_pop_predict
    [MODDAT]=rsbs_extract_data_trigger_popact_model_predict(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
    TOT.MODDAT=MODDAT;
%     
%     %dbstop if error;
%     %[DPOP]=rsbs_extra_data_trigger_pop_predict(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
%     [DPOP]=rsbs_extra_data_trigger_pop_predict_simple(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
end