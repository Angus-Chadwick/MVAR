function [TOT]=rsbs_extract_data_trigger_cell_model(out,ADDINF,idinf,NID,EYE,EYEINF,LN,RIXSES,ADDTRG)

[sel,BTRGLABEL,MOT,ODOURDELAY,gratingonnormal,gratingon2]=rsbs_extra_data_trigger_behavdata(out,ADDINF,idinf,NID,EYE,EYEINF);

[SIG,SIGPOP]=rsbs_extra_data_trigger_VR_main_model(out,EYE, sel,BTRGLABEL,MOT,LN,RIXSES);

if 0,
    for chn=5:8
        sel=find(SIGPOP.sel.cond==1);
        figure;plot(SIGPOP.tax,squeeze(nanmean(SIGPOP.DAT{1}(sel,chn,:),1)),'k');
        title(sprintf('%s',SIGPOP.label{chn}));
        pause;close;
    end



%--------------------------------------------------------------------------
% 1. rsbs_placecorr (see also rsbs_extract_data_trigger_popact)

% if 0,
% see rsbs_popactivity
out.POOLLAYER=1;
out.POOLCV=[];
for layer=1:4,
    dum=out.CRSPR{layer}.conversion(:,2);
    out.POOLCV  = cat(1,out.POOLCV,[repmat(layer,[size(dum,1),1]),dum]);
end % for layer=1:length(TOT{ses}.SIG.DAT),

out.POOLDAT=[];
for layer=1:4,
    DUM=interp1(out.frmtim(out.sel{layer}),out.dFj{1}{layer}',out.frmtim);
    out.POOLDAT=cat(1,out.POOLDAT,DUM');
end
TCOR=rsbs_placecorr(out,ADDINF,LN,RIXSES);

end
%--------------------------------------------------------------------------
% 2. rsbs_predict_resp_corrlag
if 1,
    
    [CORLAG]=rsbs_extra_data_trigger_cell_model_corlag(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2,ADDTRG);
    %[CORRLAG]=rsbs_extra_data_trigger_popact_model_corlag(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
    
    %[DPOP]=rsbs_extra_data_trigger_pop_predict_allsmp(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
    TOT.CORLAG=CORLAG;
end

%--------------------------------------------------------------------------
% 3. model

% % if 0,
% % how much does population activity of one cell type explain of the other cell type? (in addition to its own history)
% % rsbs_random_forest_stats
% % rsbs_predict_resp_corrlag
% % rsbs_predict_resp_granger

% 
% [MOD]=rsbs_extra_data_trigger_cell_model_corlag(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2,ADDTRG);
% TOT.MOD=MOD;
% 
% rsbs_extra_data_trigger_pop_predictcell_simple
% %[DPOP]=rsbs_extra_data_trigger_pop_predict_allsmp(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);
% 
% %dbstop if error; rsbs_extra_data_trigger_pop_predict
% %[DPOP]=rsbs_extra_data_trigger_pop_predictcell_simple(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2,ADDTRG);
% [DPOP]=rsbs_extra_data_trigger_pop_predictcell_simple(SIG,SIGPOP,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2,ADDTRG);
% % end
