function [DPOP]=rsbs_extract_data_trigger_alltrials_switching_AC(out,ADDINF,idinf,NID,EYE,EYEINF,SW,RIXSES)

%% This function extracts behavioural data to get triggers around certain events, then extracts calcium signals for different cell types around those events, then analyses population activity (AC)

[sel,BTRGLABEL,MOT,ODOURDELAY,gratingonnormal,gratingon2]=rsbs_extra_data_trigger_behavdata(out,ADDINF,idinf,NID,EYE,EYEINF); % This function extracts all behavioural triggers (AC)

[SIG]=rsbs_extra_data_trigger_VR_main(out,EYE, sel,BTRGLABEL,MOT);  % This script extracts calcium data around the extracted triggers (AC)

%[DPOP]=rsbs_extra_data_trigger_pop(SIG,idinf,NID,ODOURDELAY,LN,RIXSES,gratingonnormal,gratingon2);

[DPOP]=rsbs_extra_data_trigger_pop_select(SIG,idinf,NID,ODOURDELAY,SW,RIXSES,gratingonnormal,gratingon2); % Added on 13/07/2016 (AC)

% 
% % median split
% out.POOLLAYER=1;
% PARAM.prc=[30];
% clear SPLIT;
% [SPLIT.PD,SPLIT.PD_okchn,SPLIT.PD_POOLCV,SPLIT.PD_Twin]=rsbs_split_prctile_addtrg(out,sel,BTRGLABEL,PARAM);
% 
% %------
% % get single frame data for later stratification
% [STRATDAT]=rsbs_extra_data_trigger_stratdat(SIG,idinf,NID,ODOURDELAY,LN,RIXSES);
% 
% 
% % stratification for running speed in 1s window
% % edit rsbs_ROI_CRSP_rsp_sdprepst_pop.m
% % edit rsbs_ROI_CRSP_rsp_sdprepst_pop_run.m
% 
% % stratification for mean response
% 
% 
% % stratification for standard deviation