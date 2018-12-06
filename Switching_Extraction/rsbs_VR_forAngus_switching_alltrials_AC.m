clear all;
% cd('C:\Jasper\OneDrive\OneDrive - University College London\forAngus\Angus_Switching_Code\')
location = 1;
if location == 2,
    % otherwise problem with nansum nanstd etc
    rmpath(genpath('C:\Jasper\Matlab\fieldtrip-20120821\')); 
end

if location == 1,
    ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/';
    RSG.savedirforAngus = [ROOT,'Saved_Data/'];

    if ~exist('TSK'), TSK=[]; end
    if ~isfield(TSK,'SW'),
        TSK.SW=load([RSG.savedirforAngus,'SW_TLTPDPOP.mat']);
    end
end

% INCLUDE_SD = { % rsbs_sd_pop_prepro_cutoffs
%     'M70_20141106_B1'
%     'M73_20141101_B1'
%     'M75_20141102_B1'
%     'M75_20141107_B1'
%     'M80_20141031_B1'
%     'M81_20141031_B1'
%     'M87_20141108_B1'
%     'M89_20141030_B1'
%     'M93_20141111_B1'
%     } % M93_20141023_B1


INCLUDE_SD = { % rsbs_sd_pop_prepro_cutoffs
'M70_20141104_B1'
'M70_20141106_B1'
'M70_20141115_B1'
'M71_20141104_B1'
'M73_20141101_B1'
'M73_20141104_B1'
'M75_20141102_B1'
'M75_20141105_B1'
'M75_20141107_B1'
'M75_20141117_B1'
'M80_20141031_B1'
'M80_20141103_B2'
'M80_20141108_B1'
'M80_20141114_B1'
'M81_20141031_B1'
'M81_20141105_B1'
'M81_20141108_B1'
'M81_20141113_B1'
'M81_20141117_B1'
'M87_20141105_B1'
'M87_20141108_B1'
'M87_20141110_B1'
'M89_20141030_B1'
'M89_20141103_B1'
'M89_20141113_B1'
'M89_20141115_B1'
'M93_20141103_B1'
'M93_20141107_B1'
'M93_20141111_B1'

 } 


for rix=1:size(INCLUDE_SD,1),
    
    
    clear TOT;
    
    clear idinf;NID = 1;
    name = INCLUDE_SD{rix}
    ix1=strfind(name,'_B');
    
    idinf{NID}.id    = name(1:ix1-1);
    idinf{NID}.block = str2num(name(ix1+2:end));
    
    clear EXTINF;
    fname=sprintf('%s%s_B%d_extinf',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'EXTINF');
    nam    = EXTINF.nam;
    if strcmp(nam,'LN'), idinf{NID}.type='SD'; elseif strcmp(nam,'SW'), idinf{NID}.type='SWITCH'; end
    
    RIXSES = EXTINF.RIXSES;
    EYE    = EXTINF.EYE;
    EYEINF = EXTINF.EYEINF;
    
    clear out ADDINF;
    fname=sprintf('%s%s_B%d_dat',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'out','ADDINF');
    
    if 0,
        ADDINF.sel.cond==1       
        ADDINF.gratingonnormal'
        
        tmp=TSK.SW.TOTALLredux{EXTINF.RIXSES(1)}.TRL{RIXSES(2)};
        tmp.sel.trl(tmp.CN{1})'
    end
    %------------------------------------------------------------------
    % run analysis
    
    clear ADDTRG;
    extract_ADDTRG=1;
    if extract_ADDTRG, % rsbs_extract_data_trigger
        dbstop if error;
        ADDTRG=rsbs_extract_data_trigger_alltrials_switching_AC(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES);
        fname=sprintf('%sADDTRG2_alltrials_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
        save(fname,'ADDTRG','-v7.3');
    else
        fname=sprintf('%sADDTRG2_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
        load(fname,'ADDTRG');
    end
    
    % save data needed for later
    TOT.RIXSES = RIXSES;
    TOT.ADDTRG = ADDTRG;
    
end % for NID=1:length(idinf),
