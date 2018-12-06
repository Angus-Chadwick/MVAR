function [ctot_rel, ctot_irrel, ctot_rel_rw, ctot_irrel_rw, ctot_rel_shuffres, ctot_irrel_shuffres, CELLLAB_TOT]  = model_noisecorr_switching(changeI, changeW)


fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;

%clear all;
ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/';
RSG.savedirforAngus = [ROOT,'Saved_Data/'];

if ~exist('TSK'), TSK=[]; end
if ~isfield(TSK,'SW'),
    TSK.SW=load([RSG.savedirforAngus,'SW_TLTPDPOP.mat']);
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
          
        clear ADDTRG;
        extract_ADDTRG=0;
        if extract_ADDTRG, % rsbs_extract_data_trigger
            dbstop if error;
            ADDTRG=rsbs_extract_data_trigger(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES);
            fname=sprintf('%sADDTRG2_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            save(fname,'ADDTRG','-v7.3');
        else
            fname=sprintf('%sADDTRG2_alltrials_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            load(fname,'ADDTRG');
        end
       
        
 
    %% Get cell labels

     POOLCV=[];
        for layer=1:4,
            dum=out.CRSPR{layer}.conversion(:,2);
            POOLCV  = cat(1,POOLCV,[repmat(layer,[size(dum,1),1]),dum]);
        end % for layer=1:length(TOT{ses}.SIG.DAT),
        
        % find cell type of each cell
        CELLLAB=NaN(size(POOLCV,1),1); % index into TL
        for chn=1:size(POOLCV,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCV(chn,1))&(TSK.(nam).TL.B==POOLCV(chn,2)) );  
            if not(isempty(ix)),
                CELLLAB(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end

                
CELLLAB(isnan(CELLLAB)) = [];

CELLLAB0 = CELLLAB;

CELLLAB_ALL{rix} = CELLLAB;

TOT.ADDTRG = ADDTRG;
TOT.RIXSES = RIXSES;


% Choose behavioural conditions to include for model fitting

    
       Conds_rel = [1,2,12];  % Vertical, angled, grey 
        
       Conds_irrel = [5,6,13]; % Vertical, angled, grey
       
       Nconds = length(Conds_rel);

for Cond=1:Nconds

    Condi_rel = Conds_rel(Cond);
    Condi_irrel = Conds_irrel(Cond);

    rmat_rel{rix, Cond}  = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.group1;
    rmat_irrel{rix, Cond} = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_irrel}.group1;

    drmat_rel{rix, Cond} = diff(rmat_rel{rix, Cond},1,3);
    drmat_irrel{rix, Cond} = diff(rmat_irrel{rix, Cond},1,3);

    vmatrel{rix, Cond} = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.group2;
    vmatirrel{rix, Cond} = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_irrel}.group2;

tsamples = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.t;
rngSTM = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.rngSTM;
rngBSL = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.rngBSL;
rngtot = [rngBSL, rngSTM];
Ntrialsrel{Cond} = size(drmat_rel{rix, Cond},1);
Ntrialsirrel{Cond} = size(drmat_irrel{rix, Cond},1);
Nsamples = length(rngtot);
tMOT{Cond} = TOT.ADDTRG{1}.PL{1}.DCOR{Condi_rel}.tMOT;

clear Nvrel
clear Nvirrel

for Trl=1:Ntrialsrel{Cond}

    vmatrel0{rix, Cond}(Trl,:)  = interp1(tMOT{Cond}, vmatrel{rix, Cond}(Trl,:), tsamples);
    Nvrel(Trl) = sum(isnan(vmatrel0{rix, Cond}(Trl,:)));

    if Nvrel(Trl) > 0
        
        V = vmatrel0{rix, Cond}(Trl,:);
        
        vmatrel0{rix, Cond}(Trl, isnan(V)) = vmatrel0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
    
end

for Trl=1:Ntrialsirrel{Cond}
    
    vmatirrel0{rix, Cond}(Trl,:) = interp1(tMOT{Cond}, vmatirrel{rix, Cond}(Trl,:), tsamples);
    Nvirrel(Trl) = sum(isnan(vmatirrel0{rix, Cond}(Trl,:)));
    
    if Nvirrel(Trl) > 0
        
        V = vmatirrel0{rix, Cond}(Trl,:);
        
        vmatirrel0{rix, Cond}(Trl, isnan(V)) = vmatirrel0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
        
end

if instantaneous

    drmatTOT_rel{Cond} = drmat_rel{rix, Cond}(:,:,rngtot-1); % or put -1 in this one?
    rmatTOT_rel{Cond} = rmat_rel{rix, Cond}(:,:,rngtot);
    vmatTOT_rel{Cond} = vmatrel0{rix, Cond}(:, rngtot);

    drmatTOT_irrel{Cond} = drmat_irrel{rix, Cond}(:,:,rngtot-1);
    rmatTOT_irrel{Cond} = rmat_irrel{rix, Cond}(:,:,rngtot);
    vmatTOT_irrel{Cond} = vmatirrel0{rix, Cond}(:, rngtot);

elseif instantaneous_vonly
    
       
    drmatTOT_rel{Cond} = drmat_rel{rix, Cond}(:,:,rngtot);
    rmatTOT_rel{Cond} = rmat_rel{rix, Cond}(:,:,rngtot);
    vmatTOT_rel{Cond} = vmatrel0{rix, Cond}(:, rngtot+1);

    drmatTOT_irrel{Cond} = drmat_irrel{rix, Cond}(:,:,rngtot);
    rmatTOT_irrel{Cond} = rmat_irrel{rix, Cond}(:,:,rngtot);
    vmatTOT_irrel{Cond} = vmatirrel0{rix, Cond}(:, rngtot+1);

    
    
else
    
    drmatTOT_rel{Cond} = drmat_rel{rix, Cond}(:,:,rngtot);
    rmatTOT_rel{Cond} = rmat_rel{rix, Cond}(:,:,rngtot);
    vmatTOT_rel{Cond} = vmatrel0{rix, Cond}(:, rngtot);

    drmatTOT_irrel{Cond} = drmat_irrel{rix, Cond}(:,:,rngtot);
    rmatTOT_irrel{Cond} = rmat_irrel{rix, Cond}(:,:,rngtot);
    vmatTOT_irrel{Cond} = vmatirrel0{rix, Cond}(:, rngtot);

    
end
    
clear Nrel
clear Nirrel

for i=1:Ntrialsrel{Cond}

    M = squeeze(drmatTOT_rel{Cond}(i,:,:));
    Nrel(i) = sum(sum(isnan(M)));

end

for i=1:Ntrialsirrel{Cond}

    M = squeeze(drmatTOT_irrel{Cond}(i,:,:));
    Nirrel(i) = sum(sum(isnan(M)));

end

Trialsrel{Cond} = find(Nrel == 0);
Trialsirrel{Cond} = find(Nirrel == 0);

drmatTOT_rel{Cond} = drmatTOT_rel{Cond}(Trialsrel{Cond},:,:);
rmatTOT_rel{Cond} = rmatTOT_rel{Cond}(Trialsrel{Cond},:,:);
vmatTOT_rel{Cond} = vmatTOT_rel{Cond}(Trialsrel{Cond},:);

drmatTOT_irrel{Cond} = drmatTOT_irrel{Cond}(Trialsirrel{Cond},:,:);
rmatTOT_irrel{Cond} = rmatTOT_irrel{Cond}(Trialsirrel{Cond},:,:);
vmatTOT_irrel{Cond} = vmatTOT_irrel{Cond}(Trialsirrel{Cond},:);

Ntrialsrel{Cond} = length(Trialsrel{Cond});
Ntrialsirrel{Cond} = length(Trialsirrel{Cond});

end



%% Fit full (V & A, BSL & STM) simultaneously

clear rmat0rel
clear drmat0rel
clear rmat0irrel
clear drmat0irrel


for Cond = 1:Nconds

    drmatSTM_rel = drmatTOT_rel{Cond};
    rmatSTM_rel = rmatTOT_rel{Cond};
    vmatSTM_rel = vmatTOT_rel{Cond};

    drmatSTM_irrel = drmatTOT_irrel{Cond};
    rmatSTM_irrel = rmatTOT_irrel{Cond};
    vmatSTM_irrel = vmatTOT_irrel{Cond};

    rmatSTM_irrel_trialmean{rix, Cond} = squeeze(mean(rmatSTM_irrel));
    rmatSTM_rel_trialmean{rix, Cond} = squeeze(mean(rmatSTM_rel));

    M = permute(drmatSTM_rel, [3,1,2]);
    drmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_rel, [3,1,2]);
    rmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]
    M = permute(vmatSTM_rel, [2, 1]);
    vmatSTM_rel0{Cond} = reshape(M, [size(M,1) * size(M,2),1]);

    M = permute(drmatSTM_irrel, [3,1,2]);
    drmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_irrel, [3,1,2]);
    rmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(vmatSTM_irrel, [2, 1]);
    vmatSTM_irrel0{Cond} = reshape(M, [size(M,1) * size(M,2), 1]);

    vmatrel_all{rix,Cond}   = mean(vmatTOT_rel{Cond}); 
    vmatirrel_all{rix, Cond} = mean(vmatTOT_irrel{Cond}); 
    
    vmatrel0{rix, Cond} = vmatTOT_rel{Cond};
    vmatirrel0{rix, Cond} = vmatTOT_irrel{Cond};

    CLS = [1,3,4,5];


     rmat0rel{Cond}= rmatSTM_rel0.';
     drmat0rel{Cond} = drmatSTM_rel0.';
 
     rmat0irrel{Cond} = rmatSTM_irrel0.';
     drmat0irrel{Cond} = drmatSTM_irrel0.';


end



%% Set up regression matrix equations


[xrel{rix}, xirrel{rix}, Measurement_rel{rix}, Measurement_irrel{rix}, Prediction_rel{rix}, Prediction_irrel{rix}] = fitLDS_Switching(fixvars, changeW, changeI, drmat0rel, drmat0irrel, rmat0rel, rmat0irrel, rngtot, Ntrialsrel, Ntrialsirrel, vmatSTM_rel0, vmatSTM_irrel0, CELLLAB_ALL{rix});

residualrel{rix}  = Measurement_rel{rix} - Prediction_rel{rix};
residualirrel{rix} = Measurement_irrel{rix} - Prediction_irrel{rix};

%% Calculate pooled variance and selectivity of model fits


if and(ismember(1, Conds_rel), ismember(2,Conds_rel))  % get selectivity measures
    
    for T=1:Ntrialsrel{1}

      ResidualV_rel{rix}(:,T) = mean(squeeze(drmatTOT_rel{1}(T,:,9:end)) - xrel{rix}(:, [1:size(xrel{rix},1), (size(xrel{rix},1) + 9):(size(xrel{rix},1) + 17), end]) * [squeeze(rmatTOT_rel{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_rel{1}(T,9:end))], 2);
      ResidualV_rel_TOT{rix}(:,:,T) = squeeze(drmatTOT_rel{1}(T,:,:))    - xrel{rix}(:, [1:size(xrel{rix},1), (size(xrel{rix},1) + 1):(size(xrel{rix},1) + 17), end]) * [squeeze(rmatTOT_rel{1}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_rel{1}(T,:))];
      
    end

    for T=1:Ntrialsrel{2}

        ResidualA_rel{rix}(:,T) = mean(squeeze(drmatTOT_rel{2}(T,:,9:end)) - xrel{rix}(:, [1:size(xrel{rix},1), (size(xrel{rix},1) + 9 + 17):(size(xrel{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_rel{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_rel{2}(T,9:end))], 2);
        ResidualA_rel_TOT{rix}(:,:,T) = squeeze(drmatTOT_rel{2}(T,:,:))    - xrel{rix}(:, [1:size(xrel{rix},1), (size(xrel{rix},1) + 1 + 17):(size(xrel{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_rel{2}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_rel{2}(T,1:end))];

    end

    for T=1:Ntrialsirrel{1}

      ResidualV_irrel{rix}(:,T) = mean(squeeze(drmatTOT_irrel{1}(T,:,9:end)) - xirrel{rix}(:, [1:size(xirrel{rix},1), (size(xirrel{rix},1) + 9):(size(xirrel{rix},1) + 17), end]) * [squeeze(rmatTOT_irrel{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_irrel{1}(T,9:end))], 2);
      ResidualV_irrel_TOT{rix}(:,:,T) = squeeze(drmatTOT_irrel{1}(T,:,:))    - xirrel{rix}(:, [1:size(xirrel{rix},1), (size(xirrel{rix},1) + 1):(size(xirrel{rix},1) + 17), end]) * [squeeze(rmatTOT_irrel{1}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_irrel{1}(T,1:end))];
     
    end

    for T=1:Ntrialsirrel{2}

      ResidualA_irrel{rix}(:,T) = mean(squeeze(drmatTOT_irrel{2}(T,:,9:end)) - xirrel{rix}(:, [1:size(xirrel{rix},1), (size(xirrel{rix},1) + 9 + 17):(size(xirrel{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_irrel{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_irrel{2}(T,9:end))], 2);
      ResidualA_irrel_TOT{rix}(:,:,T) = squeeze(drmatTOT_irrel{2}(T,:,:))    - xirrel{rix}(:, [1:size(xirrel{rix},1), (size(xirrel{rix},1) + 1 + 17):(size(xirrel{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_irrel{2}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_irrel{2}(T,1:end))];
      
    end

      nrelA = size(ResidualA_rel{rix},2);
      nrelV = size(ResidualV_rel{rix},2);
      nirrelA = size(ResidualA_irrel{rix},2);
      nirrelV = size(ResidualV_irrel{rix},2);

      PoolVar_rel_Out{rix}  = sqrt(( var(mean(rmatTOT_rel{1}(:,:,9:end),3))  * (nrelV-1 ) + var(mean(rmatTOT_rel{2}(:,:,9:end),3)) * (nrelA-1) ) / ( (nrelV - 1) + (nrelA - 1)));
      PoolVar_irrel_Out{rix} = sqrt(( var(mean(rmatTOT_irrel{1}(:,:,9:end),3)) * (nirrelV-1) + var(mean(rmatTOT_irrel{2}(:,:,9:end),3)) * (nirrelA-1) ) / ( (nirrelV - 1) + (nirrelA - 1)));

     
end

% find selectively responding cells

    for i=1:size(rmat_rel{rix, 1}, 2)

        psel_irrel{rix}(i) = ranksum(nanmean(rmat_irrel{rix, 1}(:,i,rngtot(9:end)), 3), nanmean(rmat_irrel{rix,2}(:,i,rngtot(9:end)),3));
        psel_rel{rix}(i) = ranksum(nanmean(rmat_rel{rix, 1}(:,i,rngtot(9:end)), 3), nanmean(rmat_rel{rix,2}(:,i,rngtot(9:end)),3));

    end


Trialsrel_all{rix} = Trialsrel;
Trialsirrel_all{rix} = Trialsirrel;

drmatrel_all{rix} = horzcat(drmat0rel{:});
drmatirrel_all{rix} = horzcat(drmat0irrel{:});

rmatrel_all{rix} = horzcat(rmat0rel{:});
rmatirrel_all{rix} = horzcat(rmat0irrel{:});


end


CLS = [1,3,4,5];
CELLLAB_TOT = vertcat(CELLLAB_ALL{:});
RXS = 1:length(CELLLAB_ALL);


%% Analytically calculate the contribution to noise correlations arising from each source

tsims = 9:17;
% 
 [ctot_rel, ctot_irrel]       = noisecorr_analytic_decomp_relirrelshuff(xrel,    xirrel,    rmat_rel, rmat_irrel, ResidualV_rel_TOT, ResidualV_irrel_TOT, vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot, CELLLAB_ALL,RXS);

% Delete weights

for rix = 1:length(RXS)
    
    xrel_rw{rix} = diag(diag(xrel{RXS(rix)}(:,1:size(xrel{RXS(rix)},1))));
    xirrel_rw{rix} = diag(diag(xirrel{RXS(rix)}(:,1:size(xirrel{RXS(rix)},1))));

end
     
 [ctot_rel_rw, ctot_irrel_rw] = noisecorr_analytic_decomp_relirrelshuff(xrel_rw, xirrel_rw, rmat_rel, rmat_irrel, ResidualV_rel_TOT, ResidualV_irrel_TOT, vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot, CELLLAB_ALL, RXS);

% 
shufflemethod = 'trials';

RXS = 1:length(CELLLAB_ALL);

Nshuff = 100;

for S= 1:Nshuff

for rix = RXS

LrelV = size(ResidualV_rel_TOT{rix}, 3); LirrelV = size(ResidualV_irrel_TOT{rix}, 3);
for L = 1:(LrelV)
ResidualV_TOT{rix}(:,:,L) = ResidualV_rel_TOT{rix}(:,:,L);
end
for L = 1:(LirrelV)
ResidualV_TOT{rix}(:,:,LrelV + L) = ResidualV_irrel_TOT{rix}(:,:,L);
end

if strmatch(shufflemethod, 'full')

    for i=1:size(ResidualV_TOT{rix},1)

        ResidualV_TOT_shuff{rix}(i,:,:) = ResidualV_TOT{rix}(i,:,randperm(LrelV + LirrelV));

    end
                 
elseif strmatch(shufflemethod, 'trials')
    
    ResidualV_TOT_shuff{rix} = ResidualV_TOT{rix}(:,:,randperm(LrelV + LirrelV));

end

ResidualV_rel_TOT_shuff{rix} = ResidualV_TOT_shuff{rix}(:,:, 1:size(ResidualV_rel_TOT{rix}, 3));
ResidualV_irrel_TOT_shuff{rix} = ResidualV_TOT_shuff{rix}(:,:, (size(ResidualV_rel_TOT{rix}, 3)+1):end);
% 
end

[ctot_rel_shuffres{S}, ctot_irrel_shuffres{S}] = noisecorr_analytic_decomp_relirrelshuff(xrel, xirrel, rmat_rel, rmat_irrel, ResidualV_rel_TOT_shuff, ResidualV_irrel_TOT_shuff, vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot, CELLLAB_ALL,RXS);

end



