function [SI_model_pre_full, SI_model_post_full, SI_model_pre, SI_model_post, Model_pre, Model_post, Model_pre_full, Model_post_full, CELLLAB_TOT] = model_SI(changeI, changeW)
%function [SI, poolvar, CELLLAB_TOT] = model_SI(changeI, changeW)
%% Equations to do regression on linear dynamical system from individual trials

fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;
zscore_data = 0;

ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/Learning/';
RSG.savedirforAngus = [ROOT,'Saved_Data/'];

if ~exist('TSK'), TSK=[]; end
if ~isfield(TSK,'LN'),
    TSK.LN=load([RSG.savedirforAngus,'LN_TLTPDPOP.mat']);
end

INCLUDE_SD = { % rsbs_sd_pop_prepro_cutoffs
    'M70_20141022_B1'	'M70_20141028_B1' % M70_20141022_B1
    'M71_20141020_B1'	'M71_20141029_B1' % M71_20141020_B1 NB small number of trials
    'M72_20141020_B1'	'M72_20141101_B1' % M72_20141020_B1
    'M73_20141020_B1'	'M73_20141028_B1' % M73_20141020_B1
    'M75_20141021_B1'	'M75_20141029_B1' % M75_20141021_B1
    'M80_20141023_B1'	'M80_20141028_B1' % M80_20141023_B1
    'M81_20141021_B1'	'M81_20141029_B1' % M81_20141021_B1 NB small number of trials
    'M87_20141021_B1'	'M87_20141028_B1' % M87_20141021_B1
    'M89_20141021_B1'	'M89_20141029_B1' % M89_20141021_B1
    'M93_20141023_B1'	'M93_20141028_B1'} % M93_20141023_B1


for rix=1:size(INCLUDE_SD,1),
    
    if rix ~= 2 & rix ~= 7
    
    clear TOT;
    for PREPST=1:2,
        clear idinf;NID = 1;
        name = INCLUDE_SD{rix,PREPST}
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
        
        % save data needed for later
        TOT{PREPST}.RIXSES = RIXSES;
        TOT{PREPST}.ADDTRG = ADDTRG;        
    end % for PREPST=1:2,
   
    
    
            % find same cells pre and post
        CELL{1} = TOT{1}.ADDTRG{1}.POOLCV; % PRE
        CELL{2} = TOT{2}.ADDTRG{1}.POOLCV; % PST
        CELLUNI = unique([CELL{1};CELL{2}],'rows');
        IX=NaN(size(CELLUNI,1),2);
        for i=1:size(CELLUNI,1),
            for PREPST=1:2,
                ix=find(CELL{PREPST}(:,1)==CELLUNI(i,1)&CELL{PREPST}(:,2)==CELLUNI(i,2));
                if ~isempty(ix),
                    IX(i,PREPST)=ix;
                end
            end
        end
        IX = IX(all(~isnan(IX),2),:); % all cells with pre and post data
        
    

%% Get cell labels

 POOLCVPRE = CELL{1};
 POOLCVPOST = CELL{2};
 
        % find cell type of each cell
        CELLLABPRE=NaN(size(POOLCVPRE,1),1); % index into TL
        RIXSES = TOT{1}.RIXSES;
        for chn=1:size(POOLCVPRE,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCVPRE(chn,1))&(TSK.(nam).TL.B==POOLCVPRE(chn,2)) );  
            if not(isempty(ix)),
                CELLLABPRE(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end

                % find cell type of each cell
        CELLLABPOST=NaN(size(POOLCVPOST,1),1); % index into TL
                RIXSES = TOT{2}.RIXSES;
        for chn=1:size(POOLCVPOST,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCVPOST(chn,1))&(TSK.(nam).TL.B==POOLCVPOST(chn,2)) );  
            if not(isempty(ix)),
                CELLLABPOST(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end

        CELLLAB = CELLLABPRE(IX(:,1)); % common cell labels
        
        CELLLAB_ALL{rix} = CELLLAB;
                
% Inputs:
% rmat: matrix of time samples x trials x cells for the dF/F signal
% drmat: matrix of time samples x trials x cells for the change in dF/F between consecutive time samples
% Nsamples: number of time samples
% Ntrials: number of trials

% Choose behavioural conditions to include for model fitting

%Conds = 7; % grey onset only
%Conds = [1,2]; % gratings only
Conds = [1,2,7];  % Vertical, angled and grey onset
Nconds = length(Conds);


for Cond=1:Nconds

    Condi = Conds(Cond);

rmat_pre{rix, Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.group1(:, IX(:,1), :);
rmat_post{rix, Cond} = TOT{2}.ADDTRG{1}.PL{1}.DCOR{Condi}.group1(:, IX(:,2), :);

drmat_pre{rix, Cond} = diff(rmat_pre{rix, Cond},1,3);
drmat_post{rix, Cond} = diff(rmat_post{rix, Cond},1,3);

vmatpre{rix, Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.group2;
vmatpost{rix, Cond} = TOT{2}.ADDTRG{1}.PL{1}.DCOR{Condi}.group2;

% Outputs:
% Inputsmat: matrix of time samples x cells for the average input to each cell
% CONmat: matrix of cells x cells for the connectivity between cells

% Express dynamical systems equation as a matrix product

tsamples = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.t;
rngSTM = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.rngSTM;
rngBSL = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.rngBSL;
rngtot = [rngBSL, rngSTM];
Ntrialspre{Cond} = size(drmat_pre{rix, Cond},1);
Ntrialspost{Cond} = size(drmat_post{rix, Cond},1);
Nsamples = length(rngtot);
tMOT{Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.tMOT;

clear Nvpre
clear Nvpost

for Trl=1:Ntrialspre{Cond}

    vmatpre0{rix, Cond}(Trl,:)  = interp1(tMOT{Cond}, vmatpre{rix, Cond}(Trl,:), tsamples);
    Nvpre(Trl) = sum(isnan(vmatpre0{rix, Cond}(Trl,:)));

    if Nvpre(Trl) > 0
        
        V = vmatpre0{rix, Cond}(Trl,:);
        
        vmatpre0{rix, Cond}(Trl, isnan(V)) = vmatpre0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
    
end

for Trl=1:Ntrialspost{Cond}
    
    vmatpost0{rix, Cond}(Trl,:) = interp1(tMOT{Cond}, vmatpost{rix, Cond}(Trl,:), tsamples);
    Nvpost(Trl) = sum(isnan(vmatpost0{rix, Cond}(Trl,:)));
    
    if Nvpost(Trl) > 0
        
        V = vmatpost0{rix, Cond}(Trl,:);
        
        vmatpost0{rix, Cond}(Trl, isnan(V)) = vmatpost0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
        
end

if instantaneous

    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot-1); % or put -1 in this one?
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot-1);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot);

elseif instantaneous_vonly
    
       
    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot);
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot+1);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot+1);

    
    
else
    
    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot);
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot);

    
end
    
clear Npre
clear Npost

for i=1:Ntrialspre{Cond}

    M = squeeze(drmatTOT_pre{Cond}(i,:,:));
    Npre(i) = sum(sum(isnan(M)));

end

for i=1:Ntrialspost{Cond}

    M = squeeze(drmatTOT_post{Cond}(i,:,:));
    Npost(i) = sum(sum(isnan(M)));

end

Trialspre{Cond} = find(Npre == 0);
Trialspost{Cond} = find(Npost == 0);

drmatTOT_pre{Cond} = drmatTOT_pre{Cond}(Trialspre{Cond},:,:);
rmatTOT_pre{Cond} = rmatTOT_pre{Cond}(Trialspre{Cond},:,:);
vmatTOT_pre{Cond} = vmatTOT_pre{Cond}(Trialspre{Cond},:);

drmatTOT_post{Cond} = drmatTOT_post{Cond}(Trialspost{Cond},:,:);
rmatTOT_post{Cond} = rmatTOT_post{Cond}(Trialspost{Cond},:,:);
vmatTOT_post{Cond} = vmatTOT_post{Cond}(Trialspost{Cond},:);

Ntrialspre{Cond} = length(Trialspre{Cond});
Ntrialspost{Cond} = length(Trialspost{Cond});

end



%% Fit full (V & A, BSL & STM) simultaneously

clear rmat0pre
clear drmat0pre
clear rmat0post
clear drmat0post


for Cond = 1:Nconds

    drmatSTM_pre = drmatTOT_pre{Cond};
    rmatSTM_pre = rmatTOT_pre{Cond};
    vmatSTM_pre = vmatTOT_pre{Cond};

    drmatSTM_post = drmatTOT_post{Cond};
    rmatSTM_post = rmatTOT_post{Cond};
    vmatSTM_post = vmatTOT_post{Cond};

    rmatSTM_post_trialmean{rix, Cond} = squeeze(mean(rmatSTM_post));
    rmatSTM_pre_trialmean{rix, Cond} = squeeze(mean(rmatSTM_pre));

    M = permute(drmatSTM_pre, [3,1,2]);
    drmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_pre, [3,1,2]);
    rmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]
    M = permute(vmatSTM_pre, [2, 1]);
    vmatSTM_pre0{Cond} = reshape(M, [size(M,1) * size(M,2),1]);

    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(vmatSTM_post, [2, 1]);
    vmatSTM_post0{Cond} = reshape(M, [size(M,1) * size(M,2), 1]);

    vmatpre_all{rix,Cond}   = mean(vmatTOT_pre{Cond}); 
    vmatpost_all{rix, Cond} = mean(vmatTOT_post{Cond}); 
    
    vmatpre0{rix, Cond} = vmatTOT_pre{Cond};
    vmatpost0{rix, Cond} = vmatTOT_post{Cond};

    CLS = [1,3,4,5];


     rmat0pre{Cond}= rmatSTM_pre0.';
     drmat0pre{Cond} = drmatSTM_pre0.';
 
     rmat0post{Cond} = rmatSTM_post0.';
     drmat0post{Cond} = drmatSTM_post0.';


end


%% Set up regression matrix equations
       
[xpre{rix}, xpost{rix}, Measurement_pre{rix}, Measurement_post{rix}, Prediction_pre{rix}, Prediction_post{rix}] = fitLDS(fixvars, changeW, changeI, drmat0pre, drmat0post, rmat0pre, rmat0post, rngtot, Ntrialspre, Ntrialspost, vmatSTM_pre0, vmatSTM_post0, CELLLAB_ALL{rix});
      
residualpre{rix}  = Measurement_pre{rix} - Prediction_pre{rix};
residualpost{rix} = Measurement_post{rix} - Prediction_post{rix};
   
%% Calculate pooled variance and selectivity of model fits


if and(ismember(1, Conds), ismember(2,Conds))  % get selectivity measures
    
    for T=1:Ntrialspre{1}

      ResidualV_pre{rix}(:,T) = mean(squeeze(drmatTOT_pre{1}(T,:,9:end)) - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 9):(size(xpre{rix},1) + 17), end]) * [squeeze(rmatTOT_pre{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_pre{1}(T,9:end))], 2);
      ResidualV_pre_TOT{rix}(:,:,T) = squeeze(drmatTOT_pre{1}(T,:,:))    - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 1):(size(xpre{rix},1) + 17), end]) * [squeeze(rmatTOT_pre{1}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_pre{1}(T,:))];
      
    end

    for T=1:Ntrialspre{2}

        ResidualA_pre{rix}(:,T) = mean(squeeze(drmatTOT_pre{2}(T,:,9:end)) - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 9 + 17):(size(xpre{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_pre{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_pre{2}(T,9:end))], 2);
        ResidualA_pre_TOT{rix}(:,:,T) = squeeze(drmatTOT_pre{2}(T,:,:))    - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 1 + 17):(size(xpre{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_pre{2}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_pre{2}(T,1:end))];

    end

    for T=1:Ntrialspost{1}

      ResidualV_post{rix}(:,T) = mean(squeeze(drmatTOT_post{1}(T,:,9:end)) - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 9):(size(xpost{rix},1) + 17), end]) * [squeeze(rmatTOT_post{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_post{1}(T,9:end))], 2);
      ResidualV_post_TOT{rix}(:,:,T) = squeeze(drmatTOT_post{1}(T,:,:))    - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 1):(size(xpost{rix},1) + 17), end]) * [squeeze(rmatTOT_post{1}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_post{1}(T,1:end))];
     
    end

    for T=1:Ntrialspost{2}

      ResidualA_post{rix}(:,T) = mean(squeeze(drmatTOT_post{2}(T,:,9:end)) - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 9 + 17):(size(xpost{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_post{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_post{2}(T,9:end))], 2);
      ResidualA_post_TOT{rix}(:,:,T) = squeeze(drmatTOT_post{2}(T,:,:))    - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 1 + 17):(size(xpost{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_post{2}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_post{2}(T,1:end))];
      
    end

      npreA = size(ResidualA_pre{rix},2);
      npreV = size(ResidualV_pre{rix},2);
      npostA = size(ResidualA_post{rix},2);
      npostV = size(ResidualV_post{rix},2);

      PoolVar_pre_Out{rix}  = sqrt(( var(mean(rmatTOT_pre{1}(:,:,9:end),3))  * (npreV-1 ) + var(mean(rmatTOT_pre{2}(:,:,9:end),3)) * (npreA-1) ) / ( (npreV - 1) + (npreA - 1)));
      PoolVar_post_Out{rix} = sqrt(( var(mean(rmatTOT_post{1}(:,:,9:end),3)) * (npostV-1) + var(mean(rmatTOT_post{2}(:,:,9:end),3)) * (npostA-1) ) / ( (npostV - 1) + (npostA - 1)));

     
end
     
Trialspre_all{rix} = Trialspre;
Trialspost_all{rix} = Trialspost;

drmatpre_all{rix} = horzcat(drmat0pre{:});
drmatpost_all{rix} = horzcat(drmat0post{:});

    end
end

CELLLAB_TOT = vertcat(CELLLAB_ALL{:});

%% Analytically calculate the contribution to selectivity arising from each source

var_method = 'fullvar';
manipulation = 'removeall';
shufflemethod = 'trials'; % shuffle residuals for each cell independently ('full') or shuffle trials across residual vectors ('trials')

clear SI
clear poolvar

tsims = 9:17;

W1 = 1; W2 = 1;

% Calculate full model with residuals

[SI, poolvar] = selectivity_analytic(manipulation, var_method, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);

SI_model_pre_full = SI.SIinputs(:,1) + SI.SIinitstate(:,1) + SI.SIvel(:,1) + SI.SIres(:,1);
SI_model_post_full = SI.SIinputs(:,2) + SI.SIinitstate(:,2) + SI.SIvel(:,2) + SI.SIres(:,2);

% Shuffle residuals across learning

RXS = [1,3,4,5,6,8,9,10];

Nshuff = 100;

for S= 1:Nshuff

for rix = RXS

LpreA = size(ResidualA_pre_TOT{rix}, 3); LpostA = size(ResidualA_post_TOT{rix}, 3);
for L = 1:(LpreA)
ResidualA_TOT{rix}(:,:,L) = ResidualA_pre_TOT{rix}(:,:,L);
end
for L = 1:(LpostA)
ResidualA_TOT{rix}(:,:,LpreA + L) = ResidualA_post_TOT{rix}(:,:,L);
end
LpreV = size(ResidualV_pre_TOT{rix}, 3); LpostV = size(ResidualV_post_TOT{rix}, 3);
for L = 1:(LpreV)
ResidualV_TOT{rix}(:,:,L) = ResidualV_pre_TOT{rix}(:,:,L);
end
for L = 1:(LpostV)
ResidualV_TOT{rix}(:,:,LpreV + L) = ResidualV_post_TOT{rix}(:,:,L);
end

if strmatch(shufflemethod, 'full')

       % dum = find(CELLLAB_ALL{rix} == 4);  % choose cell type whose residuals are to be shuffled
    
        ResidualA_TOT_shuff{rix} = ResidualA_TOT{rix};
        ResidualV_TOT_shuff{rix} = ResidualV_TOT{rix};
    
   for i=1:size(ResidualV_TOT{rix},1) 

      % if ismember(i, dum)      
        
            ResidualA_TOT_shuff{rix}(i,:,:) = ResidualA_TOT{rix}(i,:,randperm(LpreA + LpostA));
            ResidualV_TOT_shuff{rix}(i,:,:) = ResidualV_TOT{rix}(i,:,randperm(LpreV + LpostV));

     % end
        
   end
                 
elseif strmatch(shufflemethod, 'trials')
    
    ResidualA_TOT_shuff{rix} = ResidualA_TOT{rix}(:,:,randperm(LpreA + LpostA));
    ResidualV_TOT_shuff{rix} = ResidualV_TOT{rix}(:,:,randperm(LpreV + LpostV));

end
    
ResidualA_pre_TOT_shuff{rix} = ResidualA_TOT_shuff{rix}(:,:, 1:size(ResidualA_pre_TOT{rix}, 3));
ResidualA_post_TOT_shuff{rix} = ResidualA_TOT_shuff{rix}(:,:, (size(ResidualA_pre_TOT{rix}, 3)+1):end);
ResidualV_pre_TOT_shuff{rix} = ResidualV_TOT_shuff{rix}(:,:, 1:size(ResidualV_pre_TOT{rix}, 3));
ResidualV_post_TOT_shuff{rix} = ResidualV_TOT_shuff{rix}(:,:, (size(ResidualV_pre_TOT{rix}, 3)+1):end);

end


[SI, poolvar] = selectivity_analytic_prepostshuff(xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT_shuff, ResidualV_post_TOT_shuff, ResidualA_pre_TOT_shuff, ResidualA_post_TOT_shuff,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);

SI_model_pre{S} = SI.SIinputs(:,1) + SI.SIinitstate(:,1) + SI.SIvel(:,1) + SI.SIres(:,1);
SI_model_post{S} = SI.SIinputs(:,2) + SI.SIinitstate(:,2) + SI.SIvel(:,2) + SI.SIres(:,2);
% SI_model_pre_deltares = SI_model_pre_full - SI.SIdeltares(:,1)  ;  % amount generated without single-cell residuals
% SI_model_post_deltares = SI_model_post_full - SI.SIdeltares(:,2) ; % amount generated without single-cell residuals

for j=1:4
    Model_pre{S}(j) = mean(abs(SI_model_pre{S}(CELLLAB_TOT == CLS(j))));
    Model_post{S}(j) = mean(abs(SI_model_post{S}(CELLLAB_TOT == CLS(j))));
    Model_pre_full(j) = mean(abs(SI_model_pre_full(CELLLAB_TOT == CLS(j))));
    Model_post_full(j) = mean(abs(SI_model_post_full(CELLLAB_TOT == CLS(j))));
%     Model_pre_deltares(j) = mean(abs(SI_model_pre_deltares(CELLLAB_TOT == CLS(j))));
%     Model_post_deltares(j) = mean(abs(SI_model_post_deltares(CELLLAB_TOT == CLS(j))));
end
 
end
 
%Delete weights from W2 to W1
% 
% var_method = 'fixedvar_rw';
% manipulation = 'removeall';
% 
% clear SI
% clear poolvar
% 
% tsims = 9:17;
% for W1 = [1,3,4,5]
%     for W2 = [1,3,4,5]
% 
%         [SI{W1,W2}, poolvar{W1,W2}] = selectivity_analytic(manipulation, var_method, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);
%        
%     end
% end

