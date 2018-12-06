
clear all
changeW = ones(4);
changeI = [1,1,1,1];

fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;

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

end


for Cond = 1:Nconds

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

drmatTOT_post_train{2} = drmatTOT_post{2};
rmatTOT_post_train{2}  = rmatTOT_post{2};

drmatTOT_post_train{3} = drmatTOT_post{3};
rmatTOT_post_train{3}  = rmatTOT_post{3};

for Cond = 2:Nconds

    drmatSTM_post = drmatTOT_post_train{Cond};
    rmatSTM_post = rmatTOT_post_train{Cond};

    drmatSTM_pre = drmatTOT_pre{Cond};
    rmatSTM_pre = rmatTOT_pre{Cond};
    
    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]

    M = permute(drmatSTM_pre, [3,1,2]);
    drmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_pre, [3,1,2]);
    rmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series


     rmat0post{Cond}= rmatSTM_post0.';
     drmat0post{Cond} = drmatSTM_post0.';
 
     rmat0pre{Cond} = rmatSTM_pre0.';
     drmat0pre{Cond} = drmatSTM_pre0.';


end

Cond = 1;

    drmatSTM_pre = drmatTOT_pre{Cond};
    rmatSTM_pre = rmatTOT_pre{Cond};
    
    M = permute(drmatSTM_pre, [3,1,2]);
    drmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_pre, [3,1,2]);
    rmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series

    rmat0pre{Cond} = rmatSTM_pre0.';
    drmat0pre{Cond} = drmatSTM_pre0.';
    
drmatTOT_post_s1 = drmatTOT_post{1};
rmatTOT_post_s1  = rmatTOT_post{1};

Error_LDS_temp = zeros([Ntrialspost{1}, length(CELLLAB_ALL{rix}), 16]);
Error_meanrate_temp = zeros([Ntrialspost{1}, length(CELLLAB_ALL{rix}), 16]);

parfor i=1:Ntrialspost{1}


% clear rmat0rel
% clear drmat0rel
% clear rmat0pre
% clear drmat0pre

TrainInds = [1:(i-1), (i+1):Ntrialspost{1}];

drmatTOT_post_train_s1 = drmatTOT_post_s1(TrainInds,:,:);
rmatTOT_post_train_s1  = rmatTOT_post_s1(TrainInds,:,:);

drmatTOT_post_test_s1 = drmatTOT_post_s1(i,:,:);
rmatTOT_post_test_s1  = rmatTOT_post_s1(i,:,:);


    drmatSTM_post = drmatTOT_post_train_s1;
    rmatSTM_post = rmatTOT_post_train_s1;


    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]


    CLS = [1,3,4,5];

    rmat0post_dum = {rmatSTM_post0.', rmat0post{2}, rmat0post{3}};
    drmat0post_dum = {drmatSTM_post0.', drmat0post{2}, drmat0post{3}};
%     rmat0rel{1} = rmatSTM_rel0.';
%     drmat0rel{1} = drmatSTM_rel0.';
%  
%     rmat0pre{1} = rmatSTM_pre0.';
%     drmat0pre{1} = drmatSTM_pre0.';



Ntrialspost_train = Ntrialspost;
Ntrialspost_train{1} = Ntrialspost{1}-1;


%% Set up regression matrix equations

[xpost, xpre, Measurement_post, Measurement_pre, Prediction_post, Prediction_pre] = fitLDS_nospeed(fixvars, changeW, changeI, drmat0post_dum, drmat0pre, rmat0post_dum, rmat0pre, rngtot, Ntrialspost_train, Ntrialspre, CELLLAB_ALL{rix});

residualpost_train  = Measurement_post - Prediction_post;
residualpre_train = Measurement_pre - Prediction_pre;

%% Get validation errors under each model


r_trials = cell([3,1]);
r_trials_TOT = cell([3,1]);
MeanTOT_train_r = cell([3,1]);

    
    MeanTOT_train_r{1} = mean(rmatTOT_post_train_s1);

    for T=1:size(drmatTOT_post_train_s1,1)

        dum = rmatTOT_post_train_s1(T,:,:) - MeanTOT_train_r{1};
        r_trials_TOT{1}(T,:,:) = dum;
        r_trials{1}(:,T) = dum(:); 
        
    end



TRL_data = squeeze(rmatTOT_post_test_s1(1,:,:));
TRL_mean = squeeze(MeanTOT_train_r{1});

Amat = xpost(:, 1:size(xpost,1)) + eye(size(xpost,1));
Inputs = xpost(:, [(size(xpost,1)+1):(size(xpost,1)+17)]);


Data = TRL_data(:,2:17) - Amat * TRL_data(:,1:16) - Inputs(:,1:16);
Error_LDS_temp(i,:,:) = Data;

Data = TRL_data(:,2:17) - repmat(mean(TRL_mean,2), [1,16]);
Error_meanrate_temp(i,:,:) = Data;

end

Error_LDS{rix} = Error_LDS_temp;
Error_MeanRate{rix} = Error_meanrate_temp;

    end

end

    
   
    



%% Get final outputs

CELLLAB_TOT = vertcat(CELLLAB_ALL{:});

RXS = [1,3,4,5,6,8,9,10];

for rix = RXS

     
    LDS_squarederror_percell{rix} = Error_LDS{rix}.^2 ./ length(CELLLAB_ALL{rix});
 
    
    % R-squared 
  
    
    LDS_Rsq{rix} = squeeze(1 - sum(sum(Error_LDS{rix}.^2,3),1) ./ sum(sum(Error_MeanRate{rix}.^2,3),1));
  
end

MEAN = mean(vertcat(LDS_Rsq{:}));
SEM  = std(vertcat(LDS_Rsq{:}));

MEAN = MEAN';
SEM = SEM';
% 
% subplot(2,1,1)
% 
% hold on
% 
% hb = bar(1:size(MEAN,1), MEAN);
% 
% for ib = 1:numel(hb)
% 
%       % Find the centers of the bars
%       xData = get(get(hb(ib),'Children'),'XData');
%       barCenters = mean(unique(xData,'rows'));
% 
%       errorbar(barCenters,MEAN(:,ib),SEM(:,ib),'k.', 'linewidth', 3);
% 
% end
% 
% set(gca, 'fontsize', 18)
% ylabel('R^2 (mean over cells)')
% set(gca, 'xtick', 1:6)
% set(gca, 'xticklabel', {'PSTH', 'LDS', 'No Interactions', 'Time-varying Interactions', 'Self-prediction (all lags)', 'Self-prediction (one lag)'})
% 
% 
% subplot(2,2,3)
% hold on
% 
% hist(vertcat(LDS_Rsq{:}) - vertcat(Lagone_autocorr_singlemodel_Rsq{:}), 100)
% 
% 
% set(gca, 'fontsize', 18)
% xlabel('\Delta R^2 (LDS vs no interactions)')
% ylabel('Number of cells')
% 
% subplot(2,2,4)
% 
% hist(vertcat(LDS_Rsq{:}) - vertcat(Lagone_Rsq{:}), 100)
% 
% set(gca, 'fontsize', 18)
% xlabel('\Delta R^2 (LDS vs time varying interactions)')
% ylabel('Number of cells')
% 
