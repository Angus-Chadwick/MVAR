%% Equations to do regression on linear dynamical system from individual trials

clear all;

changeI = [0,0,0,0];
changeW = zeros(4);


TestRelIrrel = 'Rel'; % choose whether to predict relevant or irrelevant trials

fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;


%clear all;

% attempt to match post learning to switching trials NB: missed out one
% repeated switching session for M75

 INCLUDE_SD = { % rsbs_sd_pop_prepro_cutoffs
    'M70_20141028_B1' 'M70_20141106_B1'  
    'M73_20141028_B1' 'M73_20141101_B1'
    'M75_20141029_B1' 'M75_20141102_B1'
    'M80_20141028_B1' 'M80_20141031_B1'
    % 'M81_20141031_B1'
    'M87_20141028_B1' 'M87_20141108_B1'
    'M89_20141029_B1' 'M89_20141030_B1'
    'M93_20141028_B1' 'M93_20141111_B1'
     } % M93_20141023_B1






for rix=1:size(INCLUDE_SD,1),
       
    clear TOT;
    for PSTVIS=1:2,
        
        if PSTVIS == 1
            
            ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/Learning/';
            RSG.savedirforAngus = [ROOT,'Saved_Data/'];
            
            if ~exist('TSK'), TSK=[]; end
            if ~isfield(TSK,'LN'),
                TSK.LN=load([RSG.savedirforAngus,'LN_TLTPDPOP.mat']);
            end
            
            clear idinf;NID = 1;
            name = INCLUDE_SD{rix,PSTVIS}
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
       
                fname=sprintf('%sADDTRG2_alltrials_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
                load(fname,'ADDTRG');
            
        
            % save data needed for later
            TOT{PSTVIS}.RIXSES = RIXSES;
            TOT{PSTVIS}.ADDTRG = ADDTRG;    
            
        elseif PSTVIS == 2
            
            ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/';
            RSG.savedirforAngus = [ROOT,'Saved_Data/'];
        
            if ~exist('TSK'), TSK=[]; end
            if ~isfield(TSK,'SW'),
                TSK.SW=load([RSG.savedirforAngus,'SW_TLTPDPOP.mat']);
            end
            
            clear idinf;NID = 1;
            name = INCLUDE_SD{rix,PSTVIS}
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
            fname=sprintf('%sADDTRG2_alltrials_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            load(fname,'ADDTRG');      
       
            % save data needed for later
            TOT{PSTVIS}.RIXSES = RIXSES;
            TOT{PSTVIS}.ADDTRG = ADDTRG;    
            
        end
        
    end
        
        
     
            % find same cells pre and post
        CELL{1} = TOT{1}.ADDTRG{1}.POOLCV; % PRE
        CELL{2} = TOT{2}.ADDTRG{1}.POOLCV; % PST
        CELLUNI = unique([CELL{1};CELL{2}],'rows');
        IX=NaN(size(CELLUNI,1),2);
        for i=1:size(CELLUNI,1),
            for PSTVIS=1:2,
                ix=find(CELL{PSTVIS}(:,1)==CELLUNI(i,1)&CELL{PSTVIS}(:,2)==CELLUNI(i,2));
                if ~isempty(ix),
                    IX(i,PSTVIS)=ix;
                end
            end
        end
        IX = IX(all(~isnan(IX),2),:); % all cells with pre and post data
        
           

%% Get cell labels

 POOLCVPRE = CELL{1};
 POOLCVPOST = CELL{2};
 
        % find cell type of each cell
        nam = 'LN';
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
                nam = 'SW';
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
        
        
   
 %% 
 
 

CondsPST = [1,2,7];  % Vertical, angled and grey onset
CondsVIS = [1,2,12];

Nconds = length(CondsPST);


for Cond=1:Nconds

    CondiPST = CondsPST(Cond);
    CondiVIS = CondsVIS(Cond);

    rmat_post{rix, Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{CondiPST}.group1(:, IX(:,1), :);
    rmat_rel{rix, Cond}  = TOT{2}.ADDTRG{1}.PL{1}.DCOR{CondiVIS}.group1(:, IX(:,2), :);


    drmat_rel{rix, Cond} = diff(rmat_rel{rix, Cond},1,3);
    drmat_post{rix, Cond} = diff(rmat_post{rix, Cond},1,3);



    tsamples = TOT{1}.ADDTRG{1}.PL{1}.DCOR{CondiPST}.t;
    rngSTM = TOT{1}.ADDTRG{1}.PL{1}.DCOR{CondiPST}.rngSTM;
    rngBSL = TOT{1}.ADDTRG{1}.PL{1}.DCOR{CondiPST}.rngBSL;
    rngtot = [rngBSL, rngSTM];
    Ntrialsrel{Cond} = size(drmat_rel{rix, Cond},1);
    Ntrialspost{Cond} = size(drmat_post{rix, Cond},1);
    Nsamples = length(rngtot);
    tMOT{Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{CondiPST}.tMOT;

    drmatTOT_rel{Cond} = drmat_rel{rix, Cond}(:,:,rngtot);
    rmatTOT_rel{Cond} = rmat_rel{rix, Cond}(:,:,rngtot);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);

    
   

    
clear Nrel
clear Npost

for i=1:Ntrialsrel{Cond}

    M = squeeze(drmatTOT_rel{Cond}(i,:,:));
    Nrel(i) = sum(sum(isnan(M)));

end

for i=1:Ntrialspost{Cond}

    M = squeeze(drmatTOT_post{Cond}(i,:,:));
    Npost(i) = sum(sum(isnan(M)));

end

Trialsrel{Cond} = find(Nrel == 0);
Trialspost{Cond} = find(Npost == 0);

drmatTOT_rel{Cond} = drmatTOT_rel{Cond}(Trialsrel{Cond},:,:);
rmatTOT_rel{Cond} = rmatTOT_rel{Cond}(Trialsrel{Cond},:,:);

drmatTOT_post{Cond} = drmatTOT_post{Cond}(Trialspost{Cond},:,:);
rmatTOT_post{Cond} = rmatTOT_post{Cond}(Trialspost{Cond},:,:);

Ntrialsrel{Cond} = length(Trialsrel{Cond});
Ntrialspost{Cond} = length(Trialspost{Cond});

end




%% Fit full (V & A, BSL & STM) simultaneously

drmatTOT_rel_train{2} = drmatTOT_rel{2};
rmatTOT_rel_train{2}  = rmatTOT_rel{2};

drmatTOT_rel_train{3} = drmatTOT_rel{3};
rmatTOT_rel_train{3}  = rmatTOT_rel{3};

for Cond = 2:Nconds

    drmatSTM_rel = drmatTOT_rel_train{Cond};
    rmatSTM_rel = rmatTOT_rel_train{Cond};

    drmatSTM_post = drmatTOT_post{Cond};
    rmatSTM_post = rmatTOT_post{Cond};
    
    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]

    M = permute(drmatSTM_rel, [3,1,2]);
    drmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_rel, [3,1,2]);
    rmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series


     rmat0post{Cond}= rmatSTM_post0.';
     drmat0post{Cond} = drmatSTM_post0.';
 
     rmat0rel{Cond} = rmatSTM_rel0.';
     drmat0rel{Cond} = drmatSTM_rel0.';


end

Cond = 1;

    drmatSTM_post = drmatTOT_post{Cond};
    rmatSTM_post = rmatTOT_post{Cond};
    
    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series

    rmat0post{Cond} = rmatSTM_post0.';
    drmat0post{Cond} = drmatSTM_post0.';
    
drmatTOT_rel_s1 = drmatTOT_rel{1};
rmatTOT_rel_s1  = rmatTOT_rel{1};


%%
                
    
Error_LDS_temp = zeros([Ntrialsrel{1}, length(CELLLAB_ALL{rix}), 16]);
Error_meanrate_temp = zeros([Ntrialsrel{1}, length(CELLLAB_ALL{rix}), 16]);

parfor i=1:Ntrialsrel{1}


% clear rmat0rel
% clear drmat0rel
% clear rmat0irrel
% clear drmat0irrel

TrainInds = [1:(i-1), (i+1):Ntrialsrel{1}];

drmatTOT_rel_train_s1 = drmatTOT_rel_s1(TrainInds,:,:);
rmatTOT_rel_train_s1  = rmatTOT_rel_s1(TrainInds,:,:);

drmatTOT_rel_test_s1 = drmatTOT_rel_s1(i,:,:);
rmatTOT_rel_test_s1  = rmatTOT_rel_s1(i,:,:);


    drmatSTM_rel = drmatTOT_rel_train_s1;
    rmatSTM_rel = rmatTOT_rel_train_s1;


    M = permute(drmatSTM_rel, [3,1,2]);
    drmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_rel, [3,1,2]);
    rmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]


    CLS = [1,3,4,5];

    rmat0rel_dum = {rmatSTM_rel0.', rmat0rel{2}, rmat0rel{3}};
    drmat0rel_dum = {drmatSTM_rel0.', drmat0rel{2}, drmat0rel{3}};
%     rmat0rel{1} = rmatSTM_rel0.';
%     drmat0rel{1} = drmatSTM_rel0.';
%  
%     rmat0irrel{1} = rmatSTM_irrel0.';
%     drmat0irrel{1} = drmatSTM_irrel0.';



Ntrialsrel_train = Ntrialsrel;
Ntrialsrel_train{1} = Ntrialsrel{1}-1;


%% Set up regression matrix equations

[xrel, xpost, Measurement_rel, Measurement_post, Prediction_rel, Prediction_post] = fitLDS_nospeed(fixvars, changeW, changeI, drmat0rel_dum, drmat0post, rmat0rel_dum, rmat0post, rngtot, Ntrialsrel_train, Ntrialspost, CELLLAB_ALL{rix});

residualrel_train  = Measurement_rel - Prediction_rel;
residualpost_train = Measurement_post - Prediction_post;

%% Get validation errors under each model


r_trials = cell([3,1]);
r_trials_TOT = cell([3,1]);
MeanTOT_train_r = cell([3,1]);

    
    MeanTOT_train_r{1} = mean(rmatTOT_rel_train_s1);

    for T=1:size(drmatTOT_rel_train_s1,1)

        dum = rmatTOT_rel_train_s1(T,:,:) - MeanTOT_train_r{1};
        r_trials_TOT{1}(T,:,:) = dum;
        r_trials{1}(:,T) = dum(:); 
        
    end



TRL_data = squeeze(rmatTOT_rel_test_s1(1,:,:));
TRL_mean = squeeze(MeanTOT_train_r{1});

Amat = xrel(:, 1:size(xrel,1)) + eye(size(xrel,1));
Inputs = xrel(:, [(size(xrel,1)+1):(size(xrel,1)+17)]);


Data = TRL_data(:,2:17) - Amat * TRL_data(:,1:16) - Inputs(:,1:16);
Error_LDS_temp(i,:,:) = Data.^2;
   
Data = TRL_data(:,2:17) - repmat(mean(TRL_mean,2), [1,16]);
Error_meanrate_temp(i,:,:) = Data.^2;

end


Error_LDS{rix} = Error_LDS_temp;
Error_MeanRate{rix} = Error_meanrate_temp;
    

end

%% Get final outputs

CELLLAB_TOT = vertcat(CELLLAB_ALL{:});

RXS = 1:29;

for rix = RXS

     
    LDS_squarederror_percell{rix} = Error_LDS{rix} ./ length(CELLLAB_ALL{rix});
 
    
    % R-squared 
  
    
    LDS_Rsq{rix} = squeeze(1 - sum(sum(Error_LDS{rix},3),1) ./ sum(sum(Error_MeanRate{rix},3),1));
  
end

MEAN = mean(vertcat(LDS_Rsq{:}));
SEM  = std(vertcat(LDS_Rsq{:}));

MEAN = MEAN';
SEM = SEM';
