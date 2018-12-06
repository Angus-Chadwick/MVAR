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

    if strmatch(TestRelIrrel, 'Rel')
        
       Conds_rel = [1,2,12];  % Vertical, angled, grey        
       Conds_irrel = [5,6,13]; % Vertical, angled, grey
       
    elseif strmatch(TestRelIrrel, 'Irrel')
        
       Conds_irrel = [1,2,12];  % Vertical, angled, grey        
       Conds_rel = [5,6,13]; % Vertical, angled, grey
       
    end
       
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

drmatTOT_rel_train{2} = drmatTOT_rel{2};
rmatTOT_rel_train{2}  = rmatTOT_rel{2};

drmatTOT_rel_train{3} = drmatTOT_rel{3};
rmatTOT_rel_train{3}  = rmatTOT_rel{3};

for Cond = 2:Nconds

    drmatSTM_rel = drmatTOT_rel_train{Cond};
    rmatSTM_rel = rmatTOT_rel_train{Cond};

    drmatSTM_irrel = drmatTOT_irrel{Cond};
    rmatSTM_irrel = rmatTOT_irrel{Cond};
    
    M = permute(drmatSTM_rel, [3,1,2]);
    drmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_rel, [3,1,2]);
    rmatSTM_rel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]

    M = permute(drmatSTM_irrel, [3,1,2]);
    drmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_irrel, [3,1,2]);
    rmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series


     rmat0rel{Cond}= rmatSTM_rel0.';
     drmat0rel{Cond} = drmatSTM_rel0.';
 
     rmat0irrel{Cond} = rmatSTM_irrel0.';
     drmat0irrel{Cond} = drmatSTM_irrel0.';


end

Cond = 1;

    drmatSTM_irrel = drmatTOT_irrel{Cond};
    rmatSTM_irrel = rmatTOT_irrel{Cond};
    
    M = permute(drmatSTM_irrel, [3,1,2]);
    drmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_irrel, [3,1,2]);
    rmatSTM_irrel0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series

    rmat0irrel{Cond} = rmatSTM_irrel0.';
    drmat0irrel{Cond} = drmatSTM_irrel0.';
    
drmatTOT_rel_s1 = drmatTOT_rel{1};
rmatTOT_rel_s1  = rmatTOT_rel{1};

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

[xrel, xirrel, Measurement_rel, Measurement_irrel, Prediction_rel, Prediction_irrel] = fitLDS_Switching_nospeed(fixvars, changeW, changeI, drmat0rel_dum, drmat0irrel, rmat0rel_dum, rmat0irrel, rngtot, Ntrialsrel_train, Ntrialsirrel, CELLLAB_ALL{rix});

residualrel_train  = Measurement_rel - Prediction_rel;
residualirrel_train = Measurement_irrel - Prediction_irrel;

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