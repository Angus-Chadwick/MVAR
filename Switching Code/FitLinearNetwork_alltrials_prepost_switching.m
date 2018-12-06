%% Equations to do regression on linear dynamical system from individual trials

clear all;

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

changeW = ones(4);
changeI = [1,1,1,1];

[xrel{rix}, xirrel{rix}, Measurement_rel{rix}, Measurement_irrel{rix}, Prediction_rel{rix}, Prediction_irrel{rix}] = fitLDS_Switching(fixvars, changeW, changeI, drmat0rel, drmat0irrel, rmat0rel, rmat0irrel, rngtot, Ntrialsrel, Ntrialsirrel, vmatSTM_rel0, vmatSTM_irrel0, CELLLAB_ALL{rix});
%[xrel{rix}, xirrel{rix}, Measurement_rel{rix}, Measurement_irrel{rix}, Prediction_rel{rix}, Prediction_irrel{rix}] = fitLDS(fixvars, changeW, changeI, drmat0rel, drmat0irrel, rmat0rel, rmat0irrel, rngtot, Ntrialsrel, Ntrialsirrel, vmatSTM_rel0, vmatSTM_irrel0, CELLLAB_ALL{rix});

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

%% Generate predicted vs measured traces for vertical stimulus 
%     
% 
% for rix = 1:length(RXS)
%     for T=1:length(horzcat(Trialsrel_all{RXS(rix)}{1}))
%         
%         Predicted_rel{rix}(:,:,T) = Prediction_rel{RXS(rix)}(:, (17 * (T-1) + 1):(17*T)) + permute(rmat_rel{RXS(rix),1}(Trialsrel_all{RXS(rix)}{1}(T),:, rngtot), [2,3,1]);
%         Measured_rel{rix}(:,:,T) = permute(rmat_rel{RXS(rix),1}(Trialsrel_all{RXS(rix)}{1}(T),:, rngtot+1), [2,3,1]);
%         
% 
%     end
% end
% 
% for rix = 1:length(RXS)
%     for T=1:length(horzcat(Trialsrel_all{RXS(rix)}{1}))
%         
%         Predicted_rel{rix}(:,:,T) = Prediction_rel{RXS(rix)}(:, (17 * (T-1) + 1):(17*T)) + permute(rmat_rel{RXS(rix),1}(Trialsrel_all{RXS(rix)}{1}(T),:, rngtot), [2,3,1]);
%         Measured_rel{rix}(:,:,T) = permute(rmat_rel{RXS(rix),1}(Trialsrel_all{RXS(rix)}{1}(T),:, rngtot+1), [2,3,1]);
%         
%     end
% end
% 
% figure
% hold on
% 
% rix = 4;  % session number
% n = 31;  % neuron number
% Trl1 = 4;  % trial number
% Trl2 = 30;
% 
% plot(tsamples(rngtot + 1), mean(Measured_rel{rix}(n,:,:), 3), 'color', 'b', 'linewidth', 3)
% plot(tsamples(rngtot + 1), Measured_rel{rix}(n,:,Trl1), 'color', 'b', 'linestyle', '--', 'linewidth', 2)
% plot(tsamples(rngtot + 1), Predicted_rel{rix}(n,:,Trl1), 'color', 'r', 'linestyle', '--', 'linewidth', 2)
% 
% plot(tsamples(rngtot + 1), Measured_rel{rix}(n,:,Trl2), 'color', 'b', 'linestyle', '-.', 'linewidth', 2)
% plot(tsamples(rngtot + 1), Predicted_rel{rix}(n,:,Trl2), 'color', 'r', 'linestyle', '-.', 'linewidth', 2)
% % 
% box on
% set(gca, 'fontsize', 22)
% xlabel('Time (s)')
% ylabel('dF/F')
% legend('PSTH', 'Example trial 1 measured', 'Example trial 1 predicted', 'Example trial 2 measured', 'Example trial 2 predicted')
% 
% for rix = 1:length(RXS)
%     for i=1:length(CELLLAB_ALL{RXS(rix)})
%         
%         Rsquared_pre{rix}(i) = corr(drmatpre_all{RXS(rix)}(i,:).' - residualpre{RXS(rix)}(i,:).', drmatpre_all{RXS(rix)}(i,:).').^2;
%         Rsquared_post{rix}(i) = corr(drmatpost_all{RXS(rix)}(i,:).' - residualpost{RXS(rix)}(i,:).', drmatpost_all{RXS(rix)}(i,:).').^2;
%         
%         Rsquared_pre0{rix}(i) = 1 - sum(residualpre{RXS(rix)}(i,:).^2) / sum((drmatpre_all{RXS(rix)}(i,:) - mean(drmatpre_all{RXS(rix)}(i,:))).^2);
%         Rsquared_post0{rix}(i) = 1 - sum(residualpost{RXS(rix)}(i,:).^2) / sum((drmatpost_all{RXS(rix)}(i,:) - mean(drmatpost_all{RXS(rix)}(i,:))).^2);
%                 
%         Rsquared_pre0_r{rix}(i) = 1  - sum(residualpre{RXS(rix)}(i,:).^2) / sum((rmatpre_all{RXS(rix)}(i,:) - mean(rmatpre_all{RXS(rix)}(i,:))).^2);
%         Rsquared_post0_r{rix}(i) = 1 - sum(residualpost{RXS(rix)}(i,:).^2) / sum((rmatpost_all{RXS(rix)}(i,:) - mean(rmatpost_all{RXS(rix)}(i,:))).^2);
%         
%     end
%     
%      Rsquared_pre00(rix)  = 1 - sum(sum(residualpre{RXS(rix)}.^2)) / sum(sum((drmatpre_all{RXS(rix)} - repmat(mean(drmatpre_all{RXS(rix)},2), [1, size(drmatpre_all{RXS(rix)},2)])).^2));
%      Rsquared_post00(rix) = 1 - sum(sum(residualpost{RXS(rix)}.^2)) / sum(sum((drmatpost_all{RXS(rix)} - repmat(mean(drmatpost_all{RXS(rix)},2), [1, size(drmatpost_all{RXS(rix)},2)])).^2));
% 
%      Rsquared_pre00_r(rix)  = 1 - sum(sum(residualpre{RXS(rix)}.^2)) / sum(sum((rmatpre_all{RXS(rix)} - repmat(mean(rmatpre_all{RXS(rix)},2), [1, size(drmatpre_all{RXS(rix)},2)])).^2));
%      Rsquared_post00_r(rix) = 1 - sum(sum(residualpost{RXS(rix)}.^2)) / sum(sum((rmatpost_all{RXS(rix)} - repmat(mean(rmatpost_all{RXS(rix)},2), [1, size(drmatpost_all{RXS(rix)},2)])).^2));
% 
%     
% end
% 
% Rsquared_pre_tot = horzcat(Rsquared_pre{:});
% Rsquared_post_tot = horzcat(Rsquared_post{:});


%% Relationship between weights and selectivity, with possible explanation of earlier cell-type specific prediction analyses

% Plot inputs (weighted average over sessions)

clear Inputsrel_norm
clear Inputsirrel_norm
clear Inputsrel_v_norm
clear Inputsirrel_v_norm
clear Meansrel
clear Meansirrel
clear Meansvrel
clear Meansvirrel
clear STDrel
clear STDirrel
clear N


selective = 0;

CLSNAM = {'PYR', 'PV', 'SOM', 'VIP'};

for i=1:length(CLS)
    
    Inputsrel{i} = [];
    Inputsirrel{i} = [];
    Inputsrel_vel_coeff{i} = [];
    Inputsirrel_vel_coeff{i} = [];
    InputsrelV_vel{i} = [];
    InputsirrelV_vel{i} = [];
    InputsrelA_vel{i} = [];
    InputsirrelA_vel{i} = [];
             
     for rix=1:length(RXS)
         
         Inputsrel{i} = vertcat(Inputsrel{i}, xrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), (length(CELLLAB_ALL{RXS(rix)}) + 1):(end-1)));
         Inputsirrel{i} = vertcat(Inputsirrel{i}, xirrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), (length(CELLLAB_ALL{RXS(rix)}) + 1):(end-1)));

         Inputsrel_vel_coeff{i}  = vertcat(Inputsrel_vel_coeff{i},  xrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end));
         Inputsirrel_vel_coeff{i} = vertcat(Inputsirrel_vel_coeff{i}, xirrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end));
         
         VelmeanV_rel(rix) = mean(vmatrel_all{RXS(rix), 1}(9:end));
         VelmeanV_irrel(rix) = mean(vmatirrel_all{RXS(rix), 1}(9:end));
        
         VelmeanA_rel(rix) = mean(vmatrel_all{RXS(rix), 2}(9:end));
         VelmeanA_irrel(rix) = mean(vmatirrel_all{RXS(rix), 2}(9:end));

         InputsrelV_vel{i}  = vertcat(InputsrelV_vel{i},  xrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end)  * VelmeanV_rel(rix));
         InputsirrelV_vel{i} = vertcat(InputsirrelV_vel{i}, xirrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end) * VelmeanV_irrel(rix));
         
         InputsrelA_vel{i}  = vertcat(InputsrelA_vel{i},  xrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end)  * VelmeanA_rel(rix));
         InputsirrelA_vel{i} = vertcat(InputsirrelA_vel{i}, xirrel{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end) * VelmeanA_irrel(rix));
         
     end
                  
     if selective 
         
        Meansrel(i,:)     = mean(Inputsrel{i}(psel_rel_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        Meansirrel(i,:)    = mean(Inputsirrel{i}(psel_irrel_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        
        STDrel = std(Inputsrel{i}(psel_rel_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        STDirrel = std(Inputsirrel{i}(psel_irrel_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        SEMrel = STDrel/sqrt(sum(psel_rel_tot(CELLLAB_TOT == CLS(i)) < 0.05));
        SEMirrel = STDirrel/sqrt(sum(psel_irrel_tot(CELLLAB_TOT == CLS(i)) < 0.05));
        
     else
     
        Meansrel(i,:)     = mean(Inputsrel{i});
        Meansirrel(i,:)    = mean(Inputsirrel{i});
        
        STDrel = std(Inputsrel{i});
        STDirrel = std(Inputsirrel{i});
        SEMrel = STDrel/sqrt(size(Inputsrel{i},1));
        SEMirrel = STDirrel/sqrt(size(Inputsirrel{i},1));
    
     end
       
    Meansrel_vel(i)     = mean(Inputsrel_vel_coeff{i});
    Meansirrel_vel(i)    = mean(Inputsirrel_vel_coeff{i});
    
    Nvel(i) = length(Inputsrel_vel_coeff{i});
    
    STDrel_vel(i) = std(Inputsrel_vel_coeff{i});
    STDirrel_vel(i) = std(Inputsirrel_vel_coeff{i});
  
    subplot(2,2,i)
    hold on
      errorbar(tsamples(rngtot+1), (Meansrel(i,1:17)),  SEMrel(1:17), 'color', 'b', 'linewidth', 4, 'linestyle', '--')
      errorbar(tsamples(rngtot+1), (Meansirrel(i,1:17)), SEMirrel(1:17),'color', 'b', 'linewidth', 4, 'linestyle', '-')
      errorbar(tsamples(rngtot+1), (Meansrel(i,18:34)), SEMrel(18:34), 'color', 'r', 'linewidth', 4, 'linestyle', '--')
      errorbar(tsamples(rngtot+1), (Meansirrel(i,18:34)), SEMirrel(18:34),'color', 'r', 'linewidth', 4, 'linestyle', '-')
      
      plot( (tsamples(rngtot(8)) + tsamples(rngtot(9))) / 2 * [1,1], [0, 0.03], 'color', 'k', 'linewidth', 3, 'linestyle', '--')
            
   box on
   set(gca, 'fontsize', 22)
   xlabel('Time (s)')
   ylabel('External Input')
   legend('Vertical rel', 'Vertical irrel', 'Angled rel', 'Angled irrel')
   title(CLSNAM(i))
   axis([-0.5, 1.1, -0.1, 0.1])
   
end
   
   figure
   hold on
   errorbar([1:4], Meansrel_vel, STDrel_vel ./ sqrt(Nvel), '.k', 'linewidth', 3)
   errorbar([1:4] + 0.1, Meansirrel_vel, STDirrel_vel ./ sqrt(Nvel), '.r', 'linewidth', 3)
   set(gca, 'xtick', 1:4)
   set(gca, 'xticklabel', {'PYR', 'PV', 'SOM', 'VIP'})
   box on
   set(gca, 'fontsize', 18)
   ylabel('Velocity Input Gain')
   legend('rel', 'irrel')
   
   
%% Animal by animal output selectivity

% NOTE: We have changed this from learning so that only trials with no
% nans are included now. This convention was used in most of the model-based learning
% analyses, but it was never incorporated at this stage of raw data.

RestrictTrials = 1;  % restrict to trials without nans to match model-based analyses

for rix=1:length(RXS)

    if RestrictTrials
    
        dMrel{rix} = squeeze(nanmean(nanmean(rmat_rel{RXS(rix), 1}(Trialsrel_all{rix}{1},:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_rel{RXS(rix), 2}(Trialsrel_all{rix}{2},:,rngSTM), 3)));
        dMirrel{rix} = squeeze(nanmean(nanmean(rmat_irrel{RXS(rix), 1}(Trialsirrel_all{rix}{1},:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_irrel{RXS(rix), 2}(Trialsirrel_all{rix}{2},:,rngSTM), 3)));

    else
          
        dMrel{rix} = squeeze(nanmean(nanmean(rmat_rel{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_rel{RXS(rix), 2}(:,:,rngSTM), 3)));
        dMirrel{rix} = squeeze(nanmean(nanmean(rmat_irrel{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_irrel{RXS(rix), 2}(:,:,rngSTM), 3)));
   
    end
    
    SI_Out_rel{rix} = dMrel{rix} ./ PoolVar_rel_Out{RXS(rix)};
%     
%     for c=1:size(rmat_rel{RXS(rix)}, 2)
%     
%         SI_Out_rel_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_rel{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_rel{RXS(rix), 2}(:,c,rngSTM), 3)));
%         SI_Out_irrel_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_irrel{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_irrel{RXS(rix), 2}(:,c,rngSTM), 3)));
% 
%     end
    
    SI_Out_irrel{rix} = dMirrel{rix} ./ PoolVar_irrel_Out{RXS(rix)};
    
end


%% Animal by animal input selectivity 

for rix=1:length(RXS)

          MeanVrel = mean(xrel{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanArel = mean(xrel{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVirrel = mean(xirrel{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanAirrel = mean(xirrel{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
          
          SI_In_rel{rix} = (MeanVrel - MeanArel) ./ PoolVar_rel_Out{RXS(rix)}';
          SI_In_irrel{rix} = (MeanVirrel - MeanAirrel) ./ PoolVar_irrel_Out{RXS(rix)}';

end   

%% Calculate and Plot Mean Connectivity

CONrel = cell(4);
CONirrel = CONrel;

for rix=1:length(RXS)

    Arel = xrel{RXS(rix)}(:, 1:size(xrel{RXS(rix)},1));
    Airrel = xirrel{RXS(rix)}(:, 1:size(xirrel{RXS(rix)},1));

    Arel(1:size(Arel,1)+1:end) = NaN; % remove self-weights
    Airrel(1:size(Arel,1)+1:end) = NaN; % remove self-weights
    
    for i=1:4
        for j=1:4
    
            weights = Arel(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
            CONrel{i,j} = vertcat(CONrel{i,j}, weights(:));
            weights = Airrel(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
            CONirrel{i,j} = vertcat(CONirrel{i,j}, weights(:));
            
           
        end
    end
       
end

     

CONrel_Mean = cellfun(@nanmean, CONrel);
CONirrel_Mean = cellfun(@nanmean, CONirrel);
CONrel_SEM = cellfun(@nanstd, CONrel)./sqrt(cellfun(@length, CONrel));
CONirrel_SEM = cellfun(@nanstd, CONirrel)./sqrt(cellfun(@length, CONirrel));


hold on
errorbar(CONrel_Mean(:), CONrel_SEM(:),'.b', 'linewidth', 4)
errorbar([1:16] + 0.1, CONirrel_Mean(:), CONirrel_SEM(:),'.r', 'linewidth', 4)
box on
set(gca, 'fontsize', 22)
set(gca, 'xtick', [1:16])
set(gca, 'xticklabel', {'PYR->PYR', 'PYR->PV', 'PYR->SOM', 'PYR->VIP', 'PV->PYR', 'PV->PV', 'PV->SOM', 'PV->VIP', 'SOM->PYR', 'SOM->PV', 'SOM->SOM', 'SOM->VIP', 'VIP->PYR', 'VIP->PV', 'VIP->SOM', 'VIP->VIP'})
rotateXLabels( gca, 45)
ylabel('Weight')
legend('rel', 'irrel')
title('Interaction weights before and after learning')
hold on
plot([0, 16.5], [0 0], 'color', 'k', 'linewidth', 3)
axis([0.5, 16.5, -0.005, 0.025])

    
%% Calculate noise correlations at vertical stimulus

% Test for significance of responses

NC_method = 'samples';
NC_method2 = 'jasper';

for rix=1:length(RXS)

  BSL_rel =  nanmean(rmat_rel{RXS(rix), 1}(:,:, 13:16), 3); % -0.5 to 0 s baseline
  RSP_rel =  nanmean(rmat_rel{RXS(rix), 1}(:,:, rngSTM), 3); % 0 to 1 s response

  for i=1:size(BSL_rel,2)
      
        p_resp_rel{rix}(i) = signtest(RSP_rel(:,i) - BSL_rel(:,i));
        
  end
  
  BSL_irrel =  nanmean(rmat_irrel{RXS(rix), 1}(:,:, 13:16), 3); % -0.5 to 0 s baseline
  RSP_irrel =  nanmean(rmat_irrel{RXS(rix), 1}(:,:, rngSTM), 3); % 0 to 1 s response

  for i=1:size(BSL_irrel,2)
      
        p_resp_irrel{rix}(i) = signtest(RSP_irrel(:,i) - BSL_irrel(:,i));
        
  end
  
  
end
  
corrvec_rel = cell(4);
corrvec_irrel = cell(4);
corrvec_rel_sig = cell(4);
corrvec_irrel_sig = cell(4);

for rix=1:length(RXS)
      
   Rrel  = rmat_rel{RXS(rix), 1}(:,:, rngSTM);
   Rirrel = rmat_irrel{RXS(rix), 1}(:,:, rngSTM);

   if strmatch(NC_method,  'trialmean')
   
       Rrel = mean(Rrel, 3);
       Rirrel = mean(Rirrel,3);
       Rrel_mean = Rrel - repmat(nanmean(Rrel), [size(Rrel,1), 1]);
       Rirrel_mean = Rirrel - repmat(nanmean(Rirrel), [size(Rirrel,1), 1]);
   
   elseif strmatch(NC_method, 'samples')
       
       Rrel_meansub  = permute(Rrel  - repmat(nanmean(Rrel), [size(Rrel,1),1,1]), [3,1,2]);
       Rirrel_meansub = permute(Rirrel - repmat(nanmean(Rirrel), [size(Rirrel,1),1,1]), [3,1,2]);
   
       if strmatch(NC_method2, 'angus')  % my more efficient code

           Rrel_mean  = reshape(Rrel_meansub,  [size(Rrel,1) * length(tsims), size(Rrel,2)]);
           Rirrel_mean = reshape(Rirrel_meansub, [size(Rirrel,1) * length(tsims), size(Rirrel,2)]);

           [Irel, Jrel] = find(isnan(Rrel_mean));
           [Iirrel, Jirrel] = find(isnan(Rirrel_mean));

           Nfrac_rel(rix) = length(unique(Irel)) / size(Rrel_mean,1);
           Nfrac_irrel(rix) = length(unique(Iirrel)) / size(Rirrel_mean,1);

           Rrel_mean(Irel,:) = [];
           Rirrel_mean(Iirrel,:) = [];

           [corrmat_rel{rix}, corrmat_rel_sig{rix}] = corr(Rrel_mean);
           [corrmat_irrel{rix}, corrmat_irrel_sig{rix}] = corr(Rirrel_mean);

       elseif strmatch(NC_method2, 'jasper') % a more direct replication of Jasper's code
           
           for i=1:length(CELLLAB_ALL{RXS(rix)})
               A1_rel = squeeze(Rrel_meansub(:,:,i));
               A1_irrel = squeeze(Rirrel_meansub(:,:,i));
               for j=1:length(CELLLAB_ALL{RXS(rix)})
                   A2_rel = squeeze(Rrel_meansub(:,:,j));
                   dum = corrcov(nancov(A1_rel(:), A2_rel(:)));
                   corrmat_rel{rix}(i,j)=dum(1,2);
                   A2_irrel = squeeze(Rirrel_meansub(:,:,j));
                   dum = corrcov(nancov(A1_irrel(:), A2_irrel(:)));
                   corrmat_irrel{rix}(i,j)=dum(1,2);
               end
           end
       end
           
   end
       
  
      
   corrmat_rel{rix}(p_resp_rel{rix} > 0.05, :) = NaN;
   corrmat_rel{rix}(:,p_resp_rel{rix} > 0.05) = NaN;

   corrmat_irrel{rix}(p_resp_irrel{rix} > 0.05,:) = NaN;
   corrmat_irrel{rix}(:, p_resp_irrel{rix} > 0.05) = NaN;

   
   for i=1:length(CLS)
       for j=1:length(CLS)
           
           dumrel = corrmat_rel{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i),  CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumrel(dumrel == 1) = nan;
           dumirrel = corrmat_irrel{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumirrel(dumirrel == 1) = nan;
           corrvec_rel{i,j}  = vertcat(corrvec_rel{i,j},  dumrel(:));
           corrvec_irrel{i,j} = vertcat(corrvec_irrel{i,j}, dumirrel(:));

           dumrel = corrmat_rel_sig{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i),  CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumirrel = corrmat_irrel_sig{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
           corrvec_rel_sig{i,j}  = vertcat(corrvec_rel_sig{i,j},  dumrel(:));
           corrvec_irrel_sig{i,j} = vertcat(corrvec_irrel_sig{i,j}, dumirrel(:));
           
       end
   end
           
end


%% Analytical measure of total noise correlations with weight manipulations

manipulation = 'removeW1W2';

W1 = 1;
W2 = 1;
tsims = 9:17;

[ctot_rel_allweights{W1,W2}, ctot_irrel_allweights{W1,W2}, ctot_rel_manipW1W2{W1,W2}, ctot_irrel_manipW1W2{W1,W2}] = noisecorr_analytic(manipulation, W1, W2, xrel, xirrel, rmat_rel, rmat_irrel, ResidualV_rel_TOT, ResidualV_irrel_TOT, vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot,CELLLAB_ALL);



for i=1:length(ctot_rel_allweights)
    
    dcorr(i,1) = nanmean(ctot_rel_allweights{i} - ctot_rel_manipW1W2{i});
    dcorr(i,2) = nanmean(ctot_irrel_allweights{i} - ctot_irrel_manipW1W2{i});

    dcorr_sem(i,1) = nanstd(ctot_rel_allweights{i} - ctot_rel_manipW1W2{i}) / sqrt(length(ctot_rel_allweights{i}));
    dcorr_sem(i,2) = nanstd(ctot_irrel_allweights{i} - ctot_irrel_manipW1W2{i}) / sqrt(length(ctot_irrel_allweights{i}));

end
    
hold on

hb = bar(1:size(dcorr,1), -dcorr);

for ib = 1:numel(hb)

      % Find the centers of the bars
      xData = get(get(hb(ib),'Children'),'XData');
      barCenters = mean(unique(xData,'rows'));

      errorbar(barCenters,-dcorr(:,ib),dcorr_sem(:,ib),'k.', 'linewidth', 3);

end

set(gca, 'fontsize', 18)
set(gca, 'xtick', [1:10])
set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
rotateXLabels( gca, 45)
box on
legend('rel', 'irrel')
axis([0.5, 10.5, -0.08, 0.04])

%% Analytical decomposition of noise correlations into sources of underlying variability

removeweights = 0;
shuffleresiduals = 0;

[cvel_rel0, cvel_irrel0, cres_rel0, cres_irrel0, cinit_rel0, cinit_irrel0, cvel_res_rel0, cvel_res_irrel0, cvel_init_rel0, cvel_init_irrel0, cres_init_rel0, cres_init_irrel0, ctot_rel, ctot_irrel] = noisecorr_analytic_decomp(shuffleresiduals, removeweights, xrel, xirrel, rmat_rel, rmat_irrel, ResidualV_rel_TOT, ResidualV_irrel_TOT, vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot, CELLLAB_ALL);

Data_Means = [cres_rel0(:), cres_irrel0(:), cinit_rel0(:), cinit_irrel0(:), cvel_rel0(:), cvel_irrel0(:), cvel_res_rel0(:), cvel_res_irrel0(:), cvel_init_rel0(:), cvel_init_irrel0(:), cres_init_rel0(:), cres_init_irrel0(:)].';
hb = bar(1:size(Data_Means,2),Data_Means');

set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
set(gca, 'xtick', [1:10])
set(gca, 'fontsize', 18)
axis([0.5, 10.5, -0.04, 0.2])
legend('Residuals rel', 'Residuals irrel', 'Initial state rel', 'Initial state irrel', 'Velocity rel', 'Velocity irrel', 'Velocity-Residual rel', 'Velocity-Residual irrel', 'Velocity-Initial state rel', 'Velocity-Initial state irrel', 'Residual-Initial state rel', 'Residual-Initial state irrel')
ylabel('Contribution to total correlation')
rotateXLabels( gca, 45)

%% Analytically calculate the contribution to selectivity arising from each source

var_method = 'fixedvar_rw';  % fixedvar for using same variance pooled both pre/post and with/without weights, fixedvar_rw for using variance from data either pre or post and not using variance with weights deleted
manipulation = 'removeW1W2';

clear SI
clear poolvar

tsims = 9:17;

W1 = 5; W2 = 4;

        [SI{W1,W2}, poolvar] = selectivity_analytic_switching(manipulation, var_method, W1, W2, xrel, xirrel, rmat_rel, rmat_irrel, ResidualV_rel_TOT, ResidualV_irrel_TOT, ResidualA_rel_TOT, ResidualA_irrel_TOT,vmatrel0, vmatirrel0, Trialsrel_all, Trialsirrel_all, tsims, rngtot, CELLLAB_ALL, RXS);

SI_Out_rel_tot = horzcat(SI_Out_rel{:});
SI_Out_irrel_tot = horzcat(SI_Out_irrel{:});

SI_model_rel_full = SI{W1,W2}.SIinputs(:,1) + SI{W1,W2}.SIinitstate(:,1) + SI{W1,W2}.SIvel(:,1) + SI{W1,W2}.SIres(:,1) ;
SI_model_irrel_full = SI{W1,W2}.SIinputs(:,2) + SI{W1,W2}.SIinitstate(:,2) + SI{W1,W2}.SIvel(:,2) + SI{W1,W2}.SIres(:,2) ;


SI_model_rel_full_rw = SI{W1,W2}.SIinputs_rw(:,1) + SI{W1,W2}.SIinitstate_rw(:,1) + SI{W1,W2}.SIvel_rw(:,1) + SI{W1,W2}.SIres_rw(:,1) ;
SI_model_irrel_full_rw = SI{W1,W2}.SIinputs_rw(:,2) + SI{W1,W2}.SIinitstate_rw(:,2) + SI{W1,W2}.SIvel_rw(:,2) + SI{W1,W2}.SIres_rw(:,2) ;

Method = 'noabs';

for j=1:4

    if strmatch(Method, 'abs')
    
        Model_rel_full(j) = mean(abs(SI_model_rel_full(CELLLAB_TOT == CLS(j))));
        Model_irrel_full(j) = mean(abs(SI_model_irrel_full(CELLLAB_TOT == CLS(j))));
        Data_rel(j) = mean(abs(SI_Out_rel_tot(CELLLAB_TOT == CLS(j))));
        Data_irrel(j) = mean(abs(SI_Out_irrel_tot(CELLLAB_TOT == CLS(j))));

        Model_rel_full_rw(j) = mean(abs(SI_model_rel_full_rw(CELLLAB_TOT == CLS(j))));
        Model_irrel_full_rw(j) = mean(abs(SI_model_irrel_full_rw(CELLLAB_TOT == CLS(j))));

    elseif strmatch(Method, 'noabs')
        
        Model_rel_full(j) = mean((SI_model_rel_full(CELLLAB_TOT == CLS(j))));
        Model_irrel_full(j) = mean((SI_model_irrel_full(CELLLAB_TOT == CLS(j))));
        Data_rel(j) = mean((SI_Out_rel_tot(CELLLAB_TOT == CLS(j))));
        Data_irrel(j) = mean((SI_Out_irrel_tot(CELLLAB_TOT == CLS(j))));

        Model_rel_full_rw(j) = mean((SI_model_rel_full_rw(CELLLAB_TOT == CLS(j))));
        Model_irrel_full_rw(j) = mean((SI_model_irrel_full_rw(CELLLAB_TOT == CLS(j))));

    end
    
end

bar([Data_rel - Data_irrel; Model_rel_full - Model_irrel_full; Model_rel_full_rw - Model_irrel_full_rw]')

CellNum = 5;

SIinputs_rw = SI{W1,W2}.SIinputs_rw;
SIinputs = SI{W1,W2}.SIinputs;

Xrel = [ones(size(SIinputs_rw(CELLLAB_TOT == CellNum, 1))), SIinputs_rw(CELLLAB_TOT == CellNum, 1)];
Xirrel = [ones(size(SIinputs_rw(CELLLAB_TOT == CellNum, 2))), SIinputs_rw(CELLLAB_TOT == CellNum, 2)];
Yrel = SIinputs(CELLLAB_TOT == CellNum, 1);
Yirrel = SIinputs(CELLLAB_TOT == CellNum, 2);
Wrel = Yrel.' / Xrel.';
Wirrel = Yirrel.'/ Xirrel.';

figure
subplot(121)
hold on
scatter(SIinputs_rw(CELLLAB_TOT == CellNum,1), SIinputs(CELLLAB_TOT == CellNum,1))
plot([-6, 6], [ -6, 6], 'k')
plot([-6,6], Wrel(2) * [-6,6] + Wrel(1))
axis([-6,6,-6,6])
set(gca, 'fontsize', 18)
title('SOM cells, rel learning')
xlabel('Selectivity of inputs without weights')
ylabel('Selectivity of inputs convolved with weights')
box on
subplot(122)
hold on
scatter(SIinputs_rw(CELLLAB_TOT == CellNum,2), SIinputs(CELLLAB_TOT == CellNum,2))
plot([-6, 6], [ -6, 6], 'k')
plot([-6,6], Wirrel(2) * [-6,6] + Wirrel(1))
axis([-6,6,-6,6])
set(gca, 'fontsize', 18)
title('SOM cells, irrel learning')
xlabel('Selectivity of inputs without weights')
ylabel('Selectivity of inputs convolved with weights')
box on


%% Test specificity of interaction weights between vertical and angled preferring cells



selective = 1;

for relirrel = {'rel', 'irrel'}

    Aopp = cell(4);
    Asame = cell(4);
    
    A_AV = cell(4);
    A_VA = cell(4);
    A_AA = cell(4);
    A_VV = cell(4);
    
    for CellNum1 = 1:4
        for CellNum2 = 1:4

            for rix = 1:length(RXS)

                if strmatch(relirrel, 'rel') 
                    
                    if selective
                    
                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1) & psel_rel{RXS(rix)}.' < 0.05);
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2) & psel_rel{RXS(rix)}.' < 0.05);

                    else

                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1));
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2));

                    end
                        
                    Arel = xrel{RXS(rix)}(:, 1:size(xirrel{RXS(rix)},1));
                    Apv = Arel(I1, I2);    
                    
                    SI = SI_In_rel{rix};
                %     SI = SI_Out_rel{rix};
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' > 0), bsxfun(@and, SI(I1)<0, SI(I2)'<0));
                Samepref_rel{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                SamePref_rel{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' < 0), bsxfun(@and, SI(I1)<0, SI(I2)'>0));
                OppPref_rel{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                OppPref_rel{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
   
                elseif strmatch(relirrel, 'irrel')
                    
                    if selective
                    
                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1) & psel_irrel{RXS(rix)}.' < 0.05);
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2) & psel_irrel{RXS(rix)}.' < 0.05);

                    else

                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1));
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2));

                    end
                                        
                    Airrel = xirrel{RXS(rix)}(:, 1:size(xirrel{RXS(rix)},1));
                    Apv = Airrel(I1, I2);
                    
                    SI = SI_In_irrel{rix};
                %    SI = SI_Out_irrel{rix};

                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' > 0), bsxfun(@and, SI(I1)<0, SI(I2)'<0));
                SamePref_irrel{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                SamePref_irrel{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' < 0), bsxfun(@and, SI(I1)<0, SI(I2)'>0));
                OppPref_irrel{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                OppPref_irrel{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
   
                    
                 end

                for i=1:length(Apv)

                    Apv(i,i) = nan;

                end

                dum = Apv(SI(I1) < 0, SI(I2) > 0);
                Aopp{CellNum1, CellNum2} = vertcat(Aopp{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) < 0);
                Aopp{CellNum1, CellNum2} = vertcat(Aopp{CellNum1, CellNum2}, dum(:));

                dum = Apv(SI(I1) < 0, SI(I2) < 0);
                Asame{CellNum1, CellNum2} = vertcat(Asame{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) > 0);
                Asame{CellNum1, CellNum2} = vertcat(Asame{CellNum1, CellNum2}, dum(:));

                % all 4 types
                
                dum = Apv(SI(I1) < 0, SI(I2) > 0);
                A_AV{CellNum1, CellNum2} = vertcat(A_AV{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) < 0);
                A_VA{CellNum1, CellNum2} = vertcat(A_VA{CellNum1, CellNum2}, dum(:));

                dum = Apv(SI(I1) < 0, SI(I2) < 0);
                A_AA{CellNum1, CellNum2} = vertcat(A_AA{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) > 0);
                A_VV{CellNum1, CellNum2} = vertcat(A_VV{CellNum1, CellNum2}, dum(:));
                
            end
            
        end
    end

    if strmatch(relirrel, 'rel')
        
        Asame_rel_tot = cellfun(@nanmean, Asame);
        Asame_rel_tot_sem = cellfun(@nanstd, Asame) ./ sqrt(cellfun(@length, Asame));
        Aopp_rel_tot = cellfun(@nanmean, Aopp);
        Aopp_rel_tot_sem = cellfun(@nanstd, Aopp) ./ sqrt(cellfun(@length, Aopp));
        
        A_AA_rel_tot = cellfun(@nanmean, A_AA);
        A_VV_rel_tot = cellfun(@nanmean, A_VV);
        A_AV_rel_tot = cellfun(@nanmean, A_AV);
        A_VA_rel_tot = cellfun(@nanmean, A_VA);
        
        for i=1:4
            for j=1:4
                
                  pval_rel(i,j) = ranksum(Asame{i,j}, Aopp{i,j});
                  
            end
        end
                            
    elseif strmatch(relirrel, 'irrel')
       
        Asame_irrel_tot = cellfun(@nanmean, Asame);
        Asame_irrel_tot_sem = cellfun(@nanstd, Asame) ./ sqrt(cellfun(@length, Asame));
        Aopp_irrel_tot = cellfun(@nanmean, Aopp);
        Aopp_irrel_tot_sem = cellfun(@nanstd, Aopp) ./ sqrt(cellfun(@length, Aopp));
        
        A_AA_irrel_tot = cellfun(@nanmean, A_AA);
        A_VV_irrel_tot = cellfun(@nanmean, A_VV);
        A_AV_irrel_tot = cellfun(@nanmean, A_AV);
        A_VA_irrel_tot = cellfun(@nanmean, A_VA);
                
        for i=1:4
            for j=1:4
                
                  pval_irrel(i,j) = ranksum(Asame{i,j}, Aopp{i,j});
                  
            end
        end
        
    end
    
end

Data_Means = [Asame_irrel_tot(:), Aopp_irrel_tot(:), Asame_rel_tot(:), Aopp_rel_tot(:)];
Data_SEM = [ Asame_irrel_tot_sem(:), Aopp_irrel_tot_sem(:), Asame_rel_tot_sem(:), Aopp_rel_tot_sem(:)];
hb = bar(1:size(Data_Means,1),Data_Means);
hold on

for ib = 1:numel(hb)
    
    xData = get(get(hb(ib),'Children'),'XData'); 
    barCenters = mean(unique(xData,'rows'));
    errorbar(barCenters,Data_Means(:,ib),Data_SEM(:,ib),'k.', 'linewidth', 2);
    
end

set(gca, 'fontsize', 18)
set(gca, 'xtick', 1:16)
set(gca, 'xticklabel', {'PYR->PYR', 'PYR->PV', 'PYR->SOM', 'PYR->VIP', 'PV->PYR', 'PV->PV', 'PV->SOM', 'PV->VIP', 'SOM->PYR', 'SOM->PV', 'SOM->SOM', 'SOM->VIP', 'VIP->PYR', 'VIP->PV', 'VIP->SOM', 'VIP->VIP'})
%rotateXLabels( gca, 45)
legend('Same preference, olfactory','Opposite preference, olfactory', 'Same preference, visual',  'Opposite preference, visual')
ylabel('Interaction weight')

H = rgb2hsv([1,0,0]);

set(hb(1), 'FaceColor', hsv2rgb([H(1), 0.3*H(2), H(3)]))
set(hb(3), 'FaceColor', hsv2rgb([H(1), 1*H(2), H(3)]))

H = rgb2hsv([0,0,1]);

set(hb(2), 'FaceColor', hsv2rgb([H(1), 0.3*H(2), H(3)]))
set(hb(4), 'FaceColor', hsv2rgb([H(1), 1*H(2), H(3)]))

%% Delete specific weights and calculate selectivity

% delete same pref PYR to PV weights

CELL1 = 2;
CELL2 = 1;

% PYR to PV weights

for rix = RXS
    
    xpre_rwsame{rix} = xpre{rix};
    xpre_rwsame{rix}(logical(SamePref_pre{CELL1,CELL2,rix})) = 0;

    xpost_rwsame{rix} = xpost{rix};
    xpost_rwsame{rix}(logical(SamePref_post{CELL1,CELL2,rix})) = 0;
    
end

[SI_rwsame, poolvar_rwsame] = selectivity_analytic_prepostshuff(xpre_rwsame, xpost_rwsame, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);


% delete opposite pref PYR to PV weights

% PYR to PV weights

for rix = RXS
    
    xpre_rwopp{rix} = xpre{rix};
    xpre_rwopp{rix}(logical(OppPref_pre{CELL1,CELL2,rix})) = 0;

    xpost_rwopp{rix} = xpost{rix};
    xpost_rwopp{rix}(logical(OppPref_post{CELL1,CELL2,rix})) = 0;
    
end

[SI_rwopp, poolvar_rwopp] = selectivity_analytic_prepostshuff(xpre_rwopp, xpost_rwopp, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);


SI_rwsame_pre_tot = SI_rwsame.SIinputs(:,1) + SI_rwsame.SIinitstate(:,1) + SI_rwsame.SIvel(:,1) + SI_rwsame.SIres(:,1);
SI_rwopp_pre_tot = SI_rwopp.SIinputs(:,1) + SI_rwopp.SIinitstate(:,1) + SI_rwopp.SIvel(:,1) + SI_rwopp.SIres(:,1);
SI_rwsame_post_tot = SI_rwsame.SIinputs(:,2) + SI_rwsame.SIinitstate(:,2) + SI_rwsame.SIvel(:,2) + SI_rwsame.SIres(:,2);
SI_rwopp_post_tot = SI_rwopp.SIinputs(:,2) + SI_rwopp.SIinitstate(:,2) + SI_rwopp.SIvel(:,2) + SI_rwopp.SIres(:,2);
SI_Out_pre_tot = horzcat(SI_Out_pre{:});
SI_Out_post_tot = horzcat(SI_Out_post{:});

CellType = 3;

MEAN = [mean(abs(SI_rwopp_pre_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), mean(abs(SI_rwopp_post_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_post_tot(CELLLAB_TOT == CellType))); mean(abs(SI_rwsame_pre_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), mean(abs(SI_rwsame_post_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_post_tot(CELLLAB_TOT == CellType)))];
SEM = [std(abs(SI_rwopp_pre_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), std(abs(SI_rwopp_post_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_post_tot(CELLLAB_TOT == CellType))); std(abs(SI_rwsame_pre_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), std(abs(SI_rwsame_post_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_post_tot(CELLLAB_TOT == CellType)))] ./ sqrt(length(CELLLAB_TOT == CellType));

hold on
hb = bar(1:size(MEAN,1), MEAN);

for ib = 1:numel(hb)
    
    xData = get(get(hb(ib),'Children'),'XData'); 
    barCenters = mean(unique(xData,'rows'));
    errorbar(barCenters,MEAN(:,ib),SEM(:,ib),'k.', 'linewidth', 2);
    
end

set(gca, 'fontsize', 18)
set(gca,'xtick', 1:2)
set(gca, 'xticklabel', {'Delete Opposite', 'Delete Same'})
ylabel('\Delta abs SI')
legend('Pre Learning', 'Post Learning')
title('Effect of PYR to PV weight deletion on PV selectivity')