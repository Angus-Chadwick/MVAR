%function rsbs_VR_forAngus

clear all

% This section finds all sessions and blocks in the data folder and stores
% the file names in idinf

ROOT = '/nfs/nhome/live/angus/Documents/Interneuron_Data/';

RSG.savedirforAngus = [ROOT,'Saved_Data/'];

list=dir([RSG.savedirforAngus,'/*_dat*']);
clear idinf;
for NID=1:length(list),
    ix1=strfind(list(NID).name,'_B');  % find unique sessions, day and block number
    ix2=strfind(list(NID).name,'_dat');
    list(NID).name(1:ix1-1);
    idinf{NID}.id    = list(NID).name(1:ix1-1);
    idinf{NID}.block = str2num(list(NID).name(ix1+2:ix2-1));
end

%--------------------------------------------------------------------------

% This section extracts and plots different quantities of interest from the
% raw data files

clc;

if ~exist('TSK'), TSK=[]; end

if ~isfield(TSK,'SW'),
    TSK.SW=load([RSG.savedirforAngus,'/SW_TLTPDPOP.mat']);  % load in metadata for cells and animals
end



for NID=1:length(idinf),  % loop over sessions/blocks
    
    % get additional information about this session from the extinf file
    
    clear EXTINF;  
    fname=sprintf('%s%s_B%d_extinf',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'EXTINF');
    nam    = EXTINF.nam;
    if strcmp(nam,'LN'), idinf{NID}.type='SD'; elseif strcmp(nam,'SW'), idinf{NID}.type='SWITCH'; end % determines if learning or switching
    RIXSES = EXTINF.RIXSES;  % used later for cell type classification
    EYE    = EXTINF.EYE;    
    EYEINF = EXTINF.EYEINF;
    
    clear out ADDINF;
    fname=sprintf('%s%s_B%d_dat',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
    load(fname,'out','ADDINF');
    
    
    clear ADDTRG

 fname=sprintf('%s%s%s_B%d',RSG.savedirforAngus, 'ADDTRG2_', idinf{NID}.id,idinf{NID}.block);
 load(fname);
 
 
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
    
%% Extract mean response and selectivities, and reject outliers from normal distribution

            
                c1=strmatch('Von',ADDTRG{1}.P.BTRGLABEL, 'exact'); % condition 1: Vertical
                c2=strmatch('Aon',ADDTRG{1}.P.BTRGLABEL, 'exact'); % condition 2: Angled
                c3=strmatch('VonIrr',ADDTRG{1}.P.BTRGLABEL, 'exact'); % condition 1: Vertical
                c4=strmatch('AonIrr',ADDTRG{1}.P.BTRGLABEL, 'exact'); % condition 2: Angled
                m_vert_rel = ADDTRG{1}.PL{1}.TwinPP.m1(:,:,c1);
                m_ang_rel = ADDTRG{1}.PL{1}.TwinPP.m1(:,:,c2);
                m_vert_irrel = ADDTRG{1}.PL{1}.TwinPP.m1(:,:,c3); 
                m_ang_irrel  = ADDTRG{1}.PL{1}.TwinPP.m1(:,:,c4);
                v_vert_rel   = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,c1);
                v_ang_rel    = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,c2);
                v_vert_irrel = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,c3);
                v_ang_irrel  = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,c4);
                n_vert_rel = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,c1);
                n_ang_rel  = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,c2);
                n_vert_irrel = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,c3);
                n_ang_irrel  = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,c4);

                poolstd_rel   = sqrt( (v_vert_rel   .* (n_vert_rel - 1)   + v_ang_rel   .* (n_ang_rel - 1))   ./ (n_vert_rel   + n_ang_rel   - 2) );
                poolstd_irrel = sqrt( (v_vert_irrel .* (n_vert_irrel - 1) + v_ang_irrel .* (n_ang_irrel - 1)) ./ (n_vert_irrel + n_ang_irrel - 2) );

                Srel_all   = (m_vert_rel   - m_ang_rel)   ./ poolstd_rel;
                Sirrel_all = (m_vert_irrel - m_ang_irrel) ./ poolstd_irrel;

                
Nstd = 5; % number of standard deviations from mean to tolerate 

SI = [Sirrel_all; Srel_all].'; 
SIPYR = SI(CELLLAB == 1, :);

SIG = cov(SIPYR);  
MU  = mean(SIPYR);

MD = sqrt(diag((SI.' - repmat(MU.', 1, size(SI, 1))).' * inv(SIG) * (SI.' - repmat(MU.', 1, size(SI, 1))))); % Malahabinois distance of each cell

reject = and(MD > Nstd, CELLLAB==1); % reject only pyramidal cells whose selectivity is outwith tolerance bounds
accept = ~reject;


Srel_all = Srel_all(accept);
Sirrel_all = Sirrel_all(accept);

CELLLAB = CELLLAB(accept);

%%%
%%%

Srel = Srel_all(CELLLAB == 1);
Sirrel = Sirrel_all(CELLLAB == 1);

SrelPV = Srel_all(CELLLAB == 3);
SirrelPV = Sirrel_all(CELLLAB == 3);

SrelSOM = Srel_all(CELLLAB == 4);
SirrelSOM = Sirrel_all(CELLLAB == 4);

SrelVIP = Srel_all(CELLLAB == 5);
SirrelVIP = Sirrel_all(CELLLAB == 5);

dS = Srel_all - Sirrel_all;
dS_PYR = dS(CELLLAB == 1);
dS_PV = dS(CELLLAB == 3);
dS_SOM = dS(CELLLAB == 4);
dS_VIP = dS(CELLLAB == 5);

                
 %% get noise correlation spectrum in relevant and irrelevant conditions

rix = NID;

dF = (m_vert_rel + m_ang_rel - (m_vert_irrel + m_ang_irrel));
dF = dF(accept);
dFPYR = dF(CELLLAB == 1); % Fractional change in PYR firing between conditions
dFPV = dF(CELLLAB == 3); % Fractional change in PV firing between conditions
dFSOM = dF(CELLLAB == 4); % Fractional change in SOM firing between conditions
dFVIP = dF(CELLLAB == 5); % Fractional change in VIP firing between conditions

COND = 11;

MATrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MATrel(isnan(MATrel)) = 0;
MATrel = MATrel + triu(MATrel,1).';  %% symmetrise matrix
MATrel = MATrel(accept, accept);
MATrel = MATrel + eye(size(MATrel,1));


COND = 11;

MATirrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MATirrel(isnan(MATirrel)) = 0;
MATirrel = MATirrel + triu(MATirrel,1).';  %% symmetrise matrix
MATirrel = MATirrel(accept, accept);
MATirrel = MATirrel + eye(size(MATirrel,1));

rix = NID;

CELLS = [1,3,4,5];

[V, D] = eig(MATrel); % eigenspectrum of noise correlation matrix

for i=1:size(V,2)  
    for C = 1:4
           
    [ModeCorr{rix}(i,C), pval{rix}(i,C)] = corr(Srel_all(CELLLAB==CELLS(C)).', V(CELLLAB==CELLS(C),i));
    [ModeCorrabs{rix}(i,C), pvalabs{rix}(i,C)] = corr(abs(Srel_all(CELLLAB==CELLS(C))).', V(CELLLAB==CELLS(C),i));
       
    end
    
        V(:,i) = V(:,i) * sign(ModeCorr{rix}(i,1));
    
end
 
%      alpha = 0.01;
% 
%      selmodes = find((pval{rix}(:,1)) < alpha); % modes which determine selectivity
%      selmodes_corr = ModeCorr{rix}(selmodes,1);
%  
   %  [Y,I] = sort(abs(selmodes_corr), 'descend');
    [Y,I] = sort(pval{rix}, 'ascend');
    
%    I = 1:size(V,2);

    lambda = diag(D);
    
%% Implement Cross-validated Prediction Model

clear MODE1_TRAIN
clear MODE1_TEST
clear Srel_TRAIN
clear Sirrel_TRAIN

NPYR = sum(CELLLAB == 1);

for i=1:NPYR
     
    ind0 = find(CELLLAB == 1);
    ind1 = [1:(i-1), (i+1):NPYR];
    ind  = ind0(ind1); 
       
    ModeCorr0{rix} = zeros([size(V,1), 1]);
    pval0{rix} = zeros([size(V,1), 1]);
    
    for j=1:size(V,2)  
    
     [ModeCorr0{rix}(j,1), pval0{rix}(j,1)] = corr(Srel_all(ind).', V(ind,j));

    end 
    
    [Y0{i},I0{i}] = sort(pval0{rix}, 'ascend');

end


    
for Nmodes = 2:175
    
    MODE1_TRAIN = zeros([NPYR, NPYR - 1, Nmodes]);
    MODE1_TEST = zeros([NPYR, Nmodes]);
    MODESUM_TRAIN = zeros([NPYR, NPYR - 1]);
    MODESUM_TEST = zeros([NPYR, 1]);
    Srel_TRAIN = zeros([1, NPYR - 1]);
    Sirrel_TRAIN = zeros([1, NPYR - 1]);

for i=1:NPYR
    
    % training datasets - remove one pyramidal cell and compute prediction
    % metrics
    
    ind0 = find(CELLLAB == 1);
    ind1 = [1:(i-1), (i+1):NPYR];
    ind  = ind0(ind1); 
     
    MODE1_TRAIN(i,:,:) = V(ind,I0{i}(1:Nmodes));
    MODESUM_TRAIN(i,:) = V(ind, I0{i}(1:Nmodes)) * lambda(I0{i}(1:Nmodes));

    Srel_TRAIN(i,:) = Srel(ind1);
    Sirrel_TRAIN(i,:) = Sirrel(ind1);
    
    % test datasets - use all peer cells for prediction
          
    MODE1_TEST(i,:) = V(ind0(i), I0{i}(1:Nmodes));
    MODESUM_TEST(i) = V(ind0(i), I0{i}(1:Nmodes)) * lambda(I0{i}(1:Nmodes));  
 
end




METHOD = 2;  % 1 for 10 fold cross validation and 2 for leave-one-out cross validation

if METHOD == 1 % if N fold cross validation

elseif METHOD == 2 % leave-one-out cross validation
    
     NPYR = sum(CELLLAB == 1);
 
    clear condmean_rel_peers
    clear condmean_rel_sum
    clear condmean_rel_null
    clear condcov_rel_peers
    clear condcov_rel_sum
    clear condcov_rel_null

    clear LNULL
    clear LTEST 
    clear LSUM

    CVAR = zeros(6);
    
    condmean_rel_peers = zeros([1, NPYR]);
    condmean_rel_null  = zeros([1, NPYR]);
    condmean_rel_sum   = zeros([1, NPYR]);
    condcov_rel_peers = zeros([1, NPYR]);
    condcov_rel_sum   = zeros([1, NPYR]);
    condcov_rel_null  = zeros([1, NPYR]);
    LTEST = zeros([1, NPYR]);
    LNULL = zeros([1, NPYR]);
    LSUM  = zeros([1, NPYR]);
   
    for i=1:NPYR
            
       DATA_TEST  = [Srel(i); MODE1_TEST(i,:).'; Sirrel(i)].';
       DATA_TRAIN = [Srel_TRAIN(i,:); squeeze(MODE1_TRAIN(i,:,:)).'; Sirrel_TRAIN(i,:)].';
   
       DATA_TEST_SUM  = [Srel(i); MODESUM_TEST(i,:).'; Sirrel(i)].';
       DATA_TRAIN_SUM = [Srel_TRAIN(i,:); squeeze(MODESUM_TRAIN(i,:)); Sirrel_TRAIN(i,:)].';
       
       DATA_TEST_NULL  = [Srel(i); Sirrel(i)].';
       DATA_TRAIN_NULL = [Srel_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
  
       meansel = mean(DATA_TRAIN,1);
       meansel_sum = mean(DATA_TRAIN_SUM,1);
       covsel = cov(DATA_TRAIN);
       covsel_sum = cov(DATA_TRAIN_SUM);
    
       meansel_null = mean(DATA_TRAIN_NULL,1);
       covsel_null = cov(DATA_TRAIN_NULL);

       condmean_rel_peers(i) = meansel(1)      + (covsel(1, 2:end)      * inv(covsel(2:end, 2:end)))      * (DATA_TEST(1,2:end)      - meansel(2:end)).';
       condmean_rel_sum(i) = meansel_sum(1)    + (covsel_sum(1, 2:end)  * inv(covsel_sum(2:end, 2:end)))  * (DATA_TEST_SUM(1,2:end)  - meansel_sum(2:end)).';
       condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
       
       condcov_rel_peers(i) = covsel(1,1) - covsel(1,2:end) * inv(covsel(2:end,2:end)) * covsel(2:end,1);
       condcov_rel_sum(i) = covsel_sum(1,1) - covsel_sum(1,2:end) * inv(covsel_sum(2:end,2:end)) * covsel_sum(2:end,1);
       condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
       
       LTEST(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_peers(i).')) - 1/2 * (DATA_TEST(:,1)      - condmean_rel_peers(i).').' * (condcov_rel_peers(i).').^(-1) * (DATA_TEST(:,1)      - condmean_rel_peers(i).'));
       LSUM(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_sum(i).')) - 1/2 * (DATA_TEST_SUM(:,1)      - condmean_rel_sum(i).').' * (condcov_rel_sum(i).').^(-1) * (DATA_TEST_SUM(:,1)      - condmean_rel_sum(i).'));
       LNULL(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
        
       
    end  
end

SHUFFLE = 0;

if SHUFFLE

    pred_shuffle{rix} = condmean_rel_peers;
    
else
    
    pred_peers{rix, Nmodes} = condmean_rel_peers;
    pred_sum{rix, Nmodes}   = condmean_rel_sum;

    
end

pred_null{rix, Nmodes} = condmean_rel_null;
FracImproved(rix, Nmodes) = sum(LTEST > LNULL) / length(LTEST);

end

Sobs{rix} = Srel;

end

%% Generate 4d Scatter plots

subplot(221)
hold on
scatter(V(CELLLAB == 1, I(1)), V(CELLLAB == 1, I(2)), 500, Srel_all(CELLLAB == 1).', '.', 'linewidth', 2)
scatter(V(CELLLAB == 3, I(1)), V(CELLLAB == 3, I(2)), 500, Srel_all(CELLLAB == 3), '*', 'linewidth', 2)
scatter(V(CELLLAB == 4, I(1)), V(CELLLAB == 4, I(2)), 500, Srel_all(CELLLAB == 4).', 'o', 'linewidth', 2)
scatter(V(CELLLAB == 5, I(1)), V(CELLLAB == 5, I(2)), 500, Srel_all(CELLLAB == 5), '+', 'linewidth', 2)
set(gca, 'Clim', [min(Srel_all), max(Srel_all)])
colorbar

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 2')
box on


subplot(222)
hold on
scatter(V(CELLLAB == 1, I(3)), V(CELLLAB == 1, I(2)), 500, Srel_all(CELLLAB == 1).', '.', 'linewidth', 2)
scatter(V(CELLLAB == 3, I(3)), V(CELLLAB == 3, I(2)), 500, Srel_all(CELLLAB == 3), '*', 'linewidth', 2)
scatter(V(CELLLAB == 4, I(3)), V(CELLLAB == 4, I(2)), 500, Srel_all(CELLLAB == 4).', 'o', 'linewidth', 2)
scatter(V(CELLLAB == 5, I(3)), V(CELLLAB == 5, I(2)), 500, Srel_all(CELLLAB == 5), '+', 'linewidth', 2)
set(gca, 'Clim', [min(Srel_all), max(Srel_all)])
colorbar

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 2')
box on

subplot(223)
hold on
scatter(V(CELLLAB == 1, I(1)), V(CELLLAB == 1, I(4)), 500, Srel_all(CELLLAB == 1).', '.', 'linewidth', 2)
scatter(V(CELLLAB == 3, I(1)), V(CELLLAB == 3, I(4)), 500, Srel_all(CELLLAB == 3), '*', 'linewidth', 2)
scatter(V(CELLLAB == 4, I(1)), V(CELLLAB == 4, I(4)), 500, Srel_all(CELLLAB == 4).', 'o', 'linewidth', 2)
scatter(V(CELLLAB == 5, I(1)), V(CELLLAB == 5, I(4)), 500, Srel_all(CELLLAB == 5), '+', 'linewidth', 2)
set(gca, 'Clim', [min(Srel_all), max(Srel_all)])
colorbar

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 4')
box on

subplot(224)
hold on
scatter(V(CELLLAB == 1, I(3)), V(CELLLAB == 1, I(4)), 500, Srel_all(CELLLAB == 1).', '.', 'linewidth', 2)
scatter(V(CELLLAB == 3, I(3)), V(CELLLAB == 3, I(4)), 500, Srel_all(CELLLAB == 3), '*', 'linewidth', 2)
scatter(V(CELLLAB == 4, I(3)), V(CELLLAB == 4, I(4)), 500, Srel_all(CELLLAB == 4).', 'o', 'linewidth', 2)
scatter(V(CELLLAB == 5, I(3)), V(CELLLAB == 5, I(4)), 500, Srel_all(CELLLAB == 5), '+', 'linewidth', 2)
set(gca, 'Clim', [min(Srel_all), max(Srel_all)])
colorbar

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 4')
box on

% for i=1:100
%     
%     
%     for j=1:length(CELLLAB)
%     
%                 D{i}(j,:) = V(j, I(i)) - V(:, I(i));
%                 
%     end
%         
% end



