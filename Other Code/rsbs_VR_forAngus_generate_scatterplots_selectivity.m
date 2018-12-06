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

%% Get Noise Correlations

rix = NID;

COND = 1;

MATrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MATrel(isnan(MATrel)) = 0;
MATrel = MATrel + triu(MATrel,1).';  %% symmetrise matrix
MATrel = MATrel + eye(size(MATrel,1));


COND = 5;

MATirrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MATirrel(isnan(MATirrel)) = 0;
MATirrel = MATirrel + triu(MATirrel,1).';  %% symmetrise matrix
MATirrel = MATirrel + eye(size(MATirrel,1));

%% Extract mean response and selectivities

            
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


Spre = Srel_all(CELLLAB == 1);
Spost = Sirrel_all(CELLLAB == 1);

Spre_PV = Srel_all(CELLLAB == 3);
Spost_PV = Sirrel_all(CELLLAB == 3);

Spre_SOM = Srel_all(CELLLAB == 4);
Spost_SOM = Sirrel_all(CELLLAB == 4);

Spre_VIP = Srel_all(CELLLAB == 5);
Spost_VIP = Sirrel_all(CELLLAB == 5);               
                
%% Get Eigenspectra

% Pre Learning

[V, D] = eig(MATirrel); % eigenspectrum of noise correlation matrix


for i=1:size(V,1)
    
    if median(V(CELLLAB==1,i)) < 0
        
        V(:,i) = -V(:,i);
        
    end
end

figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')

subplot(221)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 1), 100 * abs(Spre), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 1), 100 * abs(Spre_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 1), 100 * abs(Spre_SOM),'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 1), 100 * abs(Spre_VIP), 'k','+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 2')
box on

subplot(222)

hold on

scatter(V(CELLLAB0 == 1, end-2), V(CELLLAB0 == 1, end - 1), 100 * abs(Spre), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end-2), V(CELLLAB0 == 3, end - 1), 100 * abs(Spre_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end-2), V(CELLLAB0 == 4, end - 1), 100 * abs(Spre_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end-2), V(CELLLAB0 == 5, end - 1), 100 * abs(Spre_VIP),' k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 2')
box on


subplot(223)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 3), 100 * abs(Spre), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 3), 100 * abs(Spre_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 3), 100 * abs(Spre_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 3), 100 * abs(Spre_VIP), 'k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 4')
box on


subplot(224)

hold on
 
scatter(V(CELLLAB0 == 1, end - 2), V(CELLLAB0 == 1, end-3), 100 * abs(Spre), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end - 2), V(CELLLAB0 == 3, end-3), 100 * abs(Spre_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end - 2), V(CELLLAB0 == 4, end-3), 100 * abs(Spre_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end - 2), V(CELLLAB0 == 5, end-3), 100 * abs(Spre_VIP), 'k', '+', 'linewidth', 2)
 
set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 4')
box on

set (gcf, 'paperpositionmode', 'manual','paperposition',[0 0 75 50])
print(strcat(idinf{rix}.id, '_switch_irrel'),'-dpng','-r0')

close(gcf)

% Post Learning

[V, D] = eig(MATrel); % eigenspectrum of noise correlation matrix


for i=1:size(V,1)
    
    if median(V(CELLLAB==1,i)) < 0
        
        V(:,i) = -V(:,i);
        
    end
end

figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')

subplot(221)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 1), 100 * abs(Spost), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 1), 100 * abs(Spost_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 1), 100 * abs(Spost_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 1), 100 * abs(Spost_VIP), 'k','+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 2')
box on

subplot(222)

hold on

scatter(V(CELLLAB0 == 1, end-2), V(CELLLAB0 == 1, end - 1), 100 * abs(Spost), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end-2), V(CELLLAB0 == 3, end - 1), 100 * abs(Spost_PV) , '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end-2), V(CELLLAB0 == 4, end - 1), 100 * abs(Spost_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end-2), V(CELLLAB0 == 5, end - 1), 100 * abs(Spost_VIP), 'k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
xlabel('Mode 3')
ylabel('Mode 2')
box on


subplot(223)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 3), 100 * abs(Spost), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 3), 100 * abs(Spost_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 3), 100 * abs(Spost_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 3), 100 * abs(Spost_VIP), 'k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
xlabel('Mode 1')
ylabel('Mode 4')
box on


subplot(224)

hold on
 
scatter(V(CELLLAB0 == 1, end - 2), V(CELLLAB0 == 1, end-3), 100 * abs(Spost), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end - 2), V(CELLLAB0 == 3, end-3), 100 * abs(Spost_PV), '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end - 2), V(CELLLAB0 == 4, end-3), 100 * abs(Spost_SOM), 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end - 2), V(CELLLAB0 == 5, end-3), 100 * abs(Spost_VIP), 'k', '+', 'linewidth', 2)
 
set(gca, 'fontsize', 18)
xlabel('Mode 3')
ylabel('Mode 4')
box on

set (gcf, 'paperpositionmode', 'manual','paperposition',[0 0 75 50])
print(strcat(idinf{rix}.id, '_switch_rel'),'-dpng','-r0')

close(gcf)


end
