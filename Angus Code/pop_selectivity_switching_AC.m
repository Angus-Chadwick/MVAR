
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

%% Calculate mahalanobis distance etc

CellNum = 1

MATrelV = ADDTRG{1}.PL{1}.DCOR{1}.CORNSE1;  % noise correlations for grey corridor condition
MATrelA = ADDTRG{1}.PL{1}.DCOR{2}.CORNSE1;  % noise correlations for grey corridor condition
MATirrelV = ADDTRG{1}.PL{1}.DCOR{5}.CORNSE1;  % noise correlations for grey corridor condition
MATirrelA = ADDTRG{1}.PL{1}.DCOR{6}.CORNSE1;  % noise correlations for grey corridor condition

MATrelV(isnan(MATrelV)) = 0;
MATrelV = MATrelV + triu(MATrelV, 1).';
MATrelV = MATrelV(accept, accept);
MATrelV(find(eye(size(MATrelV, 1)))) = 1;

MATirrelV(isnan(MATirrelV)) = 0;
MATirrelV = MATirrelV + triu(MATirrelV, 1).';
MATirrelV = MATirrelV(accept, accept);
MATirrelV(find(eye(size(MATirrelV, 1)))) = 1;

MATrelA(isnan(MATrelA)) = 0;
MATrelA = MATrelA + triu(MATrelA, 1).';
MATrelA = MATrelA(accept, accept);
MATrelA(find(eye(size(MATrelA, 1)))) = 1;

MATirrelA(isnan(MATirrelA)) = 0;
MATirrelA = MATirrelA + triu(MATirrelA, 1).';
MATirrelA = MATirrelA(accept, accept);
MATirrelA(find(eye(size(MATirrelA, 1)))) = 1;


%% True correlation matrix

rix = NID;

CellType = 1;

sigmaVrel = sqrt(v_vert_rel);
sigmaArel = sqrt(v_ang_rel);

sigmaVirrel = sqrt(v_vert_irrel);
sigmaAirrel = sqrt(v_ang_irrel);

CVrelA = zeros(size(MATrelV));
CVrelV = zeros(size(MATrelV));
CVirrelA = zeros(size(MATirrelV));
CVirrelV = zeros(size(MATirrelV));

for i=1:size(MATrelV)
    for j=1:size(MATirrelV)
       
        CVrelA(i,j) = sigmaArel(i)  * sigmaArel(j);
        CVrelV(i,j) =  sigmaVrel(i)  * sigmaVrel(j);
        
        CVirrelA(i,j) = sigmaAirrel(i)  * sigmaAirrel(j);
        CVirrelV(i,j) =  sigmaVirrel(i)  * sigmaVirrel(j);
        
                
    end
end

CovrelA = CVrelA .* MATrelA;
CovrelV = CVrelV .* MATrelV;
CovirrelA = CVirrelA .* MATirrelA;
CovirrelV = CVirrelV .* MATirrelV;

Mtotrel = corrcov(1/2*(CovrelA + CovrelV));
Mtotirrel = corrcov(1/2*(CovirrelA + CovirrelV));

Sirrel = Sirrel_all(CELLLAB == CellType);
Srel = Srel_all(CELLLAB == CellType);

Mtotrel_pyr = Mtotrel(CELLLAB == CellType, CELLLAB == CellType);
Mtotirrel_pyr = Mtotirrel(CELLLAB == CellType, CELLLAB == CellType);

Minvrel = pinv(Mtotrel_pyr);
%Minvrel = (Minvrel + Minvrel.')/2;
Minvirrel = pinv(Mtotirrel_pyr);
%Minvirrel = (Minvirrel + Minvirrel.')/2;

 R = randperm(length(Sirrel));
 R = R(1:100);
%R = 1:length(Sirrel);

Dmirrelp(rix) = 1/4 * Sirrel * Minvirrel * Sirrel.';
Dmrelp(rix) = 1/4 * Srel * Minvrel * Srel.';

Dmirrel(rix) = 1/4 * Sirrel(R) * inv(Mtotirrel_pyr(R,R)) * Sirrel(R).';
Dmrel(rix) = 1/4 * Srel(R) * inv(Mtotrel_pyr(R,R)) * Srel(R).';

Dmirrel_diag(rix) = 1/4 * Sirrel(R) * Sirrel(R).';
Dmrel_diag(rix)  = 1/4 * Srel(R) * Srel(R).';

Dmirrel_normed(rix) = Dmirrel(rix) / norm(Sirrel(R))^2;
Dmrel_normed(rix) = Dmrel(rix) / norm(Srel(R))^2;

[Vrel, Drel] = eig(Mtotrel_pyr);
[Virrel, Dirrel] = eig(Mtotirrel_pyr);

lambdarel{rix} = diag(Drel);
lambdairrel{rix} = diag(Dirrel);

for i=1:length(diag(Drel))

    Mrel{rix}(i) = 1/4 * (dot(Vrel(:,i), Srel))^2 / lambdarel{rix}(i);
    Mirrel{rix}(i) = 1/4 * (dot(Virrel(:,i), Sirrel))^2 / lambdairrel{rix}(i);
    
    thetarel{rix}(i) = acos(dot(Vrel(:,i), Srel) / norm(Srel));
    thetairrel{rix}(i) = acos(dot(Virrel(:,i), Sirrel) / norm(Sirrel));

    Mrel_cum{rix} = cumsum(sort(Mrel{rix}, 'descend'));
    Mirrel_cum{rix} = cumsum(sort(Mirrel{rix}, 'descend'));
    
end

sumvarrel(rix) = sum(sigmaVrel.^2 + sigmaArel.^2);
sumvarirrel(rix) = sum(sigmaVirrel.^2 + sigmaAirrel.^2);


for i=1:5
    for j=1:5
        
        FracNegrel(i,j,rix) = sum(sum(Mtotrel(CELLLAB == i, CELLLAB == j) < 0))./(sum(CELLLAB == i) * sum(CELLLAB == j));
        FracNegirrel(i,j,rix) = sum(sum(Mtotirrel(CELLLAB == i, CELLLAB == j) < 0))./(sum(CELLLAB == i) * sum(CELLLAB == j));

    end
end

M(rix) = max(max(inv(Mtotrel)));
DD(rix) = det(Mtotrel);
NN(rix) = size(Mtotrel, 1);

EErel(rix) = mean(mean(abs(inv(Mtotrel_pyr) * Mtotrel_pyr - eye(size(Mtotrel_pyr,1)))));
EEirrel(rix) = mean(mean(abs(inv(Mtotirrel_pyr) * Mtotirrel_pyr - eye(size(Mtotrel_pyr,1)))));


    end
    



