%function rsbs_VR_forAngus

clear all

for Ni = 1:1

    N = Ni;
    Ndim(Ni) = N;
    
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
    
clear ADDTRG

 fname=sprintf('%s%s%s_B%d',RSG.savedirforAngus, 'ADDTRG2_', idinf{NID}.id,idinf{NID}.block);
 load(fname);
 
 %% get noise correlation spectrum in relevant and irrelevant conditions

rix = NID;
 
COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
%MAT = (MAT + MAT.')/2;  %% symmetrise matrix
MAT = MAT + triu(MAT, 1).';

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

for i=1:size(V,1)
    
    if median(V(CELLLAB0==1,i)) < 0
        
        V(:,i) = -V(:,i);
        
    end
end


% PYR-PYR

QrelV = zeros(length(CELLLAB));
Q0relV = zeros(length(CELLLAB));


for i=1:length(CELLLAB)
    for j=1:length(CELLLAB)
        
QrelV(i,j) = norm(V(i, (end - (N-1)):end) - V(j, (end - (N-1)):end));
Q0relV(i,j) = abs(V(i, Ni) - V(j, Ni)); 

    end
end



QrelPYRPYR{rix} = QrelV(CELLLAB==1, CELLLAB==1);
QrelPYRPYR{rix}(QrelPYRPYR{rix} == 0) = [];
QrelPYRPV{rix} = QrelV(CELLLAB==1, CELLLAB==3);
QrelPYRPV{rix} = QrelPYRPV{rix}(:);
QrelPYRSOM{rix} = QrelV(CELLLAB==1, CELLLAB==4);
QrelPYRSOM{rix} = QrelPYRSOM{rix}(:);
QrelPYRVIP{rix} = QrelV(CELLLAB==1, CELLLAB==5);
QrelPYRVIP{rix} = QrelPYRVIP{rix}(:);
QrelPVPV{rix} = QrelV(CELLLAB==3, CELLLAB==3);
QrelPVPV{rix}(QrelPVPV{rix} == 0) = [];
QrelPVSOM{rix} = QrelV(CELLLAB==3, CELLLAB==4);
QrelPVSOM{rix} = QrelPVSOM{rix}(:);
QrelPVVIP{rix} = QrelV(CELLLAB==3, CELLLAB==5);
QrelPVVIP{rix} = QrelPVVIP{rix}(:);
QrelSOMSOM{rix} = QrelV(CELLLAB==4, CELLLAB==4);
QrelSOMSOM{rix}(QrelSOMSOM{rix} == 0) = [];
QrelSOMVIP{rix} = QrelV(CELLLAB==4, CELLLAB==5);
QrelSOMVIP{rix} = QrelSOMVIP{rix}(:);
QrelVIPVIP{rix} = QrelV(CELLLAB==5, CELLLAB==5);
QrelVIPVIP{rix}(QrelVIPVIP{rix} == 0) = [];

Q0relPYRPYR{rix} = Q0relV(CELLLAB==1, CELLLAB==1);
Q0relPYRPYR{rix}(Q0relPYRPYR{rix} == 0) = [];
Q0relPYRPV{rix} = Q0relV(CELLLAB==1, CELLLAB==3);
Q0relPYRPV{rix} = Q0relPYRPV{rix}(:);
Q0relPYRSOM{rix} = Q0relV(CELLLAB==1, CELLLAB==4);
Q0relPYRSOM{rix} = Q0relPYRSOM{rix}(:);
Q0relPYRVIP{rix} = Q0relV(CELLLAB==1, CELLLAB==5);
Q0relPYRVIP{rix} = Q0relPYRVIP{rix}(:);
Q0relPVPV{rix} = Q0relV(CELLLAB==3, CELLLAB==3);
Q0relPVPV{rix}(Q0relPVPV{rix} == 0) = [];
Q0relPVSOM{rix} = Q0relV(CELLLAB==3, CELLLAB==4);
Q0relPVSOM{rix} = Q0relPVSOM{rix}(:);
Q0relPVVIP{rix} = Q0relV(CELLLAB==3, CELLLAB==5);
Q0relPVVIP{rix} = Q0relPVVIP{rix}(:);
Q0relSOMSOM{rix} = Q0relV(CELLLAB==4, CELLLAB==4);
Q0relSOMSOM{rix}(Q0relSOMSOM{rix} == 0) = [];
Q0relSOMVIP{rix} = Q0relV(CELLLAB==4, CELLLAB==5);
Q0relSOMVIP{rix} = Q0relSOMVIP{rix}(:);
Q0relVIPVIP{rix} = Q0relV(CELLLAB==5, CELLLAB==5);
Q0relVIPVIP{rix}(Q0relVIPVIP{rix} == 0) = [];


% task-irrelevant

COND = 5; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
%MAT = (MAT + MAT.')/2;  %% symmetrise matrix
MAT = MAT + triu(MAT, 1).';

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

for i=1:size(V,1)
    
    if median(V(CELLLAB0==1,i)) < 0
        
        V(:,i) = -V(:,i);
        
    end
end



[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

% PYR-PYR

QirrelV = zeros(length(CELLLAB));
Q0irrelV = zeros(length(CELLLAB));


for i=1:length(CELLLAB)
    for j=1:length(CELLLAB)
        
QirrelV(i,j) = norm(V(i, (end - (N-1)):end) - V(j, (end - (N-1)):end));
Q0irrelV(i,j) = abs(V(i,Ni) - V(j,Ni)); 

    end
end


QirrelPYRPYR{rix} = QirrelV(CELLLAB==1, CELLLAB==1);
QirrelPYRPYR{rix}(QirrelPYRPYR{rix} == 0) = [];
QirrelPYRPV{rix} = QirrelV(CELLLAB==1, CELLLAB==3);
QirrelPYRPV{rix} = QirrelPYRPV{rix}(:);
QirrelPYRSOM{rix} = QirrelV(CELLLAB==1, CELLLAB==4);
QirrelPYRSOM{rix} = QirrelPYRSOM{rix}(:);
QirrelPYRVIP{rix} = QirrelV(CELLLAB==1, CELLLAB==5);
QirrelPYRVIP{rix} = QirrelPYRVIP{rix}(:);
QirrelPVPV{rix} = QirrelV(CELLLAB==3, CELLLAB==3);
QirrelPVPV{rix}(QirrelPVPV{rix} == 0) = [];
QirrelPVSOM{rix} = QirrelV(CELLLAB==3, CELLLAB==4);
QirrelPVSOM{rix} = QirrelPVSOM{rix}(:);
QirrelPVVIP{rix} = QirrelV(CELLLAB==3, CELLLAB==5);
QirrelPVVIP{rix} = QirrelPVVIP{rix}(:);
QirrelSOMSOM{rix} = QirrelV(CELLLAB==4, CELLLAB==4);
QirrelSOMSOM{rix}(QirrelSOMSOM{rix} == 0) = [];
QirrelSOMVIP{rix} = QirrelV(CELLLAB==4, CELLLAB==5);
QirrelSOMVIP{rix} = QirrelSOMVIP{rix}(:);
QirrelVIPVIP{rix} = QirrelV(CELLLAB==5, CELLLAB==5);
QirrelVIPVIP{rix}(QirrelVIPVIP{rix} == 0) = [];

Q0irrelPYRPYR{rix} = Q0irrelV(CELLLAB==1, CELLLAB==1);
Q0irrelPYRPYR{rix}(Q0irrelPYRPYR{rix} == 0) = [];
Q0irrelPYRPV{rix} = Q0irrelV(CELLLAB==1, CELLLAB==3);
Q0irrelPYRPV{rix} = Q0irrelPYRPV{rix}(:);
Q0irrelPYRSOM{rix} = Q0irrelV(CELLLAB==1, CELLLAB==4);
Q0irrelPYRSOM{rix} = Q0irrelPYRSOM{rix}(:);
Q0irrelPYRVIP{rix} = Q0irrelV(CELLLAB==1, CELLLAB==5);
Q0irrelPYRVIP{rix} = Q0irrelPYRVIP{rix}(:);
Q0irrelPVPV{rix} = Q0irrelV(CELLLAB==3, CELLLAB==3);
Q0irrelPVPV{rix}(Q0irrelPVPV{rix} == 0) = [];
Q0irrelPVSOM{rix} = Q0irrelV(CELLLAB==3, CELLLAB==4);
Q0irrelPVSOM{rix} = Q0irrelPVSOM{rix}(:);
Q0irrelPVVIP{rix} = Q0irrelV(CELLLAB==3, CELLLAB==5);
Q0irrelPVVIP{rix} = Q0irrelPVVIP{rix}(:);
Q0irrelSOMSOM{rix} = Q0irrelV(CELLLAB==4, CELLLAB==4);
Q0irrelSOMSOM{rix}(Q0irrelSOMSOM{rix} == 0) = [];
Q0irrelSOMVIP{rix} = Q0irrelV(CELLLAB==4, CELLLAB==5);
Q0irrelSOMVIP{rix} = Q0irrelSOMVIP{rix}(:);
Q0irrelVIPVIP{rix} = Q0irrelV(CELLLAB==5, CELLLAB==5);
Q0irrelVIPVIP{rix}(Q0irrelVIPVIP{rix} == 0) = [];


end % for NID=1:length(idinf),




clear D
clear D1
clear D2
clear D3
clear D4
clear D5
clear D6
clear S
clear M
clear p
clear p1

%% Mean change in distance for each cell type

D(1,1) = mean(horzcat(QirrelPYRPYR{:}) - horzcat(QrelPYRPYR{:}));
D(1,2) = mean(vertcat(QirrelPYRPV{:}) - vertcat(QrelPYRPV{:}));
D(1,3) = mean(vertcat(QirrelPYRSOM{:}) - vertcat(QrelPYRSOM{:}));
D(1,4) = mean(vertcat(QirrelPYRVIP{:}) - vertcat(QrelPYRVIP{:}));
D(2,2) = mean(horzcat(QirrelPVPV{:})  - horzcat(QrelPVPV{:}));
D(2,3) = mean(vertcat(QirrelPVSOM{:}) - vertcat(QrelPVSOM{:}));
D(2,4) = mean(vertcat(QirrelPVVIP{:}) - vertcat(QrelPVVIP{:}));
D(3,3) = mean(horzcat(QirrelSOMSOM{:}) - horzcat(QrelSOMSOM{:}));
D(3,4) = mean(vertcat(QirrelSOMVIP{:}) - vertcat(QrelSOMVIP{:}));
D(4,4) = mean(horzcat(QirrelVIPVIP{:}) - horzcat(QrelVIPVIP{:}));

%%  Mean Pre and Post distance for each cell type

D1(1,1) = mean(horzcat(QirrelPYRPYR{:}));
D1(1,2) = mean(vertcat(QirrelPYRPV{:}));
D1(1,3) = mean(vertcat(QirrelPYRSOM{:}));
D1(1,4) = mean(vertcat(QirrelPYRVIP{:}));
D1(2,2) = mean(horzcat(QirrelPVPV{:}));
D1(2,3) = mean(vertcat(QirrelPVSOM{:}));
D1(2,4) = mean(vertcat(QirrelPVVIP{:}));
D1(3,3) = mean(horzcat(QirrelSOMSOM{:}));
D1(3,4) = mean(vertcat(QirrelSOMVIP{:}));
D1(4,4) = mean(horzcat(QirrelVIPVIP{:}));

D2(1,1) = mean(horzcat(QrelPYRPYR{:}));
D2(1,2) = mean(vertcat(QrelPYRPV{:}));
D2(1,3) = mean(vertcat(QrelPYRSOM{:}));
D2(1,4) = mean(vertcat(QrelPYRVIP{:}));
D2(2,2) = mean(horzcat(QrelPVPV{:}));
D2(2,3) = mean(vertcat(QrelPVSOM{:}));
D2(2,4) = mean(vertcat(QrelPVVIP{:}));
D2(3,3) = mean(horzcat(QrelSOMSOM{:}));
D2(3,4) = mean(vertcat(QrelSOMVIP{:}));
D2(4,4) = mean(horzcat(QrelVIPVIP{:}));

%% Mean pre and post squared distance for each cell type

D3(1,1) = mean(horzcat(QirrelPYRPYR{:}).^2);
D3(1,2) = mean(vertcat(QirrelPYRPV{:}).^2);
D3(1,3) = mean(vertcat(QirrelPYRSOM{:}).^2);
D3(1,4) = mean(vertcat(QirrelPYRVIP{:}).^2);
D3(2,2) = mean(horzcat(QirrelPVPV{:}).^2);
D3(2,3) = mean(vertcat(QirrelPVSOM{:}).^2);
D3(2,4) = mean(vertcat(QirrelPVVIP{:}).^2);
D3(3,3) = mean(horzcat(QirrelSOMSOM{:}).^2);
D3(3,4) = mean(vertcat(QirrelSOMVIP{:}).^2);
D3(4,4) = mean(horzcat(QirrelVIPVIP{:}).^2);

D4(1,1) = mean(horzcat(QrelPYRPYR{:}).^2);
D4(1,2) = mean(vertcat(QrelPYRPV{:}).^2);
D4(1,3) = mean(vertcat(QrelPYRSOM{:}).^2);
D4(1,4) = mean(vertcat(QrelPYRVIP{:}).^2);
D4(2,2) = mean(horzcat(QrelPVPV{:}).^2);
D4(2,3) = mean(vertcat(QrelPVSOM{:}).^2);
D4(2,4) = mean(vertcat(QrelPVVIP{:}).^2);
D4(3,3) = mean(horzcat(QrelSOMSOM{:}).^2);
D4(3,4) = mean(vertcat(QrelSOMVIP{:}).^2);
D4(4,4) = mean(horzcat(QrelVIPVIP{:}).^2);

%% Mean pre and post distance along a single dimension for each cell type

D5(1,1) = mean(horzcat(Q0irrelPYRPYR{:}));
D5(1,2) = mean(vertcat(Q0irrelPYRPV{:}));
D5(1,3) = mean(vertcat(Q0irrelPYRSOM{:}));
D5(1,4) = mean(vertcat(Q0irrelPYRVIP{:}));
D5(2,2) = mean(horzcat(Q0irrelPVPV{:}));
D5(2,3) = mean(vertcat(Q0irrelPVSOM{:}));
D5(2,4) = mean(vertcat(Q0irrelPVVIP{:}));
D5(3,3) = mean(horzcat(Q0irrelSOMSOM{:}));
D5(3,4) = mean(vertcat(Q0irrelSOMVIP{:}));
D5(4,4) = mean(horzcat(Q0irrelVIPVIP{:}));

D6(1,1) = mean(horzcat(Q0relPYRPYR{:}));
D6(1,2) = mean(vertcat(Q0relPYRPV{:}));
D6(1,3) = mean(vertcat(Q0relPYRSOM{:}));
D6(1,4) = mean(vertcat(Q0relPYRVIP{:}));
D6(2,2) = mean(horzcat(Q0relPVPV{:}));
D6(2,3) = mean(vertcat(Q0relPVSOM{:}));
D6(2,4) = mean(vertcat(Q0relPVVIP{:}));
D6(3,3) = mean(horzcat(Q0relSOMSOM{:}));
D6(3,4) = mean(vertcat(Q0relSOMVIP{:}));
D6(4,4) = mean(horzcat(Q0relVIPVIP{:}));

%% Other quantities of interest
% 
% M(1,1) = mean(horzcat(QirrelPYRPYR{:}) + horzcat(QrelPYRPYR{:}))/2;
% M(1,2) = mean(vertcat(QirrelPYRPV{:}) + vertcat(QrelPYRPV{:}))/2;
% M(1,3) = mean(vertcat(QirrelPYRSOM{:}) + vertcat(QrelPYRSOM{:}))/2;
% M(1,4) = mean(vertcat(QirrelPYRVIP{:}) + vertcat(QrelPYRVIP{:}))/2;
% M(2,2) = mean(horzcat(QirrelPVPV{:})  + horzcat(QrelPVPV{:}))/2;
% M(2,3) = mean(vertcat(QirrelPVSOM{:}) + vertcat(QrelPVSOM{:}))/2;
% M(2,4) = mean(vertcat(QirrelPVVIP{:}) + vertcat(QrelPVVIP{:}))/2;
% M(3,3) = mean(horzcat(QirrelSOMSOM{:}) + horzcat(QrelSOMSOM{:}))/2;
% M(3,4) = mean(vertcat(QirrelSOMVIP{:}) + vertcat(QrelSOMVIP{:}))/2;
% M(4,4) = mean(horzcat(QirrelVIPVIP{:}) + horzcat(QrelVIPVIP{:}))/2;
% 
% S(1,1) = sqrt(var(horzcat(QirrelPYRPYR{:})) + var(horzcat(QrelPYRPYR{:})));
% S(1,2) = sqrt(var(vertcat(QirrelPYRPV{:}))  + var(vertcat(QrelPYRPV{:})));
% S(1,3) = sqrt(var(vertcat(QirrelPYRSOM{:})) + var(vertcat(QrelPYRSOM{:})));
% S(1,4) = sqrt(var(vertcat(QirrelPYRVIP{:})) - var(vertcat(QrelPYRVIP{:})));
% S(2,2) = sqrt(var(horzcat(QirrelPVPV{:}))   - var(horzcat(QrelPVPV{:})));
% S(2,3) = sqrt(var(vertcat(QirrelPVSOM{:}))  - var(vertcat(QrelPVSOM{:})));
% S(2,4) = sqrt(var(vertcat(QirrelPVVIP{:}))  - var(vertcat(QrelPVVIP{:})));
% S(3,3) = sqrt(var(horzcat(QirrelSOMSOM{:})) - var(horzcat(QrelSOMSOM{:})));
% S(3,4) = sqrt(var(vertcat(QirrelSOMVIP{:})) - var(vertcat(QrelSOMVIP{:})));
% S(4,4) = sqrt(var(horzcat(QirrelVIPVIP{:})) - var(horzcat(QrelVIPVIP{:})));
% 
% p(1,1) = signtest(horzcat(QirrelPYRPYR{:}) - horzcat(QrelPYRPYR{:}));
% p(1,2) = signtest(vertcat(QirrelPYRPV{:}) - vertcat(QrelPYRPV{:}));
% p(1,3) = signtest(vertcat(QirrelPYRSOM{:}) - vertcat(QrelPYRSOM{:}));
% p(1,4) = signtest(vertcat(QirrelPYRVIP{:}) - vertcat(QrelPYRVIP{:}));
% p(2,2) = signtest(horzcat(QirrelPVPV{:})  - horzcat(QrelPVPV{:}));
% p(2,3) = signtest(vertcat(QirrelPVSOM{:}) - vertcat(QrelPVSOM{:}));
% p(2,4) = signtest(vertcat(QirrelPVVIP{:}) - vertcat(QrelPVVIP{:}));
% p(3,3) = signtest(horzcat(QirrelSOMSOM{:}) - horzcat(QrelSOMSOM{:}));
% p(3,4) = signtest(vertcat(QirrelSOMVIP{:}) - vertcat(QrelSOMVIP{:}));
% p(4,4) = signtest(horzcat(QirrelVIPVIP{:}) - horzcat(QrelVIPVIP{:}));
% 
% p1(1,1) = ranksum(horzcat(QirrelPYRPYR{:}) , horzcat(QrelPYRPYR{:}));
% p1(1,2) = ranksum(vertcat(QirrelPYRPV{:}) , vertcat(QrelPYRPV{:}));
% p1(1,3) = ranksum(vertcat(QirrelPYRSOM{:}) , vertcat(QrelPYRSOM{:}));
% p1(1,4) = ranksum(vertcat(QirrelPYRVIP{:}) , vertcat(QrelPYRVIP{:}));
% p1(2,2) = ranksum(horzcat(QirrelPVPV{:})  , horzcat(QrelPVPV{:}));
% p1(2,3) = ranksum(vertcat(QirrelPVSOM{:}) , vertcat(QrelPVSOM{:}));
% p1(2,4) = ranksum(vertcat(QirrelPVVIP{:}) , vertcat(QrelPVVIP{:}));
% p1(3,3) = ranksum(horzcat(QirrelSOMSOM{:}) , horzcat(QrelSOMSOM{:}));
% p1(3,4) = ranksum(vertcat(QirrelSOMVIP{:}) , vertcat(QrelSOMVIP{:}));
% p1(4,4) = ranksum(horzcat(QirrelVIPVIP{:}) , horzcat(QrelVIPVIP{:}));
% 
% MTRX(Ni, :, :) = D./M;
% MTRX_S(Ni, :, :) = D./S; 

MTRXirrel(Ni, :, :) = D1;
MTRXrel(Ni, :, :) = D2;

MTRXirrelSQR(Ni, :, :) = D3;
MTRXrelSQR(Ni, :, :) = D4;

MTRXirrel0(Ni, :, :) = D5;
MTRXrel0(Ni, :, :) = D6;

end

 col = hsv(4);
 
 for q=1:4
     
     subplot(2,2,q)
     set(gca, 'fontsize', 18)
     hold on
     box on
     xlabel('Number of Dimensions')
     ylabel('Squared Distance / N')
     
 for i=1:4
     
     if q <= i
         
     plot(squeeze(MTRXPOSTSQR(:,q,i))./([1:N].'), 'color', col(i,:), 'linewidth', 3)
     plot(squeeze(MTRXPRESQR(:,q,i))./([1:N].'),  'color', col(i,:), 'linewidth', 3, 'linestyle', '--')

     
     
     else
         
     plot(squeeze(MTRXPOSTSQR(:,i,q))./([1:N].'), 'color', col(i,:), 'linewidth', 3)
     plot(squeeze(MTRXPRESQR(:,i,q))./([1:N].'),  'color', col(i,:), 'linewidth', 3, 'linestyle', '--')
    
     end
     
 end
 end
 
     legend('PYR ignore', 'PYR attend', 'PV ignore', 'PV attend',  'SOM ignore', 'SOM attend', 'VIP ignore', 'VIP attend')
 
     