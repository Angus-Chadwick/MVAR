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

        %%%% CHECK ABOVE CODE - IT WAS WRONG FOR LEARNING DATA!!!!!
        
if ~exist('DATA1'), DATA1 = ADDTRG; end
T = DATA1{1}.P.t;
X = DATA1{1}.PL{1}.PSTH;

% Need also CELLLAB vector, which isn't yet stored in these data files...

CELLLAB(isnan(CELLLAB)) = [];

CELLLAB0 = CELLLAB;

% Tvec = find(T>0 & T<0.5);
% Xvec = nanmean(X(:,:,Tvec),3);
% XvecPV  = Xvec(:, CELLLAB==3);
% XvecPYR = Xvec(:, CELLLAB==1);
% XvecSOM = Xvec(:, CELLLAB==4);
% XvecVIP = Xvec(:, CELLLAB==5);
% XPV  = mean(XvecPV,2);
% XPYR = mean(XvecPYR,2);
% XSOM = mean(XvecSOM,2);
% XVIP = mean(XvecVIP,2);
% 
% STD = DATA1{1}.PL{1}.PSTHERR .* sqrt(DATA1{1}.PL{1}.PSTHn); % get standard deviations from SEM and counts
% STD = mean(STD(:,:, Tvec), 3);
% STDPYR = STD(:, CELLLAB==1);
% STDPV = STD(:, CELLLAB==3);
% STDSOM = STD(:, CELLLAB == 4);
% STDVIP = STD(:, CELLLAB == 5);
% 
% for i=1:11
% 
%     CORTOTPYR{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORTOT(CELLLAB==1, CELLLAB==1);
%     CORTOTPV{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORTOT(CELLLAB==3, CELLLAB==3);
%     CORTOTSOM{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORTOT(CELLLAB==4, CELLLAB==4);
%     CORTOTVIP{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORTOT(CELLLAB==5, CELLLAB==5);
% 
%     CORTOTPYRSOM{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORTOT(CELLLAB==1, CELLLAB==4);
%     MEANCORRPYRSOM{i} = nanmean(CORTOTPYRSOM{i}(:));
%     
%     MEANCORRPYR(i) = nanmean(CORTOTPYR{i}(:));
%     MEANCORRPV(i) = nanmean(CORTOTPV{i}(:));
%     MEANCORRSOM(i) = nanmean(CORTOTSOM{i}(:));
%     MEANCORRVIP(i) = nanmean(CORTOTVIP{i}(:));
%     
%     CORNSEPYR{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORNSE1(CELLLAB==1, CELLLAB==1);
%     CORNSEPV{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORNSE1(CELLLAB==3, CELLLAB==3);
%     CORNSESOM{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORNSE1(CELLLAB==4, CELLLAB==4);
%     CORNSEVIP{i} = ADDTRG{1}.PL{1}.DCOR{i}.CORNSE1(CELLLAB==5, CELLLAB==5);
% 
%     MEANCORRNSEPYR(i) = nanmean(CORNSEPYR{i}(:));
%     MEANCORRNSEPV(i) = nanmean(CORNSEPV{i}(:));
%     MEANCORRNSESOM(i) = nanmean(CORNSESOM{i}(:));
%     MEANCORRNSEVIP(i) = nanmean(CORNSEVIP{i}(:));
%     
% end

% [YPYR, IPYR] = sort(XPYR); 
% 
% [CELLS, CELLSIND] = sort(CELLLAB);
% 
% MAT = ADDTRG{1}.PL{1}.DCOR{2}.CORNSE1;  % noise correlations for vertical relevant condition
% MAT(isnan(MAT)) = 0;
% MAT = (MAT + MAT.')/2;
% 
% I = nanstd(MAT);
% 
% [I, Isort] = sort(I);

m_vert_rel = ADDTRG{1}.PL{1}.TwinPP.m1(1,:,1);
m_ang_rel  = ADDTRG{1}.PL{1}.TwinPP.m1(1,:,2);

m_vert_irrel = ADDTRG{1}.PL{1}.TwinPP.m1(1,:,5);
m_ang_irrel = ADDTRG{1}.PL{1}.TwinPP.m1(1,:,6);

v_vert_rel   = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,1);
v_ang_rel = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,2);

v_vert_irrel   = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,5);
v_ang_irrel = ADDTRG{1}.PL{1}.TwinPP.v1(1,:,6);

n_vert_rel = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,1);
n_ang_rel  = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,2);

n_vert_irrel = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,5);
n_ang_irrel  = ADDTRG{1}.PL{1}.TwinPP.n1(1,:,6);

poolstd_rel   = sqrt( (v_vert_rel   .* (n_vert_rel - 1)   + v_ang_rel   .* (n_ang_rel - 1))   ./ (n_vert_rel   + n_ang_rel   - 2) );
poolstd_irrel = sqrt( (v_vert_irrel .* (n_vert_irrel - 1) + v_ang_irrel .* (n_ang_irrel - 1)) ./ (n_vert_irrel + n_ang_irrel - 2) );

Srel_all   = (m_vert_rel   - m_ang_rel)   ./ poolstd_rel;
Sirrel_all = (m_vert_irrel - m_ang_irrel) ./ poolstd_irrel;

%%% IMPORTANT - REJECT OUTLIERS, THIS STRONGLY AFFECTS GAUSSIAN MODEL
%%% FITTING.

Nstd = 5; % number of standard deviations from mean to tolerate 

SI = [Sirrel_all; Srel_all].'; 
SIPYR = SI(CELLLAB == 1, :);

SIG = cov(SIPYR);  
MU  = mean(SIPYR);

MD = sqrt(diag((SI.' - repmat(MU.', 1, size(SI, 1))).' * inv(SIG) * (SI.' - repmat(MU.', 1, size(SI, 1))))); % Malahabinois distance of each cell

reject = and(MD > Nstd, CELLLAB==1); % reject only pyramidal cells whose selectivity is outwith tolerance bounds
accept = ~reject;

% Nmad = 5; % number of median absolute deviations to tolerate
% 
% reject = or(abs(Sirrel_all - median(Sirrel_all)) > Nmad * mad(Sirrel_all, 1), abs(Srel_all - median(Srel_all)) > Nmad * mad(Srel_all, 1));
% accept = ~reject;
% 

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

%% Investigate spectrum of noise correlation matrix

% Plot eigenspectra

hold on

COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

lambda = flipud(diag(D));

plot(lambda, 'linewidth', 3, 'color', 'r')


COND = 2; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

lambda = flipud(diag(D));

plot(lambda, 'linewidth', 3, 'color', 'b')

COND = 5; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

lambda = flipud(diag(D));

plot(lambda, 'linewidth', 3, 'linestyle', '--', 'color', 'r')

COND = 6; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

lambda = flipud(diag(D));

plot(lambda, 'linewidth', 3, 'linestyle', '--', 'color', 'b')

set(gca, 'fontsize', 18)
legend('Relevant, Vertical', 'Relevant, Angled', 'Irrelevant, Vertical', 'Irrelevant, Angled')
xlabel('Eigenvalue #')
ylabel('Eigenvalue')
box on
axis([0 5 0 25])


% Show eigenvector entries for first 4 modes in grey corridor condition

COND = 11; % experimental condition
MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor onset condition

%COND = 1;
%MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor preceding vertical relevant condition

MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(411)

if median(V(:,end)) > 0

    hist(V(:,end), 50)

elseif median(V(:,end)) < 0
    
    hist(-V(:,end), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('First Mode')
box on
axis([-0.25 0.25 0 30])

subplot(412)

if median(V(:,end-1)) > 0

    hist(V(:,end-1), 50)

elseif median(V(:,end-1)) < 0
    
    hist(-V(:,end-1), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Second Mode')
box on
axis([-0.25 0.25 0 30])

subplot(413)

if median(V(:,end-2)) > 0

    hist(V(:,end-2), 50)

elseif median(V(:,end-2)) < 0
    
    hist(-V(:,end-2), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Third Mode')
box on
axis([-0.25 0.25 0 30])

subplot(414)

if median(V(:,end-3)) > 0

    hist(V(:,end-3), 50)

elseif median(V(:,end-3)) < 0
    
    hist(-V(:,end-3), 50)  % convention is to take positive eigenvectors
    
end
set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Fourth Mode')
box on
axis([-0.25 0.25 0 30])


% Show nth mode in different stimulus conditions

n = 1; % mode number

COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(411)

if median(V(:,end - (n-1))) > 0

    hist(V(:,end - (n-1)), 50)

elseif median(V(:,end - (n-1))) < 0
    
    hist(-V(:,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Vertical, Relevant')
box on
axis([-0.25 0.25 0 30])

COND = 2; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(412)

if median(V(:,end - (n-1))) > 0

    hist(V(:,end - (n-1)), 50)

elseif median(V(:,end - (n-1))) < 0
    
    hist(-V(:,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Angled, Relevant')
box on
axis([-0.25 0.25 0 30])

COND = 5; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(413)

if median(V(:,end - (n-1))) > 0

    hist(V(:,end - (n-1)), 50)

elseif median(V(:,end - (n-1))) < 0
    
    hist(-V(:,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Vertical, Irrelevant')
box on
axis([-0.25 0.25 0 30])

COND = 6; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(414)

if median(V(:,end - (n-1))) > 0

    hist(V(:,end - (n-1)), 50)

elseif median(V(:,end - (n-1))) < 0
    
    hist(-V(:,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('Angled, Irrelevant')
box on
axis([-0.25 0.25 0 30])

% Investigate  how different cell types engage in each mode in the grey corridor

n=1; % mode

COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

subplot(411)

if median(V(CELLLAB0==1,end - (n-1))) > 0

    hist(V(CELLLAB0==1,end - (n-1)), 50)

elseif median(V(CELLLAB0==1,end - (n-1))) < 0
    
    hist(-V(CELLLAB0==1,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end


set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('PYR')
box on
axis([-0.25 0.25 0 30])

subplot(412)

if median(V(CELLLAB0==1,end - (n-1))) > 0

    hist(V(CELLLAB0==3,end - (n-1)), 50)

elseif median(V(CELLLAB0==1,end - (n-1))) < 0
    
    hist(-V(CELLLAB0==3,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('PV')
box on
axis([-0.25 0.25 0 3])

subplot(413)

if median(V(CELLLAB0==1,end - (n-1))) > 0

    hist(V(CELLLAB0==4,end - (n-1)), 50)

elseif median(V(CELLLAB0==1,end - (n-1))) < 0
    
    hist(-V(CELLLAB0==4,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('SOM')
box on
axis([-0.25 0.25 0 3])

subplot(414)

if median(V(CELLLAB0==1,end - (n-1))) > 0

    hist(V(CELLLAB0==5,end - (n-1)), 50)

elseif median(V(CELLLAB0==1,end - (n-1))) < 0
    
    hist(-V(CELLLAB0==5,end - (n-1)), 50)  % convention is to take positive eigenvectors
    
end

set(gca, 'fontsize', 18)
xlabel('Eigenvector entry')
ylabel('Number of cells')
title('VIP')
box on
axis([-0.25 0.25 0 3])

% Scatter plots of different cells in different eigenmodes (first 4)

COND = 5; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
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


subplot(221)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 1), 100, '.', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 1), 100, '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 1), 100, 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 1), 100, 'k','+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 2')
box on

subplot(222)

hold on

scatter(V(CELLLAB0 == 1, end-2), V(CELLLAB0 == 1, end - 1), 100, '.', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end-2), V(CELLLAB0 == 3, end - 1), 100, '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end-2), V(CELLLAB0 == 4, end - 1), 100, 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end-2), V(CELLLAB0 == 5, end - 1), 100, 'k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 2')
box on


subplot(223)

hold on

scatter(V(CELLLAB0 == 1, end), V(CELLLAB0 == 1, end - 3), 100, '.', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end), V(CELLLAB0 == 3, end - 3), 100, '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end), V(CELLLAB0 == 4, end - 3), 100, 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end), V(CELLLAB0 == 5, end - 3), 100, 'k', '+', 'linewidth', 2)

set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 1')
ylabel('Mode 4')
box on


subplot(224)

hold on
 
scatter(V(CELLLAB0 == 1, end - 2), V(CELLLAB0 == 1, end-3), 100, '.', 'linewidth', 2)
scatter(V(CELLLAB0 == 3, end - 2), V(CELLLAB0 == 3, end-3), 100, '*', 'linewidth', 2)
scatter(V(CELLLAB0 == 4, end - 2), V(CELLLAB0 == 4, end-3), 100, 'o', 'linewidth', 2)
scatter(V(CELLLAB0 == 5, end - 2), V(CELLLAB0 == 5, end-3), 100, 'k', '+', 'linewidth', 2)
 
set(gca, 'fontsize', 18)
legend('PYR', 'PV', 'SOM', 'VIP')
xlabel('Mode 3')
ylabel('Mode 4')
box on

%% Look for gating by SOMs - subtract noise correlations in the angled vs vertical gratings and compare results for relevant vs irrelevant conditions

%%%%% ASK JASPER WHY THERE IS A STRONG DIFFERENCE BETWEEN VERTICAL AND
%%%%% ANGLED IN THE GREY PRECEDING CONDITION!!!!!

COND = 1;

MATVrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MATVrel(isnan(MATVrel)) = 0;
MATVrel(logical(eye(size(MATVrel)))) = 1;
MATVrel = MATVrel + triu(MATVrel,1).';


COND = 2;

MATArel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MATArel(isnan(MATArel)) = 0;
MATArel(logical(eye(size(MATArel)))) = 1;
MATArel = MATArel + triu(MATArel,1).';

Mrel = MATVrel - MATArel;

COND = 5;

MATVirrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MATVirrel(isnan(MATVirrel)) = 0;
MATVirrel(logical(eye(size(MATVirrel)))) = 1;
MATVirrel = MATVirrel + triu(MATVirrel,1).';

COND = 6;

MATAirrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MATAirrel(isnan(MATAirrel)) = 0;
MATAirrel(logical(eye(size(MATAirrel)))) = 1;
MATAirrel = MATAirrel + triu(MATAirrel,1).';

Mirrel = MATVirrel - MATAirrel;

MPVrel = Mrel(CELLLAB == 1, CELLLAB == 3);
MSOMrel = Mrel(CELLLAB == 1, CELLLAB == 4);
MVIPrel = Mrel(CELLLAB == 1, CELLLAB == 5);
MPVirrel = Mirrel(CELLLAB == 1, CELLLAB == 3);
MSOMirrel = Mirrel(CELLLAB == 1, CELLLAB == 4);
MVIPirrel = Mirrel(CELLLAB == 1, CELLLAB == 5);

MPVrel(isnan(MPVrel)) = [];
MSOMrel(isnan(MSOMrel)) = [];
MVIPrel(isnan(MVIPrel)) = [];
MPVirrel(isnan(MPVirrel)) = [];
MSOMirrel(isnan(MSOMirrel)) = [];
MVIPirrel(isnan(MVIPirrel)) = [];



%% Calculate distance between cells in eigenspace

COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

N = 4; % dimension of eigenspace

% PYR-PYR

QrelV = zeros(length(CELLLAB));

for i=1:length(CELLLAB)
    for j=1:length(CELLLAB)
        
QrelV(i,j) = mean(norm(V(i, (end - (N-1)):end) - V(j, (end - (N-1)):end)), 2);

    end
end


COND = 5; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE4;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

[V, D] = eig(MAT); % eigenspectrum of noise correlation matrix

N = 4; % dimension of eigenspace

% PYR-PYR

QirrelV = zeros(length(CELLLAB));

for i=1:length(CELLLAB)
    for j=1:length(CELLLAB)
        
QirrelV(i,j) = mean(norm(V(i, (end - (N-1)):end) - V(j, (end - (N-1)):end)), 2);

    end
end

QrelPYRPYR = QrelV(CELLLAB==1, CELLLAB==1);
QrelPYRPYR(QrelPYRPYR == 0) = [];
QrelPYRPV = QrelV(CELLLAB==1, CELLLAB==3);
QrelPYRPV = QrelPYRPV(:);
QrelPYRSOM = QrelV(CELLLAB==1, CELLLAB==4);
QrelPYRSOM = QrelPYRSOM(:);
QrelPYRVIP = QrelV(CELLLAB==1, CELLLAB==5);
QrelPYRVIP = QrelPYRVIP(:);
QrelPVPV = QrelV(CELLLAB==3, CELLLAB==3);
QrelPVPV(QrelPVPV == 0) = [];
QrelPVSOM = QrelV(CELLLAB==3, CELLLAB==4);
QrelPVSOM = QrelPVSOM(:);
QrelPVVIP = QrelV(CELLLAB==3, CELLLAB==5);
QrelPVVIP = QrelPVVIP(:);
QrelSOMSOM = QrelV(CELLLAB==4, CELLLAB==4);
QrelSOMSOM(QrelSOMSOM == 0) = [];
QrelSOMVIP = QrelV(CELLLAB==4, CELLLAB==5);
QrelSOMVIP = QrelSOMVIP(:);
QrelVIPVIP = QrelV(CELLLAB==5, CELLLAB==5);
QrelVIPVIP(QrelVIPVIP == 0) = [];

QirrelPYRPYR = QirrelV(CELLLAB==1, CELLLAB==1);
QirrelPYRPYR(QirrelPYRPYR == 0) = [];
QirrelPYRPV = QirrelV(CELLLAB==1, CELLLAB==3);
QirrelPYRPV = QirrelPYRPV(:);
QirrelPYRSOM = QirrelV(CELLLAB==1, CELLLAB==4);
QirrelPYRSOM = QirrelPYRSOM(:);
QirrelPYRVIP = QirrelV(CELLLAB==1, CELLLAB==5);
QirrelPYRVIP = QirrelPYRVIP(:);
QirrelPVPV = QirrelV(CELLLAB==3, CELLLAB==3);
QirrelPVPV(QirrelPVPV == 0) = [];
QirrelPVSOM = QirrelV(CELLLAB==3, CELLLAB==4);
QirrelPVSOM = QirrelPVSOM(:);
QirrelPVVIP = QirrelV(CELLLAB==3, CELLLAB==5);
QirrelPVVIP = QirrelPVVIP(:);
QirrelSOMSOM = QirrelV(CELLLAB==4, CELLLAB==4);
QirrelSOMSOM(QirrelSOMSOM == 0) = [];
QirrelSOMVIP = QirrelV(CELLLAB==4, CELLLAB==5);
QirrelSOMVIP = QirrelSOMVIP(:);
QirrelVIPVIP = QirrelV(CELLLAB==5, CELLLAB==5);
QirrelVIPVIP(QirrelVIPVIP == 0) = [];

clear FracChange

% PYR distances

FracChange(1,1) = median(QrelPYRPYR)/median(QirrelPYRPYR);
FracChange(1,2) = median(QrelPYRPV)/median(QirrelPYRPV);
FracChange(1,3) = median(QrelPYRSOM)/median(QirrelPYRSOM);
FracChange(1,4) = median(QrelPYRVIP)/median(QirrelPYRVIP);

% PV distances;

FracChange(2,2) = median(QrelPVPV)/median(QirrelPVPV);
FracChange(2,3) = median(QrelPVSOM)/median(QirrelPVSOM);
FracChange(2,4) = median(QrelPVVIP)/median(QirrelPVVIP);

% SOM distances

FracChange(3,3) = median(QrelSOMSOM)/median(QirrelSOMSOM);
FracChange(3,4) = median(QrelSOMVIP)/median(QirrelSOMVIP);

% VIP distances

FracChange(4,4) = median(QrelVIPVIP)/median(QirrelVIPVIP);

FracChange = FracChange + triu(FracChange,1).';
%% Plot changes in selectivity for cells with different preferences in different conditions

X = linspace(min(dS), max(dS), 100);

subplot(221)

hist(dS(Sirrel < 0), X)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with negative selectivity in irrelevant condition')

subplot(223)

hist(dS(Sirrel > 0), X)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with positive selectivity in irrelevant condition')

subplot(222)

hist(dS(Srel < 0), X)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with negative selectivity in relevant condition')

subplot(224)

hist(dS(Srel > 0), X)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with positive selectivity in relevant condition')

%% Plot changes in response to different stimuli for cells with different preferences

subplot(221)

hist(dFang(Sirrel < 0), 25)

axis([-1 1 0 35])
xlabel('Change in response')
ylabel('Number of cells')
title('Angled stimulus, cells preferring angled stimulus')

subplot(223)

hist(dFang(Sirrel > 0), 25)

axis([-1 1 0 35])
xlabel('Change in response')
ylabel('Number of cells')
title('Angled stimulus, cells preferring vertical stimulus')

subplot(222)

hist(dFvert(Sirrel < 0), 25)

axis([-1 1 0 35])
xlabel('Change in response')
ylabel('Number of cells')
title('Vertical stimulus, cells preferring angled stimulus')

subplot(224)

hist(dFang(Sirrel > 0), 25)

axis([-1 1 0 35])
xlabel('Change in response')
ylabel('Number of cells')
title('Vertical stimulus, cells preferring vertical stimulus')


%% Look at a linear regression model to predict selectivity from noise correlations

% 
% XX = [XPYR, XPV, XSOM, XVIP];
% 
% betatot = (XX.' * XX)^(-1) * XX.' * (Srel - mean(Srel)).';
% 
% ModPred = XX * betatot;
% resid = (Srel - mean(Srel)) - ModPred.';
% 
% SStot = sum((Srel - mean(Srel)).^2);
% SSres = sum(resid.^2);
% 
% R2 = 1 - SSres/SStot;
% 
% 
% MAT = ADDTRG{1}.PL{1}.DCOR{11}.CORNSE1;  % noise correlations for grey corridor condition
% MAT(isnan(MAT)) = 0;
% MAT = (MAT + MAT.')/2;
% 
% M = MAT(CELLLAB==1, CELLLAB==5);
% XVIP = M * SrelVIP.';
% XVIP = XVIP - mean(XVIP);
% betaVIP = (XVIP.' * XVIP)^(-1) * XVIP.' * (Srel - mean(Srel)).';
% 
% [cVIP, p] = corr(Srel.', XVIP);
% 
% M = MAT(CELLLAB==1, CELLLAB==4);
% XSOM = M * SrelSOM.';
% XSOM = XSOM - mean(XSOM);
% betaSOM = (XSOM.' * XSOM)^(-1) * XSOM.' * (Srel - mean(Srel)).';
% 
% [cSOM, p] = corr(Srel.', XSOM);
% 
% M = MAT(CELLLAB==1, CELLLAB==3);
% XPV = M * SrelPV.';
% XPV = XPV - mean(XPV);
% betaPV = (XPV.' * XPV)^(-1) * XPV.' * (Srel - mean(Srel)).';
% 
% [cPV, p] = corr(Srel.', XPV);
% 
% M = MAT(CELLLAB==1, CELLLAB==1);
% XPYR = M * Srel.';
% XPYR = XPYR - mean(XPYR);
% betaPYR = (XPYR.' * XPYR)^(-1) * XPYR.' * (Srel - mean(Srel)).';
% 
% [cPYR, p] = corr(Srel.', XPYR);
% 
% beta = [betaPYR, betaPV, betaSOM, betaVIP];
% 
% XX = [XPYR, XPV, XSOM, XVIP, Sirrel.' - mean(Sirrel)];
% 
% betatot = (XX.' * XX)^(-1) * XX.' * (Srel - mean(Srel)).';
% 
% ModPred = XX * betatot;
% resid = (Srel - mean(Srel)) - ModPred.';
% 
% SStot = sum((Srel - mean(Srel)).^2);
% SSres = sum(resid.^2);
% 
% R2 = 1 - SSres/SStot;

%% Look at a linear regression model to predict selectivity CHANGES from noise correlations and peer selectivity CHANGES

% MAT = ADDTRG{1}.PL{1}.DCOR{11}.CORNSE1;  % noise correlations for grey corridor
% MAT(isnan(MAT)) = 0;
% MAT = (MAT + MAT.')/2;
% 
% dSSOM = (SrelSOM - SirrelSOM);
% dSPV  = (SrelPV  - SirrelPV);
% dSVIP = (SrelVIP - SirrelVIP);
% 
% M = MAT(CELLLAB==1, CELLLAB==5);
% XVIP = M * dSVIP.';
% XVIP = XVIP - mean(XVIP);
% betaVIP = (XVIP.' * XVIP)^(-1) * XVIP.' * (dS - mean(dS)).';
% 
% [cVIP, p] = corr(dS.', XVIP);
% 
% M = MAT(CELLLAB==1, CELLLAB==4);
% XSOM = M * dSSOM.';
% XSOM = XSOM - mean(XSOM);
% betaSOM = (XSOM.' * XSOM)^(-1) * XSOM.' * (dS - mean(dS)).';
% 
% [cSOM, p] = corr(dS.', XSOM);
% 
% M = MAT(CELLLAB==1, CELLLAB==3);
% XPV = M * dSPV.';
% XPV = XPV - mean(XPV);
% betaPV = (XPV.' * XPV)^(-1) * XPV.' * (dS - mean(dS)).';
% 
% [cPV, p] = corr(dS.', XPV);
% 
% M = MAT(CELLLAB==1, CELLLAB==1);
% XPYR = M * dS.';
% XPYR = XPYR - mean(XPYR);
% betaPYR = (XPYR.' * XPYR)^(-1) * XPYR.' * (dS - mean(dS)).';
% 
% [cPYR, p] = corr(dS.', XPYR);
% 
% beta = [betaPYR, betaPV, betaSOM, betaVIP];
% 
% XX = [XPYR, XPV, XSOM, XVIP, Sirrel.' - mean(Sirrel)];
% 
% betatot = (XX.' * XX)^(-1) * XX.' * (dS - mean(dS)).';
% 
% ModPred = XX * betatot;
% resid = (dS - mean(dS)) - ModPred.';
% 
% SStot = sum((dS - mean(dS)).^2);
% SSres = sum(resid.^2);
% 
% R2 = 1 - SSres/SStot;


%% Look at joint Gaussian model with mean and covariance specified by the data

% Use grey onset condition, without separating by task-relevance

% COND = 11;
% 
% %dF = (m_vert_rel + m_ang_rel - (m_vert_irrel + m_ang_irrel)) ./ (m_vert_rel + m_ang_rel + m_vert_irrel + m_ang_irrel);
% dF = (m_vert_rel + m_ang_rel - (m_vert_irrel + m_ang_irrel));
% dF = dF(accept);
% 
% MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
% MAT(isnan(MAT)) = 0;
% MAT = (MAT + MAT.')/2;  %% symmetrise matrix
% MAT = MAT(accept, accept);
% 
% M = MAT(CELLLAB==1, CELLLAB==1);
% XPYR = M * Srel.' / length(Srel);
% dFPYR = dF(CELLLAB == 1); % Fractional change in PYR firing between conditions
% XPYR_RATE = M * dFPYR.'; % firing rate of PYR multiplied by noise correlations
% 
% M = MAT(CELLLAB==1, CELLLAB==3);
% XPV = M * SrelPV.' / length(SrelPV);
% dFPV = dF(CELLLAB == 3); % Fractional change in PV firing between conditions
% XPV_RATE = M * dFPV.'; % firing rate of PV multiplied by noise correlations
% 
% M = MAT(CELLLAB==1, CELLLAB==4);
% XSOM = M * SrelSOM.' / length(SrelSOM);
% dFSOM = dF(CELLLAB == 4); % Fractional change in SOM firing between conditions
% XSOM_RATE = M * dFSOM.'; % firing rate of SOM multiplied by noise correlations
% 
% M = MAT(CELLLAB==1, CELLLAB==5);
% XVIP = M * SrelVIP.' / length(SrelVIP);
% dFVIP = dF(CELLLAB == 5); % Fractional change in VIP firing between conditions
% XVIP_RATE = M * dFVIP.'; % firing rate of VIP multiplied by noise correlations

% Use grey preceding stimulus, separated by task-relevance

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

COND = 11;

MATirrel = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MATirrel(isnan(MATirrel)) = 0;
MATirrel = MATirrel + triu(MATirrel,1).';  %% symmetrise matrix
MATirrel = MATirrel(accept, accept);

%MATrel = (MATrel + MATirrel)/2;
%MATirrel = (MATrel + MATirrel)/2;

%% Do cross validation using the above method

clear XPYR_TRAIN
clear XPYR_RATE_TRAIN
clear XPYR_dS_TRAIN

clear XPV_TRAIN
clear XPV_RATE_TRAIN
clear XPV_dS_TRAIN

clear XSOM_TRAIN
clear XSOM_RATE_TRAIN
clear XSOM_dS_TRAIN

clear XVIP_TRAIN
clear XVIP_RATE_TRAIN
clear XVIP_dS_TRAIN

clear XPYR_TEST
clear XPYR_RATE_TEST
clear XPYR_dS_TEST

clear XPV_TEST
clear XPV_RATE_TEST
clear XPV_dS_TEST

clear XSOM_TEST
clear XSOM_RATE_TEST
clear XSOM_dS_TEST

clear XVIP_TEST
clear XVIP_RATE_TEST
clear XVIP_dS_TEST

clear Srel_TRAIN
clear Sirrel_TRAIN

NPYR = sum(CELLLAB == 1);

% PRM = randperm(size(MATrel,1));
%  PRM2 = randperm(size(MATrel,1));
%  MATrel = MATrel(PRM2, PRM2);
%  MATirrel = MATirrel(PRM2, PRM2);

%PRM = randperm(length(Srel));
%Srel = Srel(PRM);
%Sirrel = Sirrel(PRM);
% 
% for i=1:size(MATrel,1)
%     
%    M = MATrel([1:(i-1), (i+1):size(MATrel,1)],i);
%    MATrel([1:(i-1), (i+1):size(MATrel,1)],i) = M(randperm(length(M)));
%    
% end
% 
% MATirrel = MATrel;

for i=1:NPYR
    
    % training datasets - remove one pyramidal cell and compute prediction
    % metrics
    
    ind0 = find(CELLLAB == 1);
    ind1 = [1:(i-1), (i+1):NPYR];
    ind  = ind0(ind1); 
   
    %M = MAT(ind, ind);
    %XPYR_TRAIN(i,:) = M * Srel(ind1).' / length(ind1);
    %XPYR_RATE_TRAIN(i,:) = M * dFPYR(ind1).' / length(ind1); % firing rate of PYR multiplied by noise correlations
    %XPYR_dS_TRAIN(i,:) = M * dS_PYR(ind1).' / length(ind1);

    Mrel = MATrel(ind, ind);
    Mirrel = MATirrel(ind, ind);
    XPYR_TRAIN(i,:) = Mrel * Srel(ind1).' / length(ind1);
    XPYR_RATE_TRAIN(i,:) = Mrel * dFPYR(ind1).' / length(ind1); % firing rate of PYR multiplied by noise correlations
    XPYR_dS_TRAIN(i,:) = (Mrel * Srel(ind1).' - Mirrel * Sirrel(ind1).') / length(ind1);
    
%     M = MAT(ind, CELLLAB==3);
%     XPV_TRAIN(i,:) = M * SrelPV.' / length(SrelPV);
%     XPV_RATE_TRAIN(i,:) = M * dFPV.' / length(dFPV); % firing rate of PV multiplied by noise correlations
%     XPV_dS_TRAIN(i,:) = M * dS_PV.' / length(dS_PV);

    Mrel = MATrel(ind, CELLLAB==3);
    Mirrel = MATirrel(ind, CELLLAB==3);
    XPV_TRAIN(i,:) = Mrel * SrelPV.' / length(SrelPV);
    XPV_RATE_TRAIN(i,:) = Mrel * dFPV.' / length(dFPV); % firing rate of PV multiplied by noise correlations
    XPV_dS_TRAIN(i,:) = (Mrel * SrelPV.' - Mirrel * SirrelPV.') / length(dS_PV);

    
%     M = MAT(ind, CELLLAB==4);
%     XSOM_TRAIN(i,:) = M * SrelSOM.' / length(SrelSOM);
%     XSOM_RATE_TRAIN(i,:) = M * dFSOM.' / length(dFSOM); % firing rate of PV multiplied by noise correlations
%     XSOM_dS_TRAIN(i,:) = M * dS_SOM.' / length(dS_SOM);

    Mrel = MATrel(ind, CELLLAB==4);
    Mirrel = MATirrel(ind, CELLLAB==4);
    XSOM_TRAIN(i,:) = Mrel * SrelSOM.' / length(SrelSOM);
    XSOM_RATE_TRAIN(i,:) = Mrel * dFSOM.' / length(dFSOM); % firing rate of PV multiplied by noise correlations
    XSOM_dS_TRAIN(i,:) = (Mrel * SrelSOM.' - Mirrel * SirrelSOM.') / length(dS_SOM);
    
    
%     M = MAT(ind, CELLLAB==5);
%     XVIP_TRAIN(i,:) = M * SrelVIP.' / length(SrelVIP);
%     XVIP_RATE_TRAIN(i,:) = M * dFVIP.' / length(dFVIP); % firing rate of PV multiplied by noise correlations
%     XVIP_dS_TRAIN(i,:) = M * dS_VIP.' / length(dS_VIP);
       
    Mrel = MATrel(ind, CELLLAB==5);
    Mirrel = MATirrel(ind, CELLLAB==5);
    XVIP_TRAIN(i,:) = Mrel * SrelVIP.' / length(SrelVIP);
    XVIP_RATE_TRAIN(i,:) = Mrel * dFVIP.' / length(dFVIP); % firing rate of PV multiplied by noise correlations
    XVIP_dS_TRAIN(i,:) = (Mrel * SrelVIP.' - Mirrel * SirrelVIP.') / length(dS_VIP);
    
    Srel_TRAIN(i,:) = Srel(ind1);
    Sirrel_TRAIN(i,:) = Sirrel(ind1);
    
    % test datasets - use all peer cells for prediction
    
%     M = MAT(CELLLAB==1, CELLLAB==1);
%     XPYR_TEST(i) = M(i,:) * Srel.' / length(Srel);
%     XPYR_RATE_TEST(i) = M(i,:) * dFPYR.' / length(dFPYR);
%     XPYR_dS_TEST(i) = M(i,:) * dS_PYR.' / length(dS_PYR);
    
    Mrel = MATrel(CELLLAB==1, CELLLAB==1);
    Mirrel = MATirrel(CELLLAB==1, CELLLAB==1);     
    XPYR_TEST(i) = Mrel(i,:) * Srel.' / length(Srel);
    XPYR_RATE_TEST(i) = Mrel(i,:) * dFPYR.' / length(dFPYR);
    XPYR_dS_TEST(i) = (Mrel(i,:) * Srel.' - Mirrel(i,:) * Sirrel.') / length(dS_PYR);
    
    
%     M = MAT(CELLLAB==1, CELLLAB==3);
%     XPV_TEST(i) = M(i,:) * SrelPV.' / length(SrelPV);
%     XPV_RATE_TEST(i) = M(i,:) * dFPV.' / length(dFPV);
%     XPV_dS_TEST(i) = M(i,:) * dS_PV.' / length(dS_PV);
        
    Mrel = MATrel(CELLLAB==1, CELLLAB==3);
    Mirrel = MATirrel(CELLLAB==1, CELLLAB==3);      
    XPV_TEST(i) = Mrel(i,:) * SrelPV.' / length(SrelPV);
    XPV_RATE_TEST(i) = Mrel(i,:) * dFPV.' / length(dFPV);
    XPV_dS_TEST(i) = (Mrel(i,:) * SrelPV.' - Mirrel(i,:) * SirrelPV.') / length(dS_PV);

%     M = MAT(CELLLAB==1, CELLLAB==4);
%     XSOM_TEST(i) = M(i,:) * SrelSOM.' / length(SrelSOM);
%     XSOM_RATE_TEST(i) = M(i,:) * dFSOM.' / length(dFSOM);
%     XSOM_dS_TEST(i) = M(i,:) * dS_SOM.' / length(dS_SOM);
    
    Mrel = MATrel(CELLLAB==1, CELLLAB==4);
    Mirrel = MATirrel(CELLLAB==1, CELLLAB==4);      
    XSOM_TEST(i) = Mrel(i,:) * SrelSOM.' / length(SrelSOM);
    XSOM_RATE_TEST(i) = Mrel(i,:) * dFSOM.' / length(dFSOM);
    XSOM_dS_TEST(i) = (Mrel(i,:) * SrelSOM.' - Mirrel(i,:) * SirrelSOM.') / length(dS_SOM);
    
    
%     M = MAT(CELLLAB==1, CELLLAB==5);
%     XVIP_TEST(i) = M(i,:) * SrelVIP.' / length(SrelVIP);
%     XVIP_RATE_TEST(i) = M(i,:) * dFVIP.' / length(dFVIP);
%     XVIP_dS_TEST(i) = M(i,:) * dS_VIP.' / length(dS_VIP);
           
    Mrel = MATrel(CELLLAB==1, CELLLAB==5);
    Mirrel = MATirrel(CELLLAB==1, CELLLAB==5);    
    XVIP_TEST(i) = Mrel(i,:) * SrelVIP.' / length(SrelVIP);
    XVIP_RATE_TEST(i) = Mrel(i,:) * dFVIP.' / length(dFVIP);
    XVIP_dS_TEST(i) = (Mrel(i,:) * SrelVIP.' - Mirrel(i,:) * SirrelVIP.') / length(dS_VIP);

    
end




METHOD = 2;  % 1 for 10 fold cross validation and 2 for leave-one-out cross validation

if METHOD == 1 % if N fold cross validation

elseif METHOD == 2 % leave-one-out cross validation
    
     NPYR = sum(CELLLAB == 1);
 
    clear condmean_rel_peers
    clear condmean_rel_null
    clear condmean_rel_COMB
    clear condmean_rel_RATE
    clear condcov_rel_peers
    clear condcov_rel_null
    clear condcov_rel_COMB
    clear condcov_rel_RATE
    clear LNULL
    clear LTEST 
    clear LRATE
    clear LCOMB
    
    for i=1:NPYR
            
       DATA_TEST  = [Srel(i); XPYR_dS_TEST(i); XPV_dS_TEST(i); XSOM_dS_TEST(i); XVIP_dS_TEST(i); Sirrel(i)].';
       DATA_TRAIN = [Srel_TRAIN(i,:); XPYR_dS_TRAIN(i,:); XPV_dS_TRAIN(i,:); XSOM_dS_TRAIN(i,:); XVIP_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
   
       DATA_TEST_RATE  = [Srel(i); XPYR_RATE_TEST(i); XPV_RATE_TEST(i); XSOM_RATE_TEST(i); XVIP_RATE_TEST(i); Sirrel(i)].';
       DATA_TRAIN_RATE = [Srel_TRAIN(i,:); XPYR_RATE_TRAIN(i,:); XPV_RATE_TRAIN(i,:); XSOM_RATE_TRAIN(i,:); XVIP_RATE_TRAIN(i,:); Sirrel_TRAIN(i,:)].';

      %DATA_TEST_COMB  = [Srel(i); XSOM_dS_TEST(i); Sirrel(i)].';
      %DATA_TRAIN_COMB = [Srel_TRAIN(i,:); XSOM_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
        
       DATA_TEST_COMB  = [Srel(i); XPYR_dS_TEST(i); XPV_dS_TEST(i); XSOM_dS_TEST(i); XVIP_dS_TEST(i); XPYR_RATE_TEST(i); XPV_RATE_TEST(i); XSOM_RATE_TEST(i); XVIP_RATE_TEST(i); Sirrel(i)].';
       DATA_TRAIN_COMB = [Srel_TRAIN(i,:); XPYR_dS_TRAIN(i,:); XPV_dS_TRAIN(i,:); XSOM_dS_TRAIN(i,:); XVIP_dS_TRAIN(i,:); XPYR_RATE_TRAIN(i,:); XPV_RATE_TRAIN(i,:); XSOM_RATE_TRAIN(i,:); XVIP_RATE_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
              
       DATA_TEST_NULL  = [Srel(i); Sirrel(i)].';
       DATA_TRAIN_NULL = [Srel_TRAIN(i,:); Sirrel_TRAIN(i,:)].';

    
       meansel = mean(DATA_TRAIN,1);
       covsel = cov(DATA_TRAIN);
    
       meansel_null = mean(DATA_TRAIN_NULL,1);
       covsel_null = cov(DATA_TRAIN_NULL);
        
       meansel_RATE = mean(DATA_TRAIN_RATE,1);
       covsel_RATE = cov(DATA_TRAIN_RATE);
       
       meansel_COMB = mean(DATA_TRAIN_COMB,1);
       covsel_COMB = cov(DATA_TRAIN_COMB);
       
       condmean_rel_peers(i) = meansel(1)      + (covsel(1, 2:end)      * inv(covsel(2:end, 2:end)))      * (DATA_TEST(1,2:end)      - meansel(2:end)).';
       condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
       condmean_rel_RATE(i)  = meansel_RATE(1) + (covsel_RATE(1, 2:end) * inv(covsel_RATE(2:end, 2:end))) * (DATA_TEST_RATE(1,2:end) - meansel_RATE(2:end)).';
       condmean_rel_COMB(i)  = meansel_COMB(1) + (covsel_COMB(1, 2:end) * inv(covsel_COMB(2:end, 2:end))) * (DATA_TEST_COMB(1,2:end) - meansel_COMB(2:end)).';      
       
       condcov_rel_peers(i) = covsel(1,1) - covsel(1,2:end) * inv(covsel(2:end,2:end)) * covsel(2:end,1);
       condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
       condcov_rel_RATE(i)  = covsel_RATE(1,1)  - covsel_RATE(1,2:end)  * inv(covsel_RATE(2:end,2:end))  * covsel_RATE(2:end,1);
       condcov_rel_COMB(i)  = covsel_COMB(1,1)  - covsel_COMB(1,2:end)  * inv(covsel_COMB(2:end,2:end))  * covsel_COMB(2:end,1);
       
       LTEST(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_peers(i).')) - 1/2 * (DATA_TEST(:,1)      - condmean_rel_peers(i).').' * (condcov_rel_peers(i).').^(-1) * (DATA_TEST(:,1)      - condmean_rel_peers(i).'));
       LNULL(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
       LRATE(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_RATE(i).'))  - 1/2 * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i).').'  * (condcov_rel_RATE(i).').^(-1)  * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i).'));
       LCOMB(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_COMB(i).'))  - 1/2 * (DATA_TEST_COMB(:,1) - condmean_rel_COMB(i).').'  * (condcov_rel_COMB(i).').^(-1)  * (DATA_TEST_COMB(:,1) - condmean_rel_COMB(i).'));
        
    end
    
end

%% Sample same number of pyramidal cells as SOM, and test predictability per peer cell. 
% For this method, we will take Nsom pyramidal cells and use these as peers
% to fit the remaining Npyr - Nsom pyramidal cell selectivities.

for CellType = 3:5

clear PredPYR
clear PredSOM
clear MDIFF

Nloop = 500;
% CellType = 3; % celltype to make comparisons with
PredPYR = zeros([Nloop, sum(CELLLAB==1) - sum(CELLLAB==CellType)]);
PredSOM = zeros([Nloop, sum(CELLLAB==1) - sum(CELLLAB==CellType)]);

IDS = zeros([Nloop, size(PredPYR,2)]);

for ITS = 1:Nloop


clear XPYR_TRAIN
clear XPYR_RATE_TRAIN
clear XPYR_dS_TRAIN

clear XPV_TRAIN
clear XPV_RATE_TRAIN
clear XPV_dS_TRAIN

clear XSOM_TRAIN
clear XSOM_RATE_TRAIN
clear XSOM_dS_TRAIN

clear XVIP_TRAIN
clear XVIP_RATE_TRAIN
clear XVIP_dS_TRAIN

clear XPYR_TEST
clear XPYR_RATE_TEST
clear XPYR_dS_TEST

clear XPV_TEST
clear XPV_RATE_TEST
clear XPV_dS_TEST

clear XSOM_TEST
clear XSOM_RATE_TEST
clear XSOM_dS_TEST

clear XVIP_TEST
clear XVIP_RATE_TEST
clear XVIP_dS_TEST

clear Srel_TRAIN
clear Sirrel_TRAIN

% select a random subsample of pyramidal cells

PYR_INDS = (CELLLAB == 1);
PYR_INDS0 = find(PYR_INDS);

SAMPLE_PEERS = randperm(length(PYR_INDS0)); % cells whose activity will be used to predict from
SAMPLE_PEERS = SAMPLE_PEERS(1:sum(CELLLAB == CellType)); 

SAMPLE_TEST = setdiff(1:length(PYR_INDS0), SAMPLE_PEERS); % cells whose activity will be predicted

NPYR = length(PYR_INDS0);
NTEST = length(SAMPLE_TEST);  % number of cells left in test sample

IDS(ITS, :) = SAMPLE_TEST; % identities of cells
  
Srel_TEST = Srel(SAMPLE_TEST);
Sirrel_TEST = Sirrel(SAMPLE_TEST);

MAT = MATrel;

parfor i=1:NTEST
    
    % training datasets - remove one pyramidal cell and compute prediction metrics
    
    ind0 = PYR_INDS0;
    
    ind1_PEERS = SAMPLE_PEERS;
    ind_PEERS  = ind0(ind1_PEERS); 
    
    ind1_TEST = SAMPLE_TEST([1:(i-1), (i+1):length(SAMPLE_TEST)]);  % leave out ith cell
    ind_TEST  = ind0(ind1_TEST);    % indexing over all cell types
 
    M = MAT(ind_TEST, ind_PEERS);
    XPYR_TRAIN(i,:) = M * Srel(ind1_PEERS).' / length(ind1_PEERS);
    XPYR_RATE_TRAIN(i,:) = M * dFPYR(ind1_PEERS).' / length(ind1_PEERS); % firing rate of PYR multiplied by noise correlations
    XPYR_dS_TRAIN(i,:) = M * dS_PYR(ind1_PEERS).' / length(ind1_PEERS);
   
   
    M = MAT(ind_TEST, CELLLAB==3);
    XPV_TRAIN(i,:) = M * SrelPV.' / length(SrelPV);
    XPV_RATE_TRAIN(i,:) = M * dFPV.' / length(dFPV); % firing rate of PV multiplied by noise correlations
    XPV_dS_TRAIN(i,:) = M * dS_PV.' / length(dS_PV);

    
    M = MAT(ind_TEST, CELLLAB==4);
    XSOM_TRAIN(i,:) = M * SrelSOM.' / length(SrelSOM);
    XSOM_RATE_TRAIN(i,:) = M * dFSOM.' / length(dFSOM); % firing rate of PV multiplied by noise correlations
    XSOM_dS_TRAIN(i,:) = M * dS_SOM.' / length(dS_SOM);

    
    M = MAT(ind_TEST, CELLLAB==5);
    XVIP_TRAIN(i,:) = M * SrelVIP.' / length(SrelVIP);
    XVIP_RATE_TRAIN(i,:) = M* dFVIP.' / length(dFVIP); % firing rate of PV multiplied by noise correlations
    XVIP_dS_TRAIN(i,:) = M * dS_VIP.' / length(dS_VIP);
       
    Srel_TRAIN(i,:) = Srel(ind1_TEST);
    Sirrel_TRAIN(i,:) = Sirrel(ind1_TEST);
    
    % test datasets - use only cell previously left out, and predict from peer dataset
    
    ind1_TEST = SAMPLE_TEST;  % test set indices, pyramidal 
    ind_TEST  = ind0(ind1_TEST);    % test set indices in full population
    
    M = MAT(ind_TEST, ind_PEERS);
    XPYR_TEST(i) = M(i,:) * Srel(ind1_PEERS).' / length(ind1_PEERS);
    XPYR_RATE_TEST(i) = M(i,:) * dFPYR(ind1_PEERS).' / length(ind1_PEERS);
    XPYR_dS_TEST(i) = M(i,:) * dS_PYR(ind1_PEERS).' / length(ind1_PEERS);

    
    M = MAT(ind_TEST, CELLLAB==3);
    XPV_TEST(i) = M(i,:) * SrelPV.' / length(SrelPV);
    XPV_RATE_TEST(i) = M(i,:) * dFPV.' / length(dFPV);
    XPV_dS_TEST(i) = M(i,:) * dS_PV.' / length(dS_PV);
        
    M = MAT(ind_TEST, CELLLAB==4);
    XSOM_TEST(i) = M(i,:) * SrelSOM.' / length(SrelSOM);
    XSOM_RATE_TEST(i) = M(i,:) * dFSOM.' / length(dFSOM);
    XSOM_dS_TEST(i) = M(i,:) * dS_SOM.' / length(dS_SOM);

       
    M = MAT(ind_TEST, CELLLAB==5);
    XVIP_TEST(i) = M(i,:) * SrelVIP.' / length(SrelVIP);
    XVIP_RATE_TEST(i) = M(i,:) * dFVIP.' / length(dFVIP);
    XVIP_dS_TEST(i) = M(i,:) * dS_VIP.' / length(dS_VIP);

   
      
end

    clear condmean_rel_peers
    clear condmean_rel_null
    clear condmean_rel_COMB
    clear condmean_rel_RATE
    clear condcov_rel_peers
    clear condcov_rel_null
    clear condcov_rel_COMB
    clear condcov_rel_RATE
    clear LNULL
    clear LTEST 
    clear LRATE
    clear LCOMB
    
     
 for i=1:NTEST
            
   
       DATA_TEST_RATE  = [Srel_TEST(i); XPYR_dS_TEST(i); Sirrel_TEST(i)].';
       DATA_TRAIN_RATE = [Srel_TRAIN(i,:); XPYR_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
  
       if CellType == 3
       
         DATA_TEST_COMB  = [Srel_TEST(i); XPV_dS_TEST(i); Sirrel_TEST(i)].';
         DATA_TRAIN_COMB = [Srel_TRAIN(i,:); XPV_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
           
       elseif CellType == 4
       
         DATA_TEST_COMB  = [Srel_TEST(i); XSOM_dS_TEST(i); Sirrel_TEST(i)].';
         DATA_TRAIN_COMB = [Srel_TRAIN(i,:); XSOM_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
       
       elseif CellType == 5
       
         DATA_TEST_COMB  = [Srel_TEST(i); XVIP_dS_TEST(i); Sirrel_TEST(i)].';
         DATA_TRAIN_COMB = [Srel_TRAIN(i,:); XVIP_dS_TRAIN(i,:); Sirrel_TRAIN(i,:)].';
           
       end    
       
           
       DATA_TEST_NULL  = [Srel_TEST(i); Sirrel_TEST(i)].';
       DATA_TRAIN_NULL = [Srel_TRAIN(i,:); Sirrel_TRAIN(i,:)].';

    
      
       meansel_null = mean(DATA_TRAIN_NULL,1);
       covsel_null = cov(DATA_TRAIN_NULL);
        
       meansel_RATE = mean(DATA_TRAIN_RATE,1);
       covsel_RATE = cov(DATA_TRAIN_RATE);
       
       meansel_COMB = mean(DATA_TRAIN_COMB,1);
       covsel_COMB = cov(DATA_TRAIN_COMB);
       
       condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
       condmean_rel_RATE(i)  = meansel_RATE(1) + (covsel_RATE(1, 2:end) * inv(covsel_RATE(2:end, 2:end))) * (DATA_TEST_RATE(1,2:end) - meansel_RATE(2:end)).';
       condmean_rel_COMB(i)  = meansel_COMB(1) + (covsel_COMB(1, 2:end) * inv(covsel_COMB(2:end, 2:end))) * (DATA_TEST_COMB(1,2:end) - meansel_COMB(2:end)).';      
       
       condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end)   * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
       condcov_rel_RATE(i)  = covsel_RATE(1,1)  - covsel_RATE(1,2:end)  * inv(covsel_RATE(2:end,2:end))  * covsel_RATE(2:end,1);
       condcov_rel_COMB(i)  = covsel_COMB(1,1)  - covsel_COMB(1,2:end)  * inv(covsel_COMB(2:end,2:end))  * covsel_COMB(2:end,1);
       
       LNULL(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
       LRATE(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_RATE(i).'))  - 1/2 * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i).').'  * (condcov_rel_RATE(i).').^(-1)  * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i).'));
       LCOMB(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_COMB(i).'))  - 1/2 * (DATA_TEST_COMB(:,1) - condmean_rel_COMB(i).').'  * (condcov_rel_COMB(i).').^(-1)  * (DATA_TEST_COMB(:,1) - condmean_rel_COMB(i).'));
        
 end
  
 PredPYR(ITS,:) = LRATE - LNULL;
 PredSOM(ITS,:) = LCOMB - LNULL;
 
end

MDIFF = zeros([1, NPYR]);

for i=1:sum(CELLLAB == 1)

    MDIFF(i) = median(PredSOM(IDS == i) - PredPYR(IDS == i));

end

[h(CellType), p(CellType)] = signtest(MDIFF);  % test whether interneurons predict better than pyramidal cells (per cell)
FracCells(CellType) = sum(MDIFF > 0)/length(MDIFF); % fraction of cells which are better predicted by interneuron than pyramidal

clear MDIFF;

end

%% Try more complex gating models with gating values for each SOM cell 

% first order model - one gating value for each SOM, irrespective of PYR cells
        
CellNum = 4; % which kind of cell is performing gating?

if CellNum == 3
    
    dFInt = dFPV;

elseif CellNum == 4
    
    dFInt = dFSOM;
    
elseif CellNum == 5
    
    dFInt = dFVIP;
    
end
    
MAT = MATrel;

M = MAT(CELLLAB==1, CELLLAB==CellNum);

NPYR = sum(CELLLAB == 1);


clear condmean_rel_peers
clear condmean_rel_null
clear condmean_rel_RATE
clear LNULL
clear LRATE
clear XSOM_RATES

N_COMB = 2^sum(CELLLAB==CellNum); % number of possible gating combinations
gating_vectors = dec2bin(0:(N_COMB-1)) - '0'; % all possible gating vectors, stored as matrix
gating_vectors = gating_vectors * 2 - 1;

DATA_NULL = [Srel; Sirrel].';
% 
%  for i=1:size(M,1)
%      for j=1:size(M,2)
%     MM = M([1:(i-1), (i+1):NPYR],j);
%     M([1:(i-1), (i+1):NPYR],j) = MM(randperm(length(MM)));
%      end
%  end

for N=1:N_COMB
    
   XSOM_RATES(N,:) = M * (dFInt.* gating_vectors(N, :)).' / length(dFInt); % firing rate of SOM multiplied by noise correlations and gating vector

   DATA_RATE = [Srel; XSOM_RATES(N,:); Sirrel].';
 

    for i=1:NPYR
    
       DATA_TEST_RATE = DATA_RATE(i,:);
       DATA_TRAIN_RATE = DATA_RATE([1:(i-1), (i+1):end],:);
       
       DATA_TEST_NULL = DATA_NULL(i,:);
       DATA_TRAIN_NULL = DATA_NULL([1:(i-1), (i+1):end],:);
    
       meansel_null = mean(DATA_TRAIN_NULL,1);
       covsel_null = cov(DATA_TRAIN_NULL);
        
       meansel_RATE = mean(DATA_TRAIN_RATE,1);
       covsel_RATE = cov(DATA_TRAIN_RATE);
       
       condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
       condmean_rel_RATE(i, N)  = meansel_RATE(1) + (covsel_RATE(1, 2:end) * inv(covsel_RATE(2:end, 2:end))) * (DATA_TEST_RATE(1,2:end) - meansel_RATE(2:end)).';

       condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
       condcov_rel_RATE(i)  = covsel_RATE(1,1)  - covsel_RATE(1,2:end)  * inv(covsel_RATE(2:end,2:end))  * covsel_RATE(2:end,1);

       LNULL(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
       LRATE(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_RATE(i).'))  - 1/2 * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).').'  * (condcov_rel_RATE(i).').^(-1)  * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).'));
        
    end

end


%% Try gating model using noise correlations different between vertical and angled gratings as gating variable instead

CellNum = 5; % which kind of cell is performing gating?

COND = 1; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

MATV = MAT;
MV = MATV(CELLLAB==1, CELLLAB == CellNum);

COND = 2; % experimental condition

MAT = ADDTRG{1}.PL{1}.DCOR{COND}.CORNSE1;  % noise correlations for grey corridor condition
MAT(isnan(MAT)) = 0;
MAT(logical(eye(size(MAT)))) = 1;
MAT = MAT + triu(MAT,1).';  %% symmetrise matrix

MATA = MAT;
MA = MATA(CELLLAB==1, CELLLAB == CellNum);

%GATEVAR = (MV - MA)./(abs(MV) + abs(MA));
GATEVAR = (MV - MA);


dF = (m_vert_rel + m_ang_rel - (m_vert_irrel + m_ang_irrel));
dF = dF(accept);
dFPYR = dF(CELLLAB == 1); % Fractional change in PYR firing between conditions
dFPV = dF(CELLLAB == 3); % Fractional change in PV firing between conditions
dFSOM = dF(CELLLAB == 4); % Fractional change in SOM firing between conditions
dFVIP = dF(CELLLAB == 5); % Fractional change in VIP firing between conditions


if CellNum == 3
    
    dFInt = dFPV;

elseif CellNum == 4
    
    dFInt = dFSOM;
    
elseif CellNum == 5
    
    dFInt = dFVIP;
    
end
    
MAT = MATrel;

M = MAT(CELLLAB==1, CELLLAB==CellNum);

NPYR = sum(CELLLAB == 1);

clear condmean_rel_peers
clear condmean_rel_null
clear condmean_rel_RATE
clear LNULL
clear LRATE
clear XSOM_RATES

Nstep = 250;

DATA_NULL = [Srel; Sirrel].';

for N=1:Nstep
     
 MidPoint =  min(mean(GATEVAR,1)) + (N-1) * range(mean(GATEVAR,1)) / (Nstep - 1);
 %MidPoint = -1 + (N-1) * 2/(Nstep-1);
 GATEVAR_CENTRED = mean(GATEVAR,1) - MidPoint;
 XSOM_RATES(N,:) = M * (GATEVAR_CENTRED .* dFInt).' / length(dFInt); % firing rate of SOM multiplied by noise correlations and gating vector

%MidPoint = min(min(GATEVAR)) + (N-1) * (max(max(GATEVAR)) - min(min(GATEVAR))) / (Nstep - 1);
%GATEVAR_CENTRED = GATEVAR - MidPoint;
%XSOM_RATES(N,:) = (M .* (GATEVAR_CENTRED)) * dFInt.' / length(dFInt); % firing rate of SOM multiplied by noise correlations and gating matrix



 DATA_RATE = [Srel; XSOM_RATES(N,:); Sirrel].';
     
 
    for i=1:NPYR
    
       DATA_TEST_RATE = DATA_RATE(i,:);
       DATA_TRAIN_RATE = DATA_RATE([1:(i-1), (i+1):end],:);
       
       DATA_TEST_NULL = DATA_NULL(i,:);
       DATA_TRAIN_NULL = DATA_NULL([1:(i-1), (i+1):end],:);
    
       meansel_null = mean(DATA_TRAIN_NULL,1);
       covsel_null = cov(DATA_TRAIN_NULL);
        
       meansel_RATE = mean(DATA_TRAIN_RATE,1);
       covsel_RATE = cov(DATA_TRAIN_RATE);
       
       condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
       condmean_rel_RATE(i, N)  = meansel_RATE(1) + (covsel_RATE(1, 2:end) * inv(covsel_RATE(2:end, 2:end))) * (DATA_TEST_RATE(1,2:end) - meansel_RATE(2:end)).';

       condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
       condcov_rel_RATE(i)  = covsel_RATE(1,1)  - covsel_RATE(1,2:end)  * inv(covsel_RATE(2:end,2:end))  * covsel_RATE(2:end,1);

       LNULL(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
       LRATE(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_RATE(i).'))  - 1/2 * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).').'  * (condcov_rel_RATE(i).').^(-1)  * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).'));
        
    end

end
% 
% 
% % Second order model - one gating for each SOM and PYR pair, reflecting gating of a specific presynaptic pyramidal cell onto all other pyramidal cells
%             
% MSOM = MAT(CELLLAB==1, CELLLAB==4);
% MPYR = MAT(CELLLAB==1, CELLLAB==1);
% 
% NPYR = sum(CELLLAB == 1);
% 
% DATA = [Srel; XPYR.'; XPV.'; XSOM.'; XVIP.'; Sirrel].';
% DATA_NULL = [Srel; Sirrel].';
% 
% clear condmean_rel_peers
% clear condmean_rel_nullc
% clear LTEST
% clear LNULL
% clear LRATE
% 
% 
% for i=1:NPYR
%     
%        DATA_TEST_NULL = DATA_NULL(i,:);
%        DATA_TRAIN_NULL = DATA_NULL([1:(i-1), (i+1):end],:);
%     
%        meansel_null = mean(DATA_TRAIN_NULL,1);
%        covsel_null = cov(DATA_TRAIN_NULL);
%         
%        condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
% 
%        condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
% 
%        LNULL(i) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
% 
%         
% end
% 
% for N=1:sum(CELLLAB == 1)
%     for M=1:sum(CELLLAB == 4)
%         
%         gating_vectors = zeros([sum(CELLLAB == 4), sum(CELLLAB==1)]);
%         
%         gating_vectors(N,M) = 1;
%               
%         XSOM_RATES = MSOM(:,M) .* dFSOM(M) .* (MPYR(:,N) .* Srel(N)); % 2nd order gating measure
%        
%         DATA_RATE = [Srel; XSOM_RATES.'; Sirrel].';
%     
%               
%     for i=1:NPYR
%     
%        DATA_TEST_RATE = DATA_RATE(i,:);
%        DATA_TRAIN_RATE = DATA_RATE([1:(i-1), (i+1):end],:);
%        
%        DATA_TEST_NULL = DATA_NULL(i,:);
%        DATA_TRAIN_NULL = DATA_NULL([1:(i-1), (i+1):end],:);
%     
%        meansel_null = mean(DATA_TRAIN_NULL,1);
%        covsel_null = cov(DATA_TRAIN_NULL);
%         
%        meansel_RATE = mean(DATA_TRAIN_RATE,1);
%        covsel_RATE = cov(DATA_TRAIN_RATE);
%        
%        condmean_rel_null(i)  = meansel_null(1) + (covsel_null(1, 2:end) * inv(covsel_null(2:end, 2:end))) * (DATA_TEST_NULL(1,2:end) - meansel_null(2:end)).';
%        condmean_rel_RATE(i,N)  = meansel_RATE(1) + (covsel_RATE(1, 2:end) * inv(covsel_RATE(2:end, 2:end))) * (DATA_TEST_RATE(1,2:end) - meansel_RATE(2:end)).';
% 
%        condcov_rel_null(i)  = covsel_null(1,1) - covsel_null(1,2:end) * inv(covsel_null(2:end,2:end)) * covsel_null(2:end,1);
%        condcov_rel_RATE(i)  = covsel_RATE(1,1)  - covsel_RATE(1,2:end)  * inv(covsel_RATE(2:end,2:end))  * covsel_RATE(2:end,1);
% 
%        LNULL(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_null(i).'))  - 1/2 * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).').'  * (condcov_rel_null(i).').^(-1)  * (DATA_TEST_NULL(:,1) - condmean_rel_null(i).'));
%        LRATE(i, N) = sum(-1/2 * log(2*pi * abs(condcov_rel_RATE(i).'))  - 1/2 * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).').'  * (condcov_rel_RATE(i).').^(-1)  * (DATA_TEST_RATE(:,1) - condmean_rel_RATE(i,N).'));
%         
%     end
%     
%     end
% end

%% Next plot the changes in selectivity under the null model and compare to data for angled and vertical preferring cells.

accept = (abs(Srel - mean(Srel)) < 2*std(Srel) & abs(Sirrel - mean(Sirrel)) < 2*std(Sirrel));

Srel_cleaned = Srel(accept);
Sirrel_cleaned = Sirrel(accept);

meansel = [mean(Srel_cleaned), mean(Sirrel_cleaned)];
covsel = cov(Srel_cleaned, Sirrel_cleaned);

condmean_rel_null = meansel(1) + std(Srel_cleaned)/std(Sirrel_cleaned) * corr(Srel_cleaned.', Sirrel_cleaned.') * (Sirrel_cleaned - mean(Sirrel_cleaned)); % conditional mean of Srel under model given data for Sirrel
condvar_rel_null = covsel(1,1) - covsel(1,2) * (1/covsel(2,2)) * covsel(2,1);
Srel_null = sqrt(condvar_rel_null) .* randn(size(condmean_rel_null)) + condmean_rel_null;
Sirrel_null = Sirrel_cleaned;
dS_null = Srel_null - Sirrel_null;

dX = diff(X);
Xnull = min(dS_null):dX(1):max(dS_null); % same bins as for original plot

subplot(221)

hist(dS_null(Sirrel_null < 0), Xnull)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with negative selectivity in irrelevant condition')

subplot(223)

hist(dS_null(Sirrel_null > 0), Xnull)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with positive selectivity in irrelevant condition')

subplot(222)

hist(dS_null(Srel_null < 0), Xnull)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with negative selectivity in relevant condition')

subplot(224)

hist(dS_null(Srel_null > 0), Xnull)

axis([-5 5 0 45])

set(gca, 'fontsize', 18)
xlabel('Change in selectivity')
ylabel('Number of cells')
title('Cells with positive selectivity in relevant condition')

