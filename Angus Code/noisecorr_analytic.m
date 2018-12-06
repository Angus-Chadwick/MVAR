function [ctot_pre_allweights, ctot_post_allweights, ctot_pre_manipW1W2, ctot_post_manipW1W2] = noisecorr_analytic(manipulation, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL)

% Inputs:
% manipulation: removeW1W2 or swapW1W2
% W1 = postsynaptic neuron type
% W2 = presynaptic neuron type
% Outputs:
% ctot_pre_allweights: total intact noise correlation
% ctot

RXS = [1,3,4,5,6,8,9,10];
CLS = [1,3,4,5];


%% Case 1: Remove weights between specific populations

if strmatch(manipulation, 'removeW1W2')

for removeW1W2 = 0:1

for prepost = {'pre', 'post'}
        
corr_tot0 = cell(4);
        
for rix = 1:length(RXS)

if strmatch(prepost, 'pre')
    
        vel_summand = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;

        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpre{RXS(rix)}(:,end);
        Vdata = vmatpre0{RXS(rix), 1};
        residuals = ResidualV_pre_TOT{RXS(rix)};
        rmat = rmat_pre{RXS(rix),1};
        Trials = Trialspre_all{RXS(rix)}{1};
    
elseif strmatch(prepost, 'post')
    
        vel_summand = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;

        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpost{RXS(rix)}(:,end);
        Vdata = vmatpost0{RXS(rix), 1};
        residuals = ResidualV_post_TOT{RXS(rix)};
        rmat = rmat_post{RXS(rix),1};
        Trials = Trialspost_all{RXS(rix)}{1};
end
          
        
    
if removeW1W2
        
        Wmat_diag = diag(diag(Wmat));
        Wmat(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2) = 0;
        Wmat = Wmat - diag(diag(Wmat)) + Wmat_diag;

end


for k = 0:length(tsims)
    
    Wmat_pow{k + 1} = Wmat^k;  % note the index offset to allow for a zero power
    
end

for t=1:length(tsims)
   
                
        for i=1:t
            
            vel_summand(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * (permute(Vdata(:, tsims(i)-1) - mean(Vdata(:, tsims(i)-1)), [2,1]));
            res_summand(:,:,t,i) = Wmat_pow{t-i+1} * (permute(residuals(:,tsims(i)-1, :) - repmat(mean(residuals(:,tsims(i)-1, :),3), [1,1,length(Trials)]), [1,3,2]));

        end 
    
        
end
    
% 
% for t=1:length(tsims)
%     for T=1:length(Trials)
%         for i=1:t
% 
%             vel_summand(:,T,t,i) = Wmat^(t-i) * Vcoefs * (Vdata(T, tsims(i)-1) - mean(Vdata(:, tsims(i)-1)));
%             res_summand(:,T,t,i) = Wmat^(t-i) * squeeze(residuals(:,tsims(i)-1, T) - mean(residuals(:,tsims(i)-1, :),3));
% 
%         end 
%     end
% end

% next execute sums to get quantities which determine covariance
    
cov_vel0 = sum(vel_summand,4);
cov_res0 = sum(res_summand,4);
cov_initstate0 = zeros([size(cov_vel0,1), size(cov_vel0,2), length(tsims)]);

rmat_perm = permute(rmat, [2,1,3]);

for t=1:length(tsims)

         cov_initstate0(:,:, t) = Wmat_pow{t + 1} * (rmat_perm(:, Trials, rngtot(tsims(1))-1) - repmat(nanmean(rmat_perm(:, Trials, rngtot(tsims(1))-1),2), [1, length(Trials), 1]));

end


% for t=1:length(tsims)
%     for T=1:length(Trials)
% 
%          cov_initstate0(:,T, t) = Wmat^t * (rmat(Trials(T), :, rngtot(tsims(1))-1) - mean(rmat(Trials, :, rngtot(tsims(1))-1))).';
% 
%     end   
% end

% next concatenate samples and trials into a single vector and compute covariance matrices

Mvel = permute(cov_vel0, [ 3, 2, 1]);
Mres = permute(cov_res0, [ 3, 2, 1]);
Minit = permute(cov_initstate0, [ 3, 2, 1]);

Mtot = Mvel + Mres + Minit;

cov_tot = cov(reshape(Mtot, [size(Mtot,1) * size(Mtot,2), size(Mtot,3)]));

D = diag(diag(cov_tot));

% transform to find correlation matrices

corr_tot{rix} = D^(-1/2) * cov_tot * D^(-1/2);

for i=1:length(CLS)
    for j=1:length(CLS)
        
      Ctot  = corr_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      
      if i==j
      
            Ctot(1:size(Ctot,1)+1:end) = nan;

      end
            
      corr_tot0{i,j} = vertcat(corr_tot0{i,j}, Ctot(:));
        
    end
end
        
end

indsx = [1,4,3,2,2,1,2,1,1,3];
indsy = [1,4,3,2,4,2,3,4,3,4];

if strmatch(prepost, 'pre')
        
ctot_pre = corr_tot0;
  

elseif strmatch(prepost, 'post')

        
ctot_post = corr_tot0;
        

end

end

 for i=1:length(indsx)

    ctot_pre0{i} = ctot_pre{indsx(i), indsy(i)};
    ctot_post0{i} = ctot_post{indsx(i), indsy(i)};

 end
 


if removeW1W2 == 0

ctot_pre_allweights = ctot_pre0;
ctot_post_allweights = ctot_post0;
   
elseif removeW1W2
    
ctot_pre_manipW1W2 = ctot_pre0;
ctot_post_manipW1W2 = ctot_post0;

    
end

end


%% Case 2: Swap weights pre/post between specific populations


elseif strmatch(manipulation, 'swapW1W2')

for swapW1W2 = 0:1

for prepost = {'pre', 'post'}
        
corr_tot0 = cell(4);
        
for rix = 1:length(RXS)

vel_summand = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
res_summand = vel_summand;

if strmatch(prepost, 'pre')

        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpre{RXS(rix)}(:,end);
        Vdata = vmatpre0{RXS(rix), 1};
        residuals = ResidualV_pre_TOT{RXS(rix)};
        rmat = rmat_pre{RXS(rix),1};
        Trials = Trialspre_all{RXS(rix)}{1};
    
elseif strmatch(prepost, 'post')

        Wmat = (eye(size(rmat_post{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpost{RXS(rix)}(:,end);
        Vdata = vmatpost0{RXS(rix), 1};
        residuals = ResidualV_post_TOT{RXS(rix)};
        rmat = rmat_post{RXS(rix),1};
        Trials = Trialspost_all{RXS(rix)}{1};
end
   
if swapW1W2 
    
    if strmatch(prepost, 'pre')
                
            Wmat_diag = diag(Wmat);
            
             Wmat(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2) = xpost{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2);
             
             for i=1:length(Wmat)
                 
                 Wmat(i,i) = Wmat_diag(i);
                 
             end
                              
        
        % uncomment the below section to swap all weights while preserving
        % the self-interaction terms
        
%                  Wmat_diag = diag(Wmat);
%                   Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
% 
%                 for i=1:length(Wmat)
%                  
%                       Wmat(i,i) = Wmat_diag(i);
%                  
%                 end              




    elseif strmatch(prepost, 'post')
        
        
              Wmat_diag = diag(Wmat);
            
              Wmat(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2) = xpre{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2);
              
             for i=1:length(Wmat)
                 
                 Wmat(i,i) = Wmat_diag(i);
                 
             end           
              
        
        % uncomment the below section to swap all weights while preserving
        % the self-interaction terms
        
%                 Wmat_diag = diag(Wmat);
%                  Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
%                      
%                 for i=1:length(Wmat)
%                  
%                        Wmat(i,i) = Wmat_diag(i);
%                  
%                 end
                   
    end
   
    
end
                
        


for t=1:length(tsims)
    for T=1:length(Trials)
        for i=1:t

            vel_summand(:,T,t,i) = Wmat^(t-i) * Vcoefs * (Vdata(T, tsims(i)-1) - mean(Vdata(:, tsims(i)-1)));
            res_summand(:,T,t,i) = Wmat^(t-i) * squeeze(residuals(:,tsims(i)-1, T) - mean(residuals(:,tsims(i)-1, :),3));

        end 
    end
end

% next execute sums to get quantities which determine covariance
      
cov_vel0 = sum(vel_summand,4);
cov_res0 = sum(res_summand,4);
cov_initstate0 = zeros([size(cov_vel0,1), size(cov_vel0,2), length(tsims)]);

for t=1:length(tsims)
    for T=1:length(Trials)

         cov_initstate0(:,T, t) = Wmat^t * (rmat(Trials(T), :, rngtot(tsims(1))-1) - mean(rmat(Trials, :, rngtot(tsims(1))-1))).';

    end   
end

Mvel = permute(cov_vel0, [ 3, 2, 1]);
Mres = permute(cov_res0, [ 3, 2, 1]);
Minit = permute(cov_initstate0, [ 3, 2, 1]);

Mtot = Mvel + Mres + Minit;

cov_tot = cov(reshape(Mtot, [size(Mtot,1) * size(Mtot,2), size(Mtot,3)]));
D = diag(diag(cov_tot));

% transform to find correlation matrices

corr_tot{rix} = D^(-1/2) * cov_tot * D^(-1/2);

for i=1:length(CLS)
    for j=1:length(CLS)
        
      Ctot  = corr_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      
      if i==j
      
            Ctot(1:size(Ctot,1)+1:end) = nan;

      end
            
        corr_tot0{i,j} = vertcat(corr_tot0{i,j}, Ctot(:));
        
    end
end
        
end

indsx = [1,4,3,2,2,1,2,1,1,3];
indsy = [1,4,3,2,4,2,3,4,3,4];

if strmatch(prepost, 'pre')
        
ctot_pre = corr_tot0;
  

elseif strmatch(prepost, 'post')

        
ctot_post = corr_tot0;
        

end

end

 for i=1:length(indsx)

    ctot_pre0{i} = ctot_pre{indsx(i), indsy(i)};
    ctot_post0{i} = ctot_post{indsx(i), indsy(i)};

 end
 


if swapW1W2 == 0

ctot_pre_allweights = ctot_pre0;
ctot_post_allweights = ctot_post0;

elseif swapW1W2
  
ctot_pre_manipW1W2 = ctot_pre0;
ctot_post_manipW1W2 = ctot_post0;
    
end

end

%% Case 3: Swap residuals pre/post

elseif strmatch(manipulation, 'swapresids')

for swapresids = 0:1

for prepost = {'pre', 'post'}
        
corr_tot0 = cell(4);
       

for rix = 1:length(RXS)

if strmatch(prepost, 'pre')
    
        vel_summand = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;

        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpre{RXS(rix)}(:,end);
        Vdata = vmatpre0{RXS(rix), 1};
        rmat = rmat_pre{RXS(rix),1};
        Trials = Trialspre_all{RXS(rix)}{1};
    
elseif strmatch(prepost, 'post')
    
        vel_summand = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;

        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpost{RXS(rix)}(:,end);
        Vdata = vmatpost0{RXS(rix), 1};
        rmat = rmat_post{RXS(rix),1};
        Trials = Trialspost_all{RXS(rix)}{1};
end


residuals_pre0 = ResidualV_pre_TOT{RXS(rix)};
residuals_post0 = ResidualV_post_TOT{RXS(rix)};

if swapresids  % generate additional residual samples with same mean and covariance

ResMean_pre = nanmean(residuals_pre0,3);
ResMean_post = nanmean(residuals_post0,3);

ResCov_pre = zeros([size(residuals_pre0,1), size(residuals_pre0,1), 17]);
ResCov_post = ResCov_pre;

for t=1:17
    
ResCov_pre(:,:,t) = nancov(squeeze(residuals_pre0(:,t,:))');
ResCov_post(:,:,t) = nancov(squeeze(residuals_post0(:,t,:))');

end

Trialsdiff = length(Trialspost_all{RXS(rix)}{1}) - length(Trialspre_all{RXS(rix)}{1});

if Trialsdiff > 0
    
    for T=1:abs(Trialsdiff)   
        for t=1:17
        
             residuals_pre0(:,t,length(Trialspre_all{RXS(rix)}{1}) + T) =  ResMean_pre(:,t) +  sqrtm(squeeze(ResCov_pre(:,:,t))) * rand([size(ResMean_pre, 1),1]);

        end
    end
    
    residuals_post = residuals_pre0;  % swap pre and post
    residuals_pre = residuals_post0(:,:, 1:length(Trialspre_all{RXS(rix)}{1}));
    
elseif Trialsdiff < 0
    
    for T=1:abs(Trialsdiff)    
        for t=1:17
            
            residuals_post0(:,t,length(Trialspost_all{RXS(rix)}{1}) + T) =  ResMean_post(:,t) + sqrtm(squeeze(ResCov_post(:,:,t))) * rand([size(ResMean_post, 1),1]);
       
        end
    end
    
    residuals_post = residuals_pre0(:,:, 1:length(Trialspost_all{RXS(rix)}{1}));  % swap pre and post
    residuals_pre = residuals_post0;
    
end

else 
    
    residuals_post = residuals_post0;
    residuals_pre = residuals_pre0;
    
end
   
if strmatch(prepost, 'pre')
    
    residuals = residuals_pre;
    
elseif strmatch(prepost, 'post')
    
    residuals = residuals_post;
    
end


for k = 0:length(tsims)
    
    Wmat_pow{k + 1} = Wmat^k;  % note the index offset to allow for a zero power
    
end


    for t=1:length(tsims)   
                
        for i=1:t
            
            vel_summand(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * permute(Vdata(:, tsims(i)-1) - mean(Vdata(:, tsims(i)-1)), [2,1]);
            res_summand(:,:,t,i) = Wmat_pow{t-i+1} * permute(residuals(:,tsims(i)-1, :) - repmat(mean(residuals(:,tsims(i)-1, :),3), [1,1,length(Trials)]), [1,3,2]);

        end    
        
    end
    
    
% next execute sums to get quantities which determine covariance
    
cov_vel0 = sum(vel_summand,4);
cov_res0 = sum(res_summand,4);
cov_initstate0 = zeros([size(cov_vel0,1), size(cov_vel0,2), length(tsims)]);

rmat_perm = permute(rmat, [2,1,3]);

for t=1:length(tsims)

         cov_initstate0(:,:, t) = Wmat_pow{t + 1} * (rmat_perm(:, Trials, rngtot(tsims(1))-1) - repmat(nanmean(rmat_perm(:, Trials, rngtot(tsims(1))-1),2), [1, length(Trials), 1]));

end


% next concatenate samples and trials into a single vector and compute covariance matrices

Mvel = permute(cov_vel0, [ 3, 2, 1]);
Mres = permute(cov_res0, [ 3, 2, 1]);
Minit = permute(cov_initstate0, [ 3, 2, 1]);

Mtot = Mvel + Mres + Minit;

cov_tot = cov(reshape(Mtot, [size(Mtot,1) * size(Mtot,2), size(Mtot,3)]));

D = diag(diag(cov_tot));

% transform to find correlation matrices

corr_tot{rix} = D^(-1/2) * cov_tot * D^(-1/2);

for i=1:length(CLS)
    for j=1:length(CLS)
        
      Ctot  = corr_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      
      if i==j
      
            Ctot(1:size(Ctot,1)+1:end) = nan;

      end
            
      corr_tot0{i,j} = vertcat(corr_tot0{i,j}, Ctot(:));
        
    end
end
        
end

indsx = [1,4,3,2,2,1,2,1,1,3];
indsy = [1,4,3,2,4,2,3,4,3,4];

if strmatch(prepost, 'pre')
        
ctot_pre = corr_tot0;
  

elseif strmatch(prepost, 'post')

        
ctot_post = corr_tot0;
        

end

end

 for i=1:length(indsx)

    ctot_pre0{i} = ctot_pre{indsx(i), indsy(i)};
    ctot_post0{i} = ctot_post{indsx(i), indsy(i)};

 end
 


if swapresids == 0

ctot_pre_allweights = ctot_pre0;
ctot_post_allweights = ctot_post0;
   
elseif swapresids
    
ctot_pre_manipW1W2 = ctot_pre0;
ctot_post_manipW1W2 = ctot_post0;

    
end

end

end

