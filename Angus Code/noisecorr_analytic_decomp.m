function [ctot_pre, ctot_post, ctot_pre_shuffres, ctot_post_shuffres] = noisecorr_analytic_decomp(removeweights, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL)

 
RXS = [1,3,4,5,6,8,9,10];
CLS = [1,3,4,5];

for prepost = {'pre', 'post'}
        
corr_vel_tot0 = cell(4);
corr_res_tot0 = cell(4);
corr_res_shuff_tot0 = cell(4);
corr_initstate_tot0 = cell(4);

corr_vel_init_tot0 = cell(4);
corr_vel_res_tot0 = cell(4);
corr_res_init_tot0 = cell(4);
corr_vel_res_shuff_tot0 = cell(4);
corr_res_shuff_init_tot0 = cell(4);      

for rix = 1:length(RXS)

if strmatch(prepost, 'pre')

        vel_summand = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;
        res_shuff_summand = res_summand;
        
        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpre{RXS(rix)}(:,end);
        Vdata = vmatpre0{RXS(rix), 1};
        residuals = ResidualV_pre_TOT{RXS(rix)};
        rmat = rmat_pre{RXS(rix),1};
        Trials = Trialspre_all{RXS(rix)}{1};
    
elseif strmatch(prepost, 'post')

        vel_summand = zeros([size(xpre{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summand = vel_summand;
        res_shuff_summand = res_summand;
        
        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpost{RXS(rix)}(:,end);
        Vdata = vmatpost0{RXS(rix), 1};
        residuals = ResidualV_post_TOT{RXS(rix)};
        rmat = rmat_post{RXS(rix),1};
        Trials = Trialspost_all{RXS(rix)}{1};
end
   
                
        
if removeweights
   
        Wmat = diag(diag(Wmat));      

end

    
residuals_shuff = zeros(size(residuals));
    
for i=1:size(residuals,1)
    
        residuals_shuff(i,:,:) = residuals(i,:,randperm(size(residuals,3)));
    
end
    


for k = 0:length(tsims)
    
    Wmat_pow{k + 1} = Wmat^k;  % note the index offset to allow for a zero power
    
end


for t=1:length(tsims)
   
                
        for i=1:t
            
            vel_summand(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * (permute(Vdata(:, tsims(i)-1) - mean(Vdata(:, tsims(i)-1)), [2,1]));
            res_summand(:,:,t,i) = Wmat_pow{t-i+1} * (permute(residuals(:,tsims(i)-1, :) - repmat(mean(residuals(:,tsims(i)-1, :),3), [1,1,length(Trials)]), [1,3,2]));
            res_shuff_summand(:,:,t,i) = Wmat_pow{t-i+1} * (permute(residuals_shuff(:,tsims(i)-1, :) - repmat(mean(residuals_shuff(:,tsims(i)-1, :),3), [1,1,length(Trials)]), [1,3,2]));

        end 
    
        
end

% next execute sums to get quantities which determine covariance
 
cov_vel0 = sum(vel_summand,4);
cov_res0 = sum(res_summand,4);
cov_res_shuff0 = sum(res_shuff_summand,4);
cov_initstate0 = zeros([size(cov_vel0,1), size(cov_vel0,2), length(tsims)]);

rmat_perm = permute(rmat, [2,1,3]);

for t=1:length(tsims)

         cov_initstate0(:,:, t) = Wmat_pow{t + 1} * (rmat_perm(:, Trials, rngtot(tsims(1))-1) - repmat(nanmean(rmat_perm(:, Trials, rngtot(tsims(1))-1),2), [1, length(Trials), 1]));

end

% next concatenate samples and trials into a single vector and compute covariance matrices

Mvel = permute(cov_vel0, [ 3, 2, 1]);
Mvel0 = reshape(Mvel, [size(Mvel,1) * size(Mvel,2), size(Mvel,3)]);
cov_vel = cov(Mvel0);
Mres = permute(cov_res0, [ 3, 2, 1]);
Mres0 = reshape(Mres, [size(Mres,1) * size(Mres,2), size(Mres,3)]);
cov_res = cov(Mres0);
Mres_shuff = permute(cov_res_shuff0, [ 3, 2, 1]);
Mres_shuff0 = reshape(Mres_shuff, [size(Mres,1) * size(Mres,2), size(Mres,3)]);
cov_res_shuff = cov(Mres_shuff0);
Minit = permute(cov_initstate0, [ 3, 2, 1]);
Minit0 = reshape(Minit, [size(Minit,1) * size(Minit,2), size(Minit,3)]);
cov_initstate = cov(Minit0);

cov_vel_init = zeros(size(cov_vel));
cov_vel_res = zeros(size(cov_vel));
cov_res_init = zeros(size(cov_vel));
cov_vel_res_shuff = zeros(size(cov_vel));
cov_res_shuff_init = zeros(size(cov_vel));


for i=1:length(cov_vel)
    for j=1:length(cov_vel)

        C = cov(Mvel0(:,i), Minit0(:,j)) + cov(Minit0(:,i), Mvel0(:,j));
        cov_vel_init(i,j) = C(1,2);
                     
        C = cov(Mvel0(:,i), Mres0(:,j)) + cov(Mres0(:,i), Mvel0(:,j));
        cov_vel_res(i,j) = C(1,2);

        C = cov(Mres0(:,i), Minit0(:,j)) + cov(Minit0(:,i), Mres0(:,j));
        cov_res_init(i,j) = C(1,2);
        
        C = cov(Mvel0(:,i), Mres_shuff0(:,j)) + cov(Mres_shuff0(:,i), Mvel0(:,j));
        cov_vel_res_shuff(i,j) = C(1,2);

        C = cov(Mres_shuff0(:,i), Minit0(:,j)) + cov(Minit0(:,i), Mres_shuff0(:,j));
        cov_res_shuff_init(i,j) = C(1,2);
        
    end
end

cov_tot = cov_vel + cov_res + cov_initstate + cov_vel_res + cov_vel_init + cov_res_init;
D = diag(diag(cov_tot));


% transform to find correlation matrices

corr_vel_tot{rix} = D^(-1/2) * cov_vel * D^(-1/2);
corr_res_tot{rix} = D^(-1/2) * cov_res * D^(-1/2);
corr_res_shuff_tot{rix} = D^(-1/2) * cov_res_shuff * D^(-1/2);
corr_initstate_tot{rix} = D^(-1/2) * cov_initstate * D^(-1/2);

corr_vel_initstate_tot{rix} = D^(-1/2) * cov_vel_init * D^(-1/2);
corr_vel_res_tot{rix} = D^(-1/2) * cov_vel_res * D^(-1/2);
corr_res_initstate_tot{rix} = D^(-1/2) * cov_res_init * D^(-1/2);
corr_vel_res_shuff_tot{rix} = D^(-1/2) * cov_vel_res_shuff * D^(-1/2);
corr_res_shuff_initstate_tot{rix} = D^(-1/2) * cov_res_shuff_init * D^(-1/2);


for i=1:length(CLS)
    for j=1:length(CLS)
        
      Cvel  = corr_vel_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cres  = corr_res_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cres_shuff  = corr_res_shuff_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cinit = corr_initstate_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      
      Cvel_init  = corr_vel_initstate_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cvel_res   = corr_vel_res_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cres_init  = corr_res_initstate_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cvel_res_shuff   = corr_vel_res_shuff_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      Cres_shuff_init  = corr_res_shuff_initstate_tot{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
      
      if i==j
      
            Cvel(1:size(Cvel,1)+1:end) = nan;
            Cres(1:size(Cres,1)+1:end) = nan;
            Cres_shuff(1:size(Cres_shuff,1)+1:end) = nan;
            Cinit(1:size(Cinit,1)+1:end) = nan;

            Cvel_init(1:size(Cvel_init,1)+1:end) = nan;
            Cvel_res(1:size(Cvel_res,1)+1:end) = nan;
            Cres_init(1:size(Cres_init,1)+1:end) = nan;
            Cvel_res_shuff(1:size(Cvel_res_shuff,1)+1:end) = nan;
            Cres_shuff_init(1:size(Cres_shuff_init,1)+1:end) = nan;
            
      end
            
        corr_vel_tot0{i,j} = vertcat(corr_vel_tot0{i,j}, Cvel(:));
        corr_res_tot0{i,j} = vertcat(corr_res_tot0{i,j}, Cres(:));
        corr_res_shuff_tot0{i,j} = vertcat(corr_res_shuff_tot0{i,j}, Cres_shuff(:));
        corr_initstate_tot0{i,j} = vertcat(corr_initstate_tot0{i,j}, Cinit(:));
        
        corr_vel_init_tot0{i,j} = vertcat(corr_vel_init_tot0{i,j}, Cvel_init(:));
        corr_vel_res_tot0{i,j} = vertcat(corr_vel_res_tot0{i,j}, Cvel_res(:));
        corr_res_init_tot0{i,j} = vertcat(corr_res_init_tot0{i,j}, Cres_init(:));
        corr_vel_res_shuff_tot0{i,j} = vertcat(corr_vel_res_shuff_tot0{i,j}, Cvel_res_shuff(:));
        corr_res_shuff_init_tot0{i,j} = vertcat(corr_res_shuff_init_tot0{i,j}, Cres_shuff_init(:));
        
    end
end
        
end

indsx = [1,4,3,2,2,1,2,1,1,3];
indsy = [1,4,3,2,4,2,3,4,3,4];

if strmatch(prepost, 'pre')
% 
% cvel_pre = cellfun(@nanmean, corr_vel_tot0);
% cres_pre = cellfun(@nanmean, corr_res_tot0);
% cinit_pre = cellfun(@nanmean, corr_initstate_tot0);
% 
% cvel_res_pre = cellfun(@nanmean, corr_vel_res_tot0);
% cvel_init_pre = cellfun(@nanmean, corr_vel_init_tot0);
% cres_init_pre = cellfun(@nanmean, corr_res_init_tot0);
% 
% cvel_pre_sem = cellfun(@nanstd, corr_vel_tot0) ./ sqrt(cellfun(@length, corr_vel_tot0));
% cres_pre_sem = cellfun(@nanstd, corr_res_tot0) ./ sqrt(cellfun(@length, corr_res_tot0));
% cinit_pre_sem = cellfun(@nanstd, corr_initstate_tot0) ./ sqrt(cellfun(@length, corr_initstate_tot0));
% 
% cvel_res_pre_sem = cellfun(@nanstd, corr_vel_res_tot0) ./ sqrt(cellfun(@length, corr_vel_res_tot0));
% cvel_init_pre_sem = cellfun(@nanstd, corr_vel_init_tot0) ./ sqrt(cellfun(@length, corr_vel_init_tot0));
% cres_init_pre_sem = cellfun(@nanstd, corr_res_init_tot0) ./ sqrt(cellfun(@length, corr_res_init_tot0));

for i=1:length(indsx)
    
    ctot_pre{i} = corr_vel_tot0{indsx(i), indsy(i)} + corr_res_tot0{indsx(i), indsy(i)} + corr_initstate_tot0{indsx(i), indsy(i)} + corr_vel_res_tot0{indsx(i), indsy(i)} + corr_vel_init_tot0{indsx(i), indsy(i)} + corr_res_init_tot0{indsx(i), indsy(i)};
    ctot_pre_shuffres{i} = corr_vel_tot0{indsx(i), indsy(i)} + corr_res_shuff_tot0{indsx(i), indsy(i)} + corr_initstate_tot0{indsx(i), indsy(i)} + corr_vel_res_shuff_tot0{indsx(i), indsy(i)} + corr_vel_init_tot0{indsx(i), indsy(i)} + corr_res_shuff_init_tot0{indsx(i), indsy(i)};
    
end


elseif strmatch(prepost, 'post')

%     
% cvel_post = cellfun(@nanmean, corr_vel_tot0);
% cres_post = cellfun(@nanmean, corr_res_tot0);
% cinit_post = cellfun(@nanmean, corr_initstate_tot0);
% 
% cvel_res_post = cellfun(@nanmean, corr_vel_res_tot0);
% cvel_init_post = cellfun(@nanmean, corr_vel_init_tot0);
% cres_init_post = cellfun(@nanmean, corr_res_init_tot0);
% 
% cvel_post_sem = cellfun(@nanstd, corr_vel_tot0) ./ sqrt(cellfun(@length, corr_vel_tot0));
% cres_post_sem = cellfun(@nanstd, corr_res_tot0) ./ sqrt(cellfun(@length, corr_res_tot0));
% cinit_post_sem = cellfun(@nanstd, corr_initstate_tot0) ./ sqrt(cellfun(@length, corr_initstate_tot0));
% 
% cvel_res_post_sem = cellfun(@nanstd, corr_vel_res_tot0) ./ sqrt(cellfun(@length, corr_vel_res_tot0));
% cvel_init_post_sem = cellfun(@nanstd, corr_vel_init_tot0) ./ sqrt(cellfun(@length, corr_vel_init_tot0));
% cres_init_post_sem = cellfun(@nanstd, corr_res_init_tot0) ./ sqrt(cellfun(@length, corr_res_init_tot0));

for i=1:length(indsx)
    
    ctot_post{i} = corr_vel_tot0{indsx(i), indsy(i)} + corr_res_tot0{indsx(i), indsy(i)} + corr_initstate_tot0{indsx(i), indsy(i)} + corr_vel_res_tot0{indsx(i), indsy(i)} + corr_vel_init_tot0{indsx(i), indsy(i)} + corr_res_init_tot0{indsx(i), indsy(i)};
    ctot_post_shuffres{i} = corr_vel_tot0{indsx(i), indsy(i)} + corr_res_shuff_tot0{indsx(i), indsy(i)} + corr_initstate_tot0{indsx(i), indsy(i)} + corr_vel_res_shuff_tot0{indsx(i), indsy(i)} + corr_vel_init_tot0{indsx(i), indsy(i)} + corr_res_shuff_init_tot0{indsx(i), indsy(i)};

end

end

end
% 
%  for i=1:length(indsx)
%     cvel_pre0(i) = cvel_pre(indsx(i), indsy(i));
%     cvel_post0(i) = cvel_post(indsx(i), indsy(i));
%     cres_pre0(i) = cres_pre(indsx(i), indsy(i));
%     cres_post0(i) = cres_post(indsx(i), indsy(i));
%     cinit_pre0(i) = cinit_pre(indsx(i), indsy(i));
%     cinit_post0(i) = cinit_post(indsx(i), indsy(i));
%         
%     cvel_res_pre0(i) = cvel_res_pre(indsx(i), indsy(i));
%     cvel_res_post0(i) = cvel_res_post(indsx(i), indsy(i));
%     cvel_init_pre0(i) = cvel_init_pre(indsx(i), indsy(i));
%     cvel_init_post0(i) = cvel_init_post(indsx(i), indsy(i));
%     cres_init_pre0(i) = cres_init_pre(indsx(i), indsy(i));
%     cres_init_post0(i) = cres_init_post(indsx(i), indsy(i));
%     
%     cvel_pre0_sem(i) = cvel_pre_sem(indsx(i), indsy(i));
%     cvel_post0_sem(i) = cvel_post_sem(indsx(i), indsy(i));
%     cres_pre0_sem(i) = cres_pre_sem(indsx(i), indsy(i));
%     cres_post0_sem(i) = cres_post_sem(indsx(i), indsy(i));
%     cinit_pre0_sem(i) = cinit_pre_sem(indsx(i), indsy(i));
%     cinit_post0_sem(i) = cinit_post_sem(indsx(i), indsy(i));
%     
%     cvel_res_pre0_sem(i) = cvel_res_pre_sem(indsx(i), indsy(i));
%     cvel_res_post0_sem(i) = cvel_res_post_sem(indsx(i), indsy(i));
%     cvel_init_pre0_sem(i) = cvel_init_pre_sem(indsx(i), indsy(i));
%     cvel_init_post0_sem(i) = cvel_init_post_sem(indsx(i), indsy(i));
%     cres_init_pre0_sem(i) = cres_init_pre_sem(indsx(i), indsy(i));
%     cres_init_post0_sem(i) = cres_init_post_sem(indsx(i), indsy(i));
%  end
%  
% 
% 
