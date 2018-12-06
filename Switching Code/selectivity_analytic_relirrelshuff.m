
function [SI, poolvar] = selectivity_analytic_relirrelshuff(xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL, RXS)

    dMean_Inputs0 = cell([4,1]);

for prepost = {'pre', 'post'}

for rix = 1:length(RXS)



if strmatch(prepost, 'pre')
% 
%         vel_summandV = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
%         res_summandV = vel_summandV;
%         Inputs_summandV = vel_summandV;
%         vel_summandA = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{2}), length(tsims), length(tsims)]);
%         res_summandA = vel_summandA;
%         Inputs_summandA = vel_summandA;

        vel_summandV = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims)]);
        res_summandV = vel_summandV;
        Inputs_summandV = vel_summandV;
        vel_summandA = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{2}), length(tsims)]);
        res_summandA = vel_summandA;
        Inputs_summandA = vel_summandA;
    
        Wmat = (eye(size(rmat_pre{RXS(rix),1}, 2)) + xpre{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpre{RXS(rix)}(:,end);
        VdataV = vmatpre0{RXS(rix), 1};
        VdataA = vmatpre0{RXS(rix), 2};
        InputsV = xpre{RXS(rix)}(:, (size(xpre{RXS(rix)},1) + 1):(size(xpre{RXS(rix)},1)+tsims(end)));
        InputsA = xpre{RXS(rix)}(:, (size(xpre{RXS(rix)},1) + 1 + 17):(size(xpre{RXS(rix)},1)+tsims(end) + 17));
        residualsV = ResidualV_pre_TOT{RXS(rix)};
        residualsA = ResidualA_pre_TOT{RXS(rix)};
        rmatV = rmat_pre{RXS(rix),1};
        rmatA = rmat_pre{RXS(rix),2};
        TrialsV = Trialspre_all{RXS(rix)}{1};
        TrialsA = Trialspre_all{RXS(rix)}{2};

elseif strmatch(prepost, 'post')

%         vel_summandV = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
%         res_summandV = vel_summandV;
%         Inputs_summandV = vel_summandV;
%         vel_summandA = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{2}), length(tsims), length(tsims)]);
%         res_summandA = vel_summandA;
%         Inputs_summandA = vel_summandA;


        vel_summandV = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims)]);
        res_summandV = vel_summandV;
        Inputs_summandV = vel_summandV;
        vel_summandA = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{2}), length(tsims)]);
        res_summandA = vel_summandA;
        Inputs_summandA = vel_summandA;

        Wmat = (eye(size(rmat_post{RXS(rix),1}, 2)) + xpost{RXS(rix)}(:,1:size(xpre{RXS(rix)},1)));
        Vcoefs = xpost{RXS(rix)}(:,end);
        VdataV = vmatpost0{RXS(rix), 1};
        VdataA = vmatpost0{RXS(rix), 2};
        InputsV = xpost{RXS(rix)}(:, (size(xpost{RXS(rix)},1) + 1):(size(xpost{RXS(rix)},1)+tsims(end)));
        InputsA = xpost{RXS(rix)}(:, (size(xpost{RXS(rix)},1) + 1 + 17):(size(xpost{RXS(rix)},1)+tsims(end) + 17));
        residualsV = ResidualV_post_TOT{RXS(rix)};
        residualsA = ResidualA_post_TOT{RXS(rix)};
        rmatV = rmat_post{RXS(rix),1};
        rmatA = rmat_post{RXS(rix),2};
        TrialsV = Trialspost_all{RXS(rix)}{1};
        TrialsA = Trialspost_all{RXS(rix)}{2};

end
        
    
% compute matrix powers:

Wmat_pow{1} = eye(size(Wmat));

for k = 1:length(tsims)

    Wmat_pow{k + 1} = Wmat_pow{k} * Wmat;  % note the index offset to allow for a zero power

end

    
residualsV_permute = permute(residualsV, [1,3,2]);
residualsA_permute = permute(residualsA, [1,3,2]);

S = zeros(size(Wmat_pow{1},1), size(Wmat_pow{1},1) * length(tsims));
SvelV = zeros(size(Wmat_pow{1},1) * length(tsims), length(TrialsV));
SvelA = zeros(size(Wmat_pow{1},1) * length(tsims), length(TrialsA));
SresV = zeros([size(Wmat_pow{1},1) * length(tsims), size(residualsV_permute,2)]);
SresA = zeros([size(Wmat_pow{1},1) * length(tsims), size(residualsA_permute,2)]);
SInputsV = zeros([size(Wmat_pow{1},1) * length(tsims), length(TrialsV)]);
SInputsA = zeros([size(Wmat_pow{1},1) * length(tsims), length(TrialsA)]);

for k=1:length(tsims)

    S(:, ((k - 1)  * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)))    = Wmat_pow{length(tsims) - k + 1};
    SvelV(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :) = Vcoefs * permute(VdataV(:, tsims(k)-1), [2,1]);
    SvelA(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :) = Vcoefs * permute(VdataA(:, tsims(k)-1), [2,1]);
    SresV(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :,:) = residualsV_permute(:,:,tsims(k)-1);
    SresA(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :,:) = residualsA_permute(:,:,tsims(k)-1);
    SInputsV(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :,:) = repmat(InputsV(:,tsims(k)-1), [1, length(TrialsV)]);
    SInputsA(((k - 1) * size(Wmat_pow{1},1) + 1):(k * size(Wmat_pow{1},1)), :,:) = repmat(InputsA(:,tsims(k)-1), [1, length(TrialsA)]);
    
end

for t=1:length(tsims)              
    
    vel_summandV(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SvelV(1:(t * size(Wmat_pow{1},1)),:);
    vel_summandA(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SvelA(1:(t * size(Wmat_pow{1},1)),:);
    res_summandV(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SresV(1:(t * size(Wmat_pow{1},1)),:,:);
    res_summandA(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SresA(1:(t * size(Wmat_pow{1},1)),:,:);
    Inputs_summandV(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SInputsV(1:(t * size(Wmat_pow{1},1)),:,:);
    Inputs_summandA(:,:,t) = S(:, (end - t * size(Wmat_pow{1},1) + 1):end) * SInputsA(1:(t * size(Wmat_pow{1},1)),:,:);

    
%         for i=1:t
%             
%             vel_summandV(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * permute(VdataV(:, tsims(i)-1), [2,1]);
%             Inputs_summandV(:,:,t,i) = repmat(Wmat_pow{t-i+1} * InputsV(:,tsims(i)-1), [1, length(TrialsV)]);
%             res_summandV(:,:,t,i) = Wmat_pow{t-i+1} * residualsV_permute(:,:,tsims(i)-1);
% 
%             vel_summandA(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * permute(VdataA(:, tsims(i)-1), [2,1]);
%             Inputs_summandA(:,:,t,i) = repmat(Wmat_pow{t-i+1} * InputsA(:,tsims(i)-1), [1, length(TrialsA)]);
%             res_summandA(:,:,t,i) = Wmat_pow{t-i+1} * residualsA_permute(:,:,tsims(i)-1);
% 
%         end    
        
end
            

initstateV0 = zeros([size(CELLLAB_ALL{RXS(rix)},1), size(vel_summandV,2), length(tsims)]);
initstateA0 = zeros([size(CELLLAB_ALL{RXS(rix)},1), size(vel_summandA,2), length(tsims)]);

for t=1:length(tsims)

         initstateV0(:,:,t) = Wmat_pow{t+1} * permute(rmatV(TrialsV, :, rngtot(tsims(1))-1), [2,1,3]);
         initstateA0(:,:,t) = Wmat_pow{t+1} * permute(rmatA(TrialsA, :, rngtot(tsims(1))-1), [2,1,3]);
    
end

nV = length(TrialsV);
nA = length(TrialsA);

% Mean_InputsV{rix} = nanmean(nanmean(sum(Inputs_summandV,4), 3),2);
% Mean_InputsA{rix} = nanmean(nanmean(sum(Inputs_summandA,4), 3), 2);
% Mean_velV{rix} = nanmean(nanmean(sum(vel_summandV,4), 3),2);
% Mean_velA{rix} = nanmean(nanmean(sum(vel_summandA,4), 3),2);
% Mean_resV{rix} = nanmean(nanmean(sum(res_summandV,4), 3),2);
% Mean_resA{rix} = nanmean(nanmean(sum(res_summandA,4), 3),2);
Mean_velV{rix} = nanmean(nanmean(vel_summandV, 3),2);
Mean_velA{rix} = nanmean(nanmean(vel_summandA, 3),2);
Mean_resV{rix} = nanmean(nanmean(res_summandV, 3),2);
Mean_resA{rix} = nanmean(nanmean(res_summandA, 3),2);
Mean_InputsV{rix} = nanmean(nanmean(Inputs_summandV, 3),2);
Mean_InputsA{rix} = nanmean(nanmean(Inputs_summandA, 3), 2);

Mean_InitstateV{rix} = nanmean(nanmean(initstateV0, 3),2) ;
Mean_InitstateA{rix} = nanmean(nanmean(initstateA0, 3), 2);


dMean_Inputs0{rix} = Mean_InputsV{rix} - Mean_InputsA{rix};
dMean_initstate0{rix} = Mean_InitstateV{rix} - Mean_InitstateA{rix};
dMean_vel0{rix} = Mean_velV{rix} - Mean_velA{rix};
dMean_res0{rix} = Mean_resV{rix} - Mean_resA{rix};

poolvar_term1{rix} = nanmean(sum(vel_summandV,4) + sum(res_summandV,4) + initstateV0, 3);
poolvar_term2{rix} = nanmean(sum(vel_summandA,4) + sum(res_summandA,4) + initstateA0,3);
nV0(rix) = nV;
nA0(rix) = nA;

end

for rix = 1:length(RXS)

    poolvar{rix} = (var(poolvar_term1{rix}, [], 2) * (nV0(rix) - 1) +  var(poolvar_term2{rix}, [], 2) * (nA0(rix) - 1)) / (nV0(rix) + nA0(rix) - 2); 

end

dMean_Inputs0_tot = vertcat(dMean_Inputs0{:});
dMean_initstate0_tot = vertcat(dMean_initstate0{:});
dMean_vel0_tot = vertcat(dMean_vel0{:});
dMean_res0_tot = vertcat(dMean_res0{:});

poolvar_tot = vertcat(poolvar{:});

if strmatch(prepost, 'pre')
    
  dMean_Inputs0_tot_pre = dMean_Inputs0_tot;
  dMean_initstate0_tot_pre = dMean_initstate0_tot;
  dMean_vel0_tot_pre = dMean_vel0_tot;
  dMean_res0_tot_pre = dMean_res0_tot;

  poolvar_tot_pre = poolvar_tot;
  
 
elseif strmatch(prepost, 'post')
    
  dMean_Inputs0_tot_post = dMean_Inputs0_tot;
  dMean_initstate0_tot_post = dMean_initstate0_tot;
  dMean_vel0_tot_post = dMean_vel0_tot;    
  dMean_res0_tot_post = dMean_res0_tot;
  
  poolvar_tot_post = poolvar_tot;

end

end
       
    poolvar_pre = poolvar_tot_pre;
    poolvar_post = poolvar_tot_post;
    Xinputs = [dMean_Inputs0_tot_pre, dMean_Inputs0_tot_post];
    Xinitstate = [dMean_initstate0_tot_pre, dMean_initstate0_tot_post];
    Xvel = [dMean_vel0_tot_pre, dMean_vel0_tot_post];
    Xres = [dMean_res0_tot_pre, dMean_res0_tot_post];
    
    
          
        SIinputs(:,1) = Xinputs(:,1)  ./ sqrt(poolvar_pre);
        SIinputs(:,2) = Xinputs(:,2)  ./ sqrt(poolvar_post);

        SIinitstate(:,1) = Xinitstate(:,1)  ./ sqrt(poolvar_pre);
        SIinitstate(:,2) = Xinitstate(:,2)  ./ sqrt(poolvar_post);

        SIvel(:,1) = Xvel(:,1)  ./ sqrt(poolvar_pre);
        SIvel(:,2) = Xvel(:,2)  ./ sqrt(poolvar_post);

        
        SIres(:,1) = Xres(:,1)  ./ sqrt(poolvar_pre);
        SIres(:,2) = Xres(:,2)  ./ sqrt(poolvar_post);

        SI = struct;
        SI.SIinputs = SIinputs;
        SI.SIinitstate = SIinitstate;
        SI.SIvel = SIvel;
        SI.SIres = SIres;
        
        poolvar = struct;
        poolvar.poolvar_pre = poolvar_pre;
        poolvar.poolvar_post = poolvar_post;

        
 
       