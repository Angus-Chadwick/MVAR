
function [SI, poolvar, TotalV] = selectivity_analytic_switching(manipulation, var_method, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL, RXS)

%%Note: This code has been updated from the learning code to include the
%%RXS vector as an argument of the function. It is otherwise identical, and
%%can handle both learning and switching data.

for removeweights = 0:1

    dMean_Inputs0 = cell([4,1]);

for prepost = {'pre', 'post'}

for rix = 1:length(RXS)



if strmatch(prepost, 'pre')

        vel_summandV = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summandV = vel_summandV;
        Inputs_summandV = vel_summandV;
        vel_summandA = zeros([size(xpre{RXS(rix)},1), length(Trialspre_all{RXS(rix)}{2}), length(tsims), length(tsims)]);
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

        vel_summandV = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{1}), length(tsims), length(tsims)]);
        res_summandV = vel_summandV;
        Inputs_summandV = vel_summandV;
        vel_summandA = zeros([size(xpost{RXS(rix)},1), length(Trialspost_all{RXS(rix)}{2}), length(tsims), length(tsims)]);
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
        
    
if removeweights
    
    if strmatch(manipulation, 'removeall')
        
        Wmat = diag(diag(Wmat));
        
    elseif strmatch(manipulation, 'removeW1W2')
    
        Wmat_diag = diag(diag(Wmat));
        Wmat(CELLLAB_ALL{RXS(rix)} == W1, CELLLAB_ALL{RXS(rix)} == W2) = 0;
        Wmat = Wmat - diag(diag(Wmat)) + Wmat_diag;
    
    end
end


% compute matrix powers:

for k = 0:length(tsims)
    
    Wmat_pow{k + 1} = Wmat^k;  % note the index offset to allow for a zero power
    
end

for t=1:length(tsims)
   
                
        for i=1:t
            
            vel_summandV(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * permute(VdataV(:, tsims(i)-1), [2,1]);
            Inputs_summandV(:,:,t,i) = repmat(Wmat_pow{t-i+1} * InputsV(:,tsims(i)-1), [1, length(TrialsV)]);
            res_summandV(:,:,t,i) = Wmat_pow{t-i+1} * permute(residualsV(:,tsims(i)-1, :), [1,3,2]);

            vel_summandA(:,:,t,i) = Wmat_pow{t-i+1} * Vcoefs * permute(VdataA(:, tsims(i)-1), [2,1]);
            Inputs_summandA(:,:,t,i) = repmat(Wmat_pow{t-i+1} * InputsA(:,tsims(i)-1), [1, length(TrialsA)]);
            res_summandA(:,:,t,i) = Wmat_pow{t-i+1} * permute(residualsA(:,tsims(i)-1, :), [1,3,2]);

        end 
    
        
end

fullresiduals = 1;  % partially conditional version - condition the effect of residuals over the residual inputs to all other cells but not over time

if fullresiduals 
    
    delta_res_summandV = zeros(size(res_summandV));
    delta_res_summandA = zeros(size(res_summandA));

    
    for t=1:length(tsims)                 
        for i=1:t

            delta_res_summandV(:,:,t,i) = diag(diag(Wmat_pow{t-i+1})) * permute(residualsV(:,tsims(i)-1, :), [1,3,2]);
            delta_res_summandA(:,:,t,i) = diag(diag(Wmat_pow{t-i+1})) * permute(residualsA(:,tsims(i)-1, :), [1,3,2]);
            
        end
    end
        
end
                 

initstateV0 = zeros([size(CELLLAB_ALL{RXS(rix)},1), size(vel_summandV,2), length(tsims)]);
initstateA0 = zeros([size(CELLLAB_ALL{RXS(rix)},1), size(vel_summandA,2), length(tsims)]);

for t=1:length(tsims)

         initstateV0(:,:,t) = Wmat_pow{t+1} * permute(rmatV(TrialsV, :, rngtot(tsims(1))-1), [2,1,3]);
         initstateA0(:,:,t) = Wmat_pow{t+1} * permute(rmatA(TrialsA, :, rngtot(tsims(1))-1), [2,1,3]);
    
end

nV = length(TrialsV);
nA = length(TrialsA);

if removeweights == 0

    TotalV{rix} = sum(Inputs_summandV,4) + initstateV0 + sum(vel_summandV,4) + sum(res_summandV,4);

end

Mean_InputsV{rix} = nanmean(nanmean(sum(Inputs_summandV,4), 3),2);
Mean_InputsA{rix} = nanmean(nanmean(sum(Inputs_summandA,4), 3), 2);
Mean_InitstateV{rix} = nanmean(nanmean(initstateV0, 3),2) ;
Mean_InitstateA{rix} = nanmean(nanmean(initstateA0, 3), 2);
Mean_velV{rix} = nanmean(nanmean(sum(vel_summandV,4), 3),2);
Mean_velA{rix} = nanmean(nanmean(sum(vel_summandA,4), 3),2);
Mean_resV{rix} = nanmean(nanmean(sum(res_summandV,4), 3),2);
Mean_resA{rix} = nanmean(nanmean(sum(res_summandA,4), 3),2);
% 
% % fully conditional version - condition effect of residuals over the state of all cells at all times
% 
% Mean_delta_resV{rix} = nanmean(nanmean(residualsV(:,tsims - 1, :), 3),2);
% Mean_delta_resA{rix} = nanmean(nanmean(residualsA(:,tsims - 1, :), 3),2);

% partially conditional version - condition over residuals of all other cells but not time

Mean_delta_resV{rix} = nanmean(nanmean(sum(delta_res_summandV,4), 3),2);
Mean_delta_resA{rix} = nanmean(nanmean(sum(delta_res_summandA,4), 3),2);

dMean_Inputs0{rix} = Mean_InputsV{rix} - Mean_InputsA{rix};
dMean_initstate0{rix} = Mean_InitstateV{rix} - Mean_InitstateA{rix};
dMean_vel0{rix} = Mean_velV{rix} - Mean_velA{rix};
dMean_res0{rix} = Mean_resV{rix} - Mean_resA{rix};
dMean_delta_res0{rix} = Mean_delta_resV{rix} - Mean_delta_resA{rix};

poolvar{rix} = (var(nanmean(sum(vel_summandV,4), 3) + nanmean(sum(res_summandV,4), 3) + nanmean(initstateV0,3), [], 2) * (nV - 1) +  var(nanmean(sum(vel_summandA,4), 3) + nanmean(sum(res_summandA,4), 3) + nanmean(initstateA0,3), [], 2) * (nA - 1)) / (nV + nA - 2); 

SI_novar{rix}      = (dMean_Inputs0{rix} + dMean_initstate0{rix} + dMean_vel0{rix})                   ./ (abs(Mean_InputsV{rix} + Mean_InitstateV{rix} + Mean_velV{rix})                  + abs(Mean_InputsA{rix} + Mean_InitstateA{rix} + Mean_velA{rix}));
SI_novar_full{rix} = (dMean_Inputs0{rix} + dMean_initstate0{rix} + dMean_vel0{rix} + dMean_res0{rix}) ./ (abs(Mean_InputsV{rix} + Mean_InitstateV{rix} + Mean_velV{rix} + Mean_resV{rix}) + abs(Mean_InputsA{rix} + Mean_InitstateA{rix} + Mean_velA{rix} + Mean_resA{rix}));

end

dMean_Inputs0_tot = vertcat(dMean_Inputs0{:});
dMean_initstate0_tot = vertcat(dMean_initstate0{:});
dMean_vel0_tot = vertcat(dMean_vel0{:});
dMean_res0_tot = vertcat(dMean_res0{:});
dMean_delta_res0_tot = vertcat(dMean_delta_res0{:});

poolvar_tot = vertcat(poolvar{:});

SI_novar_tot = vertcat(SI_novar{:});
SI_novar_full_tot = vertcat(SI_novar_full{:});


if strmatch(prepost, 'pre')
    
  dMean_Inputs0_tot_pre = dMean_Inputs0_tot;
  dMean_initstate0_tot_pre = dMean_initstate0_tot;
  dMean_vel0_tot_pre = dMean_vel0_tot;
  dMean_res0_tot_pre = dMean_res0_tot;
  dMean_delta_res0_tot_pre = dMean_delta_res0_tot;

  poolvar_tot_pre = poolvar_tot;

  SI_novar_tot_pre = SI_novar_tot;
  SI_novar_full_tot_pre = SI_novar_full_tot;

  
%  PSTH_tot_pre = PSTH;
  
elseif strmatch(prepost, 'post')
    
  dMean_Inputs0_tot_post = dMean_Inputs0_tot;
  dMean_initstate0_tot_post = dMean_initstate0_tot;
  dMean_vel0_tot_post = dMean_vel0_tot;    
  dMean_res0_tot_post = dMean_res0_tot;
  dMean_delta_res0_tot_post = dMean_delta_res0_tot;
  
  poolvar_tot_post = poolvar_tot;

  SI_novar_tot_post = SI_novar_tot;
  SI_novar_full_tot_post = SI_novar_full_tot;

%  PSTH_tot_post = PSTH;
  
end

end


if removeweights == 0;
        
    poolvar_pre = poolvar_tot_pre;
    poolvar_post = poolvar_tot_post;
    poolvar_ref = (poolvar_tot_pre + poolvar_tot_post) / 2;
    Xinputs = [dMean_Inputs0_tot_pre, dMean_Inputs0_tot_post];
    Xinitstate = [dMean_initstate0_tot_pre, dMean_initstate0_tot_post];
    Xvel = [dMean_vel0_tot_pre, dMean_vel0_tot_post];
    Xres = [dMean_res0_tot_pre, dMean_res0_tot_post];
    Xdeltares = [dMean_delta_res0_tot_pre, dMean_delta_res0_tot_post];
    
    SI_novar0 = [SI_novar_tot_pre, SI_novar_tot_post];
    SI_novar_full0 = [SI_novar_full_tot_pre, SI_novar_full_tot_post];

%    PSTH_weights = [PSTH_tot_pre, PSTH_tot_post];
    
elseif removeweights == 1;
    
    poolvar_pre_rw = poolvar_tot_pre;
    poolvar_post_rw = poolvar_tot_post;
    poolvar_ref_rw = (poolvar_tot_pre + poolvar_tot_post) / 2;
    Xinputs_rw = [dMean_Inputs0_tot_pre, dMean_Inputs0_tot_post];
    Xinitstate_rw = [dMean_initstate0_tot_pre, dMean_initstate0_tot_post];
    Xvel_rw = [dMean_vel0_tot_pre, dMean_vel0_tot_post];
    Xres_rw = [dMean_res0_tot_pre, dMean_res0_tot_post];
    Xdeltares_rw = [dMean_delta_res0_tot_pre, dMean_delta_res0_tot_post];
    
    SI_novar0_rw = [SI_novar_tot_pre, SI_novar_tot_post];
    SI_novar_full0_rw = [SI_novar_full_tot_pre, SI_novar_full_tot_post];
    
 %   PSTH_rw = [PSTH_tot_pre, PSTH_tot_post];
    
end

end

poolvar_ref_tot = (poolvar_ref + poolvar_ref_rw) / 2;

% clear PSTH
% PSTH = struct;
% PSTH.PSTH = PSTH_weights;
% PSTH.PSTH_rw = PSTH_rw;

for i=1:2
    
    if strmatch(var_method, 'fixedvar')
    
        SIinputs(:,i) = Xinputs(:,i) ./ sqrt(poolvar_ref_tot);    
        SIinputs_rw(:,i) = Xinputs_rw(:,i) ./ sqrt(poolvar_ref_tot);
        
        SI = struct;
        SI.SIinputs = SIinputs;
        SI.SIinputs_rw = SIinputs_rw;
        
        poolvar = struct;
        poolvar.poolvar_pre = poolvar_pre;
        poolvar.poolvar_post = poolvar_post;
        poolvar.poolvar_pre_rw = poolvar_pre_rw;
        poolvar.poolvar_post_rw = poolvar_post_rw;
        
    elseif strmatch(var_method, 'fixedvar_rw')
        
        SIinputs(:,1) = Xinputs(:,1)  ./ sqrt(poolvar_pre);
        SIinputs(:,2) = Xinputs(:,2)  ./ sqrt(poolvar_post);
        SIinputs_rw(:,1) = Xinputs_rw(:,1)  ./ sqrt(poolvar_pre);
        SIinputs_rw(:,2) = Xinputs_rw(:,2)  ./ sqrt(poolvar_post);
        
        SIinitstate(:,1) = Xinitstate(:,1)  ./ sqrt(poolvar_pre);
        SIinitstate(:,2) = Xinitstate(:,2)  ./ sqrt(poolvar_post);
        SIinitstate_rw(:,1) = Xinitstate_rw(:,1)  ./  sqrt(poolvar_pre);
        SIinitstate_rw(:,2) = Xinitstate_rw(:,2)  ./ sqrt(poolvar_post);

        SIvel(:,1) = Xvel(:,1)  ./ sqrt(poolvar_pre);
        SIvel(:,2) = Xvel(:,2)  ./ sqrt(poolvar_post);
        SIvel_rw(:,1) = Xvel_rw(:,1)  ./ sqrt(poolvar_pre);
        SIvel_rw(:,2) = Xvel_rw(:,2)  ./ sqrt(poolvar_post);
        
        SIres(:,1) = Xres(:,1)  ./ sqrt(poolvar_pre);
        SIres(:,2) = Xres(:,2)  ./  sqrt(poolvar_post);
        SIres_rw(:,1) = Xres_rw(:,1)  ./ sqrt(poolvar_pre);
        SIres_rw(:,2) = Xres_rw(:,2)  ./ sqrt(poolvar_post);
        
        SI = struct;
        SI.SIinputs = SIinputs;
        SI.SIinputs_rw = SIinputs_rw;
        SI.SIinitstate = SIinitstate;
        SI.SIinitstate_rw = SIinitstate_rw;
        SI.SIvel = SIvel;
        SI.SIvel_rw = SIvel_rw;
        SI.SIres = SIres;
        SI.SIres_rw = SIres_rw;

        poolvar = struct;
        poolvar.poolvar_pre = poolvar_pre;
        poolvar.poolvar_post = poolvar_post;
        poolvar.poolvar_pre_rw = poolvar_pre_rw;
        poolvar.poolvar_post_rw = poolvar_post_rw;
        
    elseif strmatch(var_method, 'fullvar')
        
        SIinputs(:,1) = Xinputs(:,1)  ./ sqrt(poolvar_pre);
        SIinputs(:,2) = Xinputs(:,2)  ./ sqrt(poolvar_post);
        SIinputs_rw(:,1) = Xinputs_rw(:,1)  ./ sqrt(poolvar_pre_rw);
        SIinputs_rw(:,2) = Xinputs_rw(:,2)  ./ sqrt(poolvar_post_rw);
        
        SIinitstate(:,1) = Xinitstate(:,1)  ./ sqrt(poolvar_pre);
        SIinitstate(:,2) = Xinitstate(:,2)  ./ sqrt(poolvar_post);
        SIinitstate_rw(:,1) = Xinitstate_rw(:,1)  ./ sqrt(poolvar_pre_rw);
        SIinitstate_rw(:,2) = Xinitstate_rw(:,2)  ./ sqrt(poolvar_post_rw);

        SIvel(:,1) = Xvel(:,1)  ./ sqrt(poolvar_pre);
        SIvel(:,2) = Xvel(:,2)  ./ sqrt(poolvar_post);
        SIvel_rw(:,1) = Xvel_rw(:,1)  ./ sqrt(poolvar_pre_rw);
        SIvel_rw(:,2) = Xvel_rw(:,2)  ./ sqrt(poolvar_post_rw);
        
        SIres(:,1) = Xres(:,1)  ./ sqrt(poolvar_pre);
        SIres(:,2) = Xres(:,2)  ./ sqrt(poolvar_post);
        SIres_rw(:,1) = Xres_rw(:,1)  ./ sqrt(poolvar_pre_rw);
        SIres_rw(:,2) = Xres_rw(:,2)  ./ sqrt(poolvar_post_rw);
        
        SIdeltares(:,1) = Xdeltares(:,1)  ./ sqrt(poolvar_pre);
        SIdeltares(:,2) = Xdeltares(:,2)  ./ sqrt(poolvar_post);
        SIdeltares_rw(:,1) = Xdeltares_rw(:,1)  ./ sqrt(poolvar_pre_rw);
        SIdeltares_rw(:,2) = Xdeltares_rw(:,2)  ./ sqrt(poolvar_post_rw);
        
        SI = struct;
        SI.SIinputs = SIinputs;
        SI.SIinputs_rw = SIinputs_rw;
        SI.SIinitstate = SIinitstate;
        SI.SIinitstate_rw = SIinitstate_rw;
        SI.SIvel = SIvel;
        SI.SIvel_rw = SIvel_rw;
        SI.SIres = SIres;
        SI.SIres_rw = SIres_rw;
        SI.SIdeltares = SIdeltares;
        SI.SIdeltares_rw = SIdeltares_rw;
        
        poolvar = struct;
        poolvar.poolvar_pre = poolvar_pre;
        poolvar.poolvar_post = poolvar_post;
        poolvar.poolvar_pre_rw = poolvar_pre_rw;
        poolvar.poolvar_post_rw = poolvar_post_rw;
        
     elseif strmatch(var_method, 'novar')
       
        SI = struct;    
        SI.SItot = SI_novar0;       
        SI.SItot_full = SI_novar_full0;       
        SI.SItot_rw = SI_novar0_rw;       
        SI.SItot_full_rw = SI_novar_full0_rw;      
        
        poolvar = struct;
        poolvar.poolvar_pre = poolvar_pre;
        poolvar.poolvar_post = poolvar_post;
        poolvar.poolvar_pre_rw = poolvar_pre_rw;
        poolvar.poolvar_post_rw = poolvar_post_rw;      
         
    end
        
end
