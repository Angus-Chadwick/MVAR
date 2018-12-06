
function [SI, poolvar] = selectivity_analytic_prepostshuff_efficient(xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL)

RXS = [1,3,4,5,6,8,9,10];

    dMean_Inputs0 = cell([4,1]);

for prepost = {'pre', 'post'}

for rix = 1:length(RXS)



if strmatch(prepost, 'pre')

    
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
     
rsimV = zeros([size(rmatV,1), size(rmatV,2), length(tsims) + 1]);
rsimA = zeros([size(rmatA,1), size(rmatA,2), length(tsims) + 1]);

rsimA(:,:,1) = rmatV(:,:,tsims(1) - 1);  % set initial state for each cell and each trial equal to data
rsimA(:,:,1) = rmatV(:,:,tsims(1) - 1);  % set initial state for each cell and each trial equal to data
   
for t=1:length(tsims)
            
    rsimA(:,:,t) = Wmat * rsimA(:,:,t-1) + InputsA(:, tsims(t) - 1) + Vcoefs * VdataA(:,tsims(t) - 1)' + permute(residualsA(:,tsims(t) - 1, :), [1,3,2]); 
    rsimV(:,:,t) = Wmat * rsimV(:,:,t-1) + InputsV(:, tsims(t) - 1) + Vcoefs * VdataV(:,tsims(t) - 1)' + permute(residualsV(:,tsims(t) - 1, :), [1,3,2]); 
            
end
            

nV = length(TrialsV);
nA = length(TrialsA);

Mean_InputsV{rix} = nanmean(nanmean(sum(Inputs_summandV,4), 3),2);
Mean_InputsA{rix} = nanmean(nanmean(sum(Inputs_summandA,4), 3), 2);
Mean_InitstateV{rix} = nanmean(nanmean(initstateV0, 3),2) ;
Mean_InitstateA{rix} = nanmean(nanmean(initstateA0, 3), 2);
Mean_velV{rix} = nanmean(nanmean(sum(vel_summandV,4), 3),2);
Mean_velA{rix} = nanmean(nanmean(sum(vel_summandA,4), 3),2);
Mean_resV{rix} = nanmean(nanmean(sum(res_summandV,4), 3),2);
Mean_resA{rix} = nanmean(nanmean(sum(res_summandA,4), 3),2);

dMean_Inputs0{rix} = Mean_InputsV{rix} - Mean_InputsA{rix};
dMean_initstate0{rix} = Mean_InitstateV{rix} - Mean_InitstateA{rix};
dMean_vel0{rix} = Mean_velV{rix} - Mean_velA{rix};
dMean_res0{rix} = Mean_resV{rix} - Mean_resA{rix};

poolvar{rix} = (var(nanmean(sum(vel_summandV,4), 3) + nanmean(sum(res_summandV,4), 3) + nanmean(initstateV0,3), [], 2) * (nV - 1) +  var(nanmean(sum(vel_summandA,4), 3) + nanmean(sum(res_summandA,4), 3) + nanmean(initstateA0,3), [], 2) * (nA - 1)) / (nV + nA - 2); 


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

        
 
       