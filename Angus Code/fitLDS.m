function [xpre, xpost, Measurement_pre, Measurement_post, Prediction_pre, Prediction_post] = fitLDS(fixvars, changeW, changeI, drmat0pre, drmat0post, rmat0pre, rmat0post, rngtot, Ntrialspre, Ntrialspost, vmatSTM_pre0, vmatSTM_post0, CELLLAB_ALL)



    rpretot = horzcat(rmat0pre{:});
    rposttot = horzcat(rmat0post{:});
    stimblockspre  = [repmat(eye(length(rngtot)),   [1, Ntrialspre{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspre{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspre{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspre{3}, 1, 1])];
    stimblockspost = [repmat(eye(length(rngtot)),   [1, Ntrialspost{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspost{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspost{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspost{3}, 1, 1])];
    vpre =  vertcat(vmatSTM_pre0{:}).';              
    vpost =  vertcat(vmatSTM_post0{:}).';                                  
   
    Lambdatot_fixI = [rpretot, zeros(size(rposttot)); zeros(size(rpretot)), rposttot; stimblockspre, stimblockspost; vpre, vpost];  % only allow weights to change                                     
    Lambdatot_fixW = [rpretot, rposttot; stimblockspre, zeros(size(stimblockspost)); zeros(size(stimblockspre)), stimblockspost; vpre, vpost];  % only allow weights to change                                     
    Lambdatot_full = [rpretot, zeros(size(rposttot)); zeros(size(rpretot)), rposttot; stimblockspre, zeros(size(stimblockspost)); zeros(size(stimblockspre)), stimblockspost; vpre, vpost];  % only allow weights to change                                     
    


drmat0tot = [horzcat(drmat0pre{:}), horzcat(drmat0post{:})];

if strmatch(fixvars, 'Inputs')

    xtot = drmat0tot / Lambdatot_fixI;
    xpre = xtot(1:end, [1:size(xtot, 1), (2*size(xtot, 1) + 1):end]);
    xpost = xtot(1:end, [(size(xtot, 1)+1):(2*size(xtot, 1)), (2*size(xtot, 1) + 1):end]);

    Prediction_tot = xtot * Lambdatot_fixI;
    
elseif strmatch(fixvars, 'Weights')

    xtot = drmat0tot / Lambdatot_fixW;
    xpre  = xtot(1:end, [1:(size(xtot, 1) + 3 * 17), end]);
    xpost = xtot(1:end, [1:(size(xtot, 1)), (size(xtot, 1) + 3 * 17 + 1):end]);
    
    Prediction_tot = xtot * Lambdatot_fixW;
    
elseif strmatch(fixvars, 'Mixed') % allow only some inputs to change
       
    CLS = [1,3,4,5];
           
    for j=1:4  % for predictions of cell type j
        
        changeWlin = CLS(find(changeW(j,:)));
        
        rpostWfixed  = rposttot;
        rpostWfixed(ismember(CELLLAB_ALL, changeWlin), :) = 0;
        rpostWchange = rposttot;
        rpostWchange(~ismember(CELLLAB_ALL, changeWlin), :) = 0;
              
        Lambdatot_partialW_ChangeI = [rpretot, rpostWfixed;  zeros(size(rpretot)), rpostWchange; stimblockspre, zeros(size(stimblockspost)); zeros(size(stimblockspre)), stimblockspost; vpre, vpost];  % only allow weights to change                                     
        Lambdatot_partialW_fixI    = [rpretot, rpostWfixed;  zeros(size(rpretot)), rpostWchange; stimblockspre, stimblockspost; vpre, vpost];
             

        
            if ~changeI(j)
                
                                     
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_fixI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xpre(CELLLAB_ALL == CLS(j),:)  = xtot{j}(:, [1:length(CELLLAB_ALL), (2*length(CELLLAB_ALL) + 1):end]);
                        xpost(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 1):end]);        

                        Prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_fixI;
                     
                    
             elseif changeI(j) 
        
                 
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_ChangeI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xpre(CELLLAB_ALL == CLS(j),:)  = xtot{j}(:, [1:length(CELLLAB_ALL),  (2*length(CELLLAB_ALL) + 1):(2*length(CELLLAB_ALL) + 51), end]);
                        xpost(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 51 + 1):end]);        

                        Prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_ChangeI;

           end
    
    end

Measurement_pre = horzcat(drmat0pre{:});
Measurement_post = horzcat(drmat0post{:});
Prediction_pre = Prediction_tot(:, 1:size(Measurement_pre,2));
Prediction_post = Prediction_tot(:, (size(Measurement_pre,2) + 1):end);

    
end