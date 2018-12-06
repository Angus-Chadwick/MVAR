function [xpre, xpost, Measurement_pre, Measurement_post, Prediction_pre, Prediction_post] = fitLDS_nospeed(fixvars, changeW, changeI, drmat0pre, drmat0post, rmat0pre, rmat0post, rngtot, Ntrialspre, Ntrialspost, CELLLAB_ALL)



    rpretot = horzcat(rmat0pre{:});
    rposttot = horzcat(rmat0post{:});
    stimblockspre  = [repmat(eye(length(rngtot)),   [1, Ntrialspre{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspre{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspre{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspre{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspre{3}, 1, 1])];
    stimblockspost = [repmat(eye(length(rngtot)),   [1, Ntrialspost{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspost{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{3}, 1, 1]); ...
                                              repmat(zeros(length(rngtot)), [1, Ntrialspost{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialspost{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialspost{3}, 1, 1])];



drmat0tot = [horzcat(drmat0pre{:}), horzcat(drmat0post{:})];

if strmatch(fixvars, 'Inputs')

elseif strmatch(fixvars, 'Weights')
    
elseif strmatch(fixvars, 'Mixed') % allow only some inputs to change
       
    CLS = [1,3,4,5];
           
    for j=1:4  % for predictions of cell type j
        
        changeWlin = CLS(find(changeW(j,:)));
        
        rpostWfixed  = rposttot;
        rpostWfixed(ismember(CELLLAB_ALL, changeWlin), :) = 0;
        rpostWchange = rposttot;
        rpostWchange(~ismember(CELLLAB_ALL, changeWlin), :) = 0;
              
        Lambdatot_partialW_ChangeI = [rpretot, rpostWfixed;  zeros(size(rpretot)), rpostWchange; stimblockspre, zeros(size(stimblockspost)); zeros(size(stimblockspre)), stimblockspost];  % only allow weights to change                                     
        Lambdatot_partialW_fixI    = [rpretot, rpostWfixed;  zeros(size(rpretot)), rpostWchange; stimblockspre, stimblockspost];
             

        
            if ~changeI(j)
                
                                     
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_fixI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xpre(CELLLAB_ALL == CLS(j),:)  = xtot{j}(:, [1:length(CELLLAB_ALL), (2*length(CELLLAB_ALL) + 1):end]);
                        xpost(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 1):end]);        

                        Prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_fixI;
                     
                    
             elseif changeI(j) 
        
                 
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_ChangeI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xpre(CELLLAB_ALL == CLS(j),:)  = xtot{j}(:, [1:length(CELLLAB_ALL),  (2*length(CELLLAB_ALL) + 1):(2*length(CELLLAB_ALL) + 51)]);
                        xpost(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 51 + 1):end]);        

                        Prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_ChangeI;

           end
    
    end

Measurement_pre = horzcat(drmat0pre{:});
Measurement_post = horzcat(drmat0post{:});
Prediction_pre = Prediction_tot(:, 1:size(Measurement_pre,2));
Prediction_post = Prediction_tot(:, (size(Measurement_pre,2) + 1):end);

    
end