function [xrel, xirrel, Measurement_rel, Measurement_irrel, prediction_rel, prediction_irrel] = fitLDS_Switching_nospeed(fixvars, changeW, changeI, drmat0rel, drmat0irrel, rmat0rel, rmat0irrel, rngtot, Ntrialsrel, Ntrialsirrel, CELLLAB_ALL)

Ncond = 3;

    rreltot = horzcat(rmat0rel{1:Ncond});
    rirreltot = horzcat(rmat0irrel{1:Ncond});

    
%     stimblocksrel  = [repmat(eye(length(rngtot)),   [1, Ntrialsrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsrel{2}, 1, 1]); ...
%                       repmat(zeros(length(rngtot)), [1, Ntrialsrel{1}, 1, 1]), repmat(eye(length(rngtot)),   [1, Ntrialsrel{2}, 1, 1])];
%                                           
%     stimblocksirrel = [repmat(eye(length(rngtot)),   [1, Ntrialsirrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsirrel{2}, 1, 1]); ...
%                        repmat(zeros(length(rngtot)), [1, Ntrialsirrel{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialsirrel{2}, 1, 1])];
  
                   
    stimblocksrel  = [repmat(eye(length(rngtot)),   [1, Ntrialsrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsrel{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsrel{3}, 1, 1]); ...
                      repmat(zeros(length(rngtot)), [1, Ntrialsrel{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialsrel{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsrel{3}, 1, 1]); ...
                      repmat(zeros(length(rngtot)), [1, Ntrialsrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsrel{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialsrel{3}, 1, 1])];
    stimblocksirrel = [repmat(eye(length(rngtot)),   [1, Ntrialsirrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsirrel{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsirrel{3}, 1, 1]); ...
                       repmat(zeros(length(rngtot)), [1, Ntrialsirrel{1}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialsirrel{2}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsirrel{3}, 1, 1]); ...
                       repmat(zeros(length(rngtot)), [1, Ntrialsirrel{1}, 1, 1]), repmat(zeros(length(rngtot)), [1, Ntrialsirrel{2}, 1, 1]), repmat(eye(length(rngtot)), [1, Ntrialsirrel{3}, 1, 1])];
                   
                   stimblocksrel   = stimblocksrel(1:(17*Ncond), 1:(sum(horzcat(Ntrialsrel{1:Ncond}))*17));
                   stimblocksirrel = stimblocksirrel(1:(17*Ncond), 1:(sum(horzcat(Ntrialsirrel{1:Ncond}))*17));

                                 
   
    Lambdatot_fixI = [rreltot, zeros(size(rirreltot)); zeros(size(rreltot)), rirreltot; stimblocksrel, stimblocksirrel];  % only allow weights to change                                     
    Lambdatot_fixW = [rreltot, rirreltot; stimblocksrel, zeros(size(stimblocksirrel)); zeros(size(stimblocksrel)), stimblocksirrel];  % only allow weights to change                                     
    


drmat0tot = [horzcat(drmat0rel{1:Ncond}), horzcat(drmat0irrel{1:Ncond})];

if strmatch(fixvars, 'Inputs')

    xtot = drmat0tot / Lambdatot_fixI;
    xrel = xtot(1:end, [1:size(xtot, 1), (2*size(xtot, 1) + 1):end]);
    xirrel = xtot(1:end, [(size(xtot, 1)+1):(2*size(xtot, 1)), (2*size(xtot, 1) + 1):end]);

    prediction_tot = xtot * Lambdatot_fixI;
    
elseif strmatch(fixvars, 'Weights')

    xtot = drmat0tot / Lambdatot_fixW;
    xrel  = xtot(1:end, [1:(size(xtot, 1) + 3 * 17), end]);
    xirrel = xtot(1:end, [1:(size(xtot, 1)), (size(xtot, 1) + 3 * 17 + 1):end]);
    
    prediction_tot = xtot * Lambdatot_fixW;
    
elseif strmatch(fixvars, 'Mixed') % allow only some inputs to change
       
    CLS = [1,3,4,5];
           
    for j=1:4  % for predictions of cell type j
        
        changeWlin = CLS(find(changeW(j,:)));
        
        rirrelWfixed  = rirreltot;
        rirrelWfixed(ismember(CELLLAB_ALL, changeWlin), :) = 0;
        rirrelWchange = rirreltot;
        rirrelWchange(~ismember(CELLLAB_ALL, changeWlin), :) = 0;
              
        Lambdatot_partialW_ChangeI = [rreltot, rirrelWfixed;  zeros(size(rreltot)), rirrelWchange; stimblocksrel, zeros(size(stimblocksirrel)); zeros(size(stimblocksrel)), stimblocksirrel];  % only allow weights to change                                     
        Lambdatot_partialW_fixI    = [rreltot, rirrelWfixed;  zeros(size(rreltot)), rirrelWchange; stimblocksrel, stimblocksirrel];
             

        
            if ~changeI(j)
                
                                     
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_fixI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xrel(CELLLAB_ALL == CLS(j),:)   = xtot{j}(:, [1:length(CELLLAB_ALL), (2*length(CELLLAB_ALL) + 1):end]);
                        xirrel(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 1):end]);        

                        prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_fixI;
                     
                    
             elseif changeI(j) 
        
                 
                        xtot{j} = drmat0tot(CELLLAB_ALL == CLS(j),:) / Lambdatot_partialW_ChangeI;
                        xtot{j}(:, (find(~ismember(CELLLAB_ALL, changeWlin))) + length(CELLLAB_ALL)) = xtot{j}(:, ~ismember(CELLLAB_ALL, changeWlin));
                        xrel(CELLLAB_ALL == CLS(j),:)   = xtot{j}(:, [1:length(CELLLAB_ALL),  (2*length(CELLLAB_ALL) + 1):(2*length(CELLLAB_ALL) + 17 * Ncond)]);
                        xirrel(CELLLAB_ALL == CLS(j),:) = xtot{j}(:, [(length(CELLLAB_ALL)+1):(2*length(CELLLAB_ALL)), (2*length(CELLLAB_ALL) + 17 * Ncond + 1):end]);   

                        prediction_tot(CELLLAB_ALL == CLS(j),:) = xtot{j} * Lambdatot_partialW_ChangeI;

           end
    
    end

Measurement_rel = horzcat(drmat0rel{1:Ncond});
Measurement_irrel = horzcat(drmat0irrel{1:Ncond});
prediction_rel = prediction_tot(:, 1:size(Measurement_rel,2));
prediction_irrel = prediction_tot(:, (size(Measurement_rel,2) + 1):end);

    
end