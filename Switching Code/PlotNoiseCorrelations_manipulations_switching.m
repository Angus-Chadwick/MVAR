
%% Plot PSTHs and single-trial measured vs predicted traces for individual cells

% ctot_irrel = noise correlations between pairs of cells of each type
% ctot_post_rw = noise correlations after weights are removed (note: self weights are left intact)
% ctot_post_shuffres = noise correlations after residuals are shuffled
% across trials for each cell (so that noise correlations in residuals are destroyed)
% 
% clear all
% load('NoiseCorrelations_modelperformance.mat')
% load('NoiseCorrelations_modelperformance_shuffweights.mat')
% load('NoiseCorr_Indices.mat')
% 
% %% get unique cell pairs
% 
% for i=1:length(ctot_post)
%         UniquePair{i} = logical(UniquePair{i});
%         ctot_pre{i} = ctot_pre{i}(UniquePair{i});
%         ctot_post{i} = ctot_post{i}(UniquePair{i});
%         ctot_pre_rw{i} = ctot_pre_rw{i}(UniquePair{i});
%         ctot_post_rw{i} = ctot_post_rw{i}(UniquePair{i});
%         
%         for j=1:length(ctot_pre_shuffres)
%                 ctot_pre_shuffweights{j}{i} = ctot_pre_shuffweights{j}{i}(UniquePair{i});
%                 ctot_post_shuffweights{j}{i} = ctot_post_shuffweights{j}{i}(UniquePair{i});
%                 ctot_pre_shuffres{j}{i} = ctot_pre_shuffres{j}{i}(UniquePair{i});
%                 ctot_post_shuffres{j}{i} = ctot_post_shuffres{j}{i}(UniquePair{i});
%         end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First show results of weight deletion

%% Plot mean noise correlations

figure

subplot(2,1,1)

plotconds = 'plot_rel'
ScatterCond = 'Rel'


SEM_rel = [cellfun(@nanstd, ctot_rel) ./ sqrt(cellfun(@length, ctot_rel)); cellfun(@nanmean, ctot_rel_rw)./ sqrt(cellfun(@length, ctot_rel)); nanmean(cellfun(@nanmean, vertcat(ctot_rel_shuffres{:}))) ./ sqrt(cellfun(@length, ctot_rel))]' ;
MEAN_rel = [cellfun(@nanmean, ctot_rel); cellfun(@nanmean, ctot_rel_rw); nanmean(cellfun(@nanmean, vertcat(ctot_rel_shuffres{:})))]';

SEM_irrel = [cellfun(@nanstd, ctot_irrel) ./ sqrt(cellfun(@length, ctot_irrel)); cellfun(@nanmean, ctot_irrel_rw)./ sqrt(cellfun(@length, ctot_irrel)); nanmean(cellfun(@nanmean, vertcat(ctot_irrel_shuffres{:}))) ./ sqrt(cellfun(@length, ctot_irrel))]' ;
MEAN_irrel = [cellfun(@nanmean, ctot_irrel); cellfun(@nanmean, ctot_irrel_rw); nanmean(cellfun(@nanmean, vertcat(ctot_irrel_shuffres{:})))]';

if strmatch(plotconds, 'plot_rel')

MEAN = [MEAN_rel];
SEM = [SEM_rel];

elseif strmatch(plotconds, 'plot_irrel')

MEAN = [MEAN_irrel];
SEM = [SEM_irrel];

end

hold on


hb = bar(1:size(MEAN,1), MEAN);

for ib = 1:numel(hb)

      % Find the centers of the bars
      xData = get(get(hb(ib),'Children'),'XData');
      barCenters = mean(unique(xData,'rows'));

      errorbar(barCenters,MEAN(:,ib),SEM(:,ib),'k.', 'linewidth', 3);

end



set(gca, 'fontsize', 22)
set(gca, 'xtick', 1:10)
set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
ylabel('Average noise correlation')

if strmatch(plotconds, 'plot_rel')

legend('Data (rel)', 'Interaction Weights Removed (rel)', 'Residuals Shuffled (rel)')

elseif strmatch(plotconds, 'plot_irrel')

legend('Data (irrel)', 'Interaction Weights Removed (irrel)', 'Residuals Shuffled (irrel)')

end
    
axis([0.5, 10.5, -0.05, 0.25])
% rotateXLabels( gca, 45)
box on

set(hb(1), 'FaceColor', 'k');
set(hb(2), 'FaceColor', 'r');
set(hb(3), 'FaceColor', 'b');

%% Scatter plots

if strmatch(ScatterCond, 'Irrel')

subplot(2,2,3)

X = vertcat(ctot_irrel{:});
%Y = vertcat(ctot_irrel_shuffres{1}{:});
Y0 = vertcat(ctot_irrel_shuffres{:});
for i=1:size(Y0,2)
    Y{i} = nanmean(horzcat(Y0{:,i}),2);
end
Y = vertcat(Y{:});
    
Y(isnan(X)) = [];  % remove nans
X(isnan(X)) = [];
B = Y' / [X, ones(size(X))]'; % Linear fit 

scatter(X,Y, 10,'.', 'b')

hold on
set(gca, 'fontsize', 22)
plot([-1,1], [-1,1], 'color', 'k', 'linewidth', 3)
% plot([-1,1], B(2) + B(1) * [-1,1], 'linewidth', 3)

box on
xlabel('True noise correlation')
ylabel('Noise correlation after shuffling')
title('Residuals shuffled')


MSE_shuffres = nanmean( (X - Y).^2);  % mean squared error of noise correlation with residuals shuffled
MeanNC = mean(X);
STDNC = std(X);
MeanNC_shuffres = mean(Y);
STDNC_shuffres = std(Y);

legend(strcat('Mean square error = ', num2str(mean(MSE_shuffres),2)), 'location', 'northwest')

subplot(2,2,4)

X = vertcat(ctot_irrel{:});
Y = vertcat(ctot_irrel_rw{:});

Y(isnan(X)) = [];  % remove nans
X(isnan(X)) = [];
B = Y' / [X, ones(size(X))]'; % Linear fit 

scatter(X,Y, 10,'.', 'r')

hold on
set(gca, 'fontsize', 22)
plot([-1,1], [-1,1], 'color', 'k', 'linewidth', 3)
% plot([-1,1], B(2) + B(1) * [-1,1], 'linewidth', 3)

box on
xlabel('True noise correlation')
ylabel('Noise correlation after weight deletion')
title('Interaction weights deleted')


MSE_rw = nanmean( (X - Y).^2);  % mean squared error of noise correlation with weights shuffled
MeanNC_rw = mean(Y);
STDNC_rw = std(Y);


legend(strcat('Mean square error = ', num2str(mean(MSE_rw),2)), 'location', 'northwest')

elseif strmatch(ScatterCond, 'Rel')


subplot(2,2,3)

X = vertcat(ctot_rel{:});
Y0 = vertcat(ctot_rel_shuffres{:});
for i=1:size(Y0,2)
    Y{i} = nanmean(horzcat(Y0{:,i}),2);
end
Y = vertcat(Y{:});
    
Y(isnan(X)) = [];  % remove nans
X(isnan(X)) = [];
B = Y' / [X, ones(size(X))]'; % Linear fit 

scatter(X,Y, 10,'.', 'b')

hold on
set(gca, 'fontsize', 22)
plot([-1,1], [-1,1], 'color', 'k', 'linewidth', 3)
% plot([-1,1], B(2) + B(1) * [-1,1], 'linewidth', 3)

box on
xlabel('True noise correlation')
ylabel('Noise correlation after shuffling')
title('Residuals shuffled')


MSE_shuffres = nanmean( (X - Y).^2);  % mean squared error of noise correlation with residuals shuffled
MeanNC = mean(X);
STDNC = std(X);
MeanNC_shuffres = mean(Y);
STDNC_shuffres = std(Y);

legend(strcat('Mean square error = ', num2str(mean(MSE_shuffres),2)), 'location', 'northwest')

subplot(2,2,4)

X = vertcat(ctot_rel{:});
Y = vertcat(ctot_rel_rw{:});

Y(isnan(X)) = [];  % remove nans
X(isnan(X)) = [];
B = Y' / [X, ones(size(X))]'; % Linear fit 

scatter(X,Y, 10,'.', 'r')

hold on
set(gca, 'fontsize', 22)
plot([-1,1], [-1,1], 'color', 'k', 'linewidth', 3)
% plot([-1,1], B(2) + B(1) * [-1,1], 'linewidth', 3)

box on
xlabel('True noise correlation')
ylabel('Noise correlation after weight deletion')
title('Interaction weights deleted')


MSE_rw = nanmean( (X - Y).^2);  % mean squared error of noise correlation with weights shuffled
MeanNC_rw = mean(Y);
STDNC_rw = std(Y);


legend(strcat('Mean square error = ', num2str(mean(MSE_rw),2)), 'location', 'northwest')

end

%%

clear X
clear Y

X = ctot_rel;
Y0 = vertcat(ctot_rel_shuffres{:});
for i=1:10
Y{i} = nanmean(horzcat(Y0{:,i}),2);
MSE_shuffres_rel(i) = nanmean((X{i} - Y{i}).^2);
Rsq_shuffres_rel(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end
Y = ctot_rel_rw;
for i=1:10
MSE_rw_rel(i) = nanmean((X{i} - Y{i}).^2);
Rsq_rw_rel(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end

clear X
clear Y

X = ctot_irrel;
Y0 = vertcat(ctot_irrel_shuffres{:});
for i=1:10
Y{i} = nanmean(horzcat(Y0{:,i}),2);
MSE_shuffres_irrel(i) = nanmean((X{i} - Y{i}).^2);
Rsq_shuffres_irrel(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end
Y = ctot_irrel_rw;
for i=1:10
MSE_rw_irrel(i) = nanmean((X{i} - Y{i}).^2);
Rsq_rw_irrel(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end

figure

% bar([MSE_rw_irrel;MSE_shuffres_irrel;MSE_rw_rel;MSE_shuffres_rel]')
% set(gca, 'fontsize', 22)
% set(gca, 'xtick', 1:10)
% set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
% legend('Delete Weights, Olfactory', 'Shuffle Residuals, Olfactory', 'Delete Weights, Visual', 'Shuffle Residuals, Visual')
% Y{i} = nanmean(horzcat(Y0{:,i}),2);
% ylabel('Mean Squared Error of Noise Correlations')
% axis([0.5, 10.5, 0, 0.04])
% title('Weights Fixed Over Switching')

bar([Rsq_rw_irrel;Rsq_shuffres_irrel;Rsq_rw_rel;Rsq_shuffres_rel]')
set(gca, 'fontsize', 22)
set(gca, 'xtick', 1:10)
set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
legend('Delete Weights, Olfactory', 'Shuffle Residuals, Olfactory', 'Delete Weights, Visual', 'Shuffle Residuals, Visual')
Y{i} = nanmean(horzcat(Y0{:,i}),2);
ylabel('R^2 of Noise Correlations')
axis([0.5, 10.5, 0, 0.04])
title('Weights Fixed Over Switching')



for i=1:10
X{i} = ctot_rel{i} - ctot_irrel{i};
end
Y0irrel = vertcat(ctot_irrel_shuffres{:});
Y0rel = vertcat(ctot_rel_shuffres{:});
for i=1:10
Y{i} = nanmean(horzcat(Y0rel{:,i}),2) -  nanmean(horzcat(Y0irrel{:,i}),2);
Rsq_DeltaNC_shuffres(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end

clear Y 
for i=1:10
Y{i} = ctot_rel_rw{i} - ctot_irrel_rw{i};
Rsq_DeltaNC_rw(i) = 1 - nansum((X{i} - Y{i}).^2) / nansum((X{i} - nanmean(X{i})).^2);
end


%% Print all stats

disp(strcat('Mean NC =', num2str(MeanNC)));
disp(strcat('Mean NC shuffled residuals =', num2str(MeanNC_shuffres)))
disp(strcat('Mean NC deleted weights =', num2str(MeanNC_rw)))
%disp(strcat('Mean NC shuffled weights =', num2str(MeanNC_shuffweights)))
disp(strcat('std NC =', num2str(STDNC)))
disp(strcat('std NC shuffled residuals =', num2str(STDNC_shuffres)))
disp(strcat('std NC deleted weights =', num2str(STDNC_rw)))
%disp(strcat('std NC shuffled weights =', num2str(STDNC_shuffweights)))
disp(strcat('Mean squared error shuffled residuals =', num2str(MSE_shuffres)))
disp(strcat('Mean squared error deleted weights =', num2str(MSE_rw)))
%disp(strcat('Mean squared error shuffled weights =', num2str(MSE_shuffweights)))
