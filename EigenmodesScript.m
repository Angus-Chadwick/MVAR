
if task='learning'

elseif task='switching'
    
xpre = xrel;
xpost = xirrel;
residualpre = residualrel;
residualpost = residualirrel;

end
    
cond = 'differential'

for rix=1:length(RXS)

          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
  
    
          if strmatch(cond, 'differential')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
    
          elseif strmatch(cond, 'Vert') 
              
            dIn_pre{rix} = MeanVpre;
            dIn_post{rix} = MeanVpost;
            
          elseif strmatch(cond, 'Ang') 
              
            dIn_pre{rix} = MeanApre;
            dIn_post{rix} = MeanApost;
            
            
          elseif strmatch(cond, 'Shuffled')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
            
            dIn_pre{rix} = dIn_pre{rix}(randperm(length(dIn_pre{rix})));
            dIn_post{rix} = dIn_post{rix}(randperm(length(dIn_post{rix})));
            
          end
    
          SNRvec = 1;
          if SNRvec
          
            Covres_pre{rix} = cov(residualpre{RXS(rix)}');
            Covres_post{rix} = cov(residualpost{RXS(rix)}');
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});

            dIn_pre{rix} = Covresinv_pre{rix} * dIn_pre{rix};
            dIn_post{rix} = Covresinv_post{rix} * dIn_post{rix};
            
          end
    
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre{rix}, Dpre{rix}] = eig(Apre{rix});
    [Vpost{rix}, Dpost{rix}] = eig(Apost{rix});
    
    Vinv_pre{rix} = inv(Vpre{rix});
    Vinv_post{rix} = inv(Vpost{rix});
    
    for i=1:size(Vinv_post{rix},1)
    
        Vinv_pre{rix}(i,:)  = Vinv_pre{rix}(i,:) ./ norm(Vinv_pre{rix}(i,:));
        Vinv_post{rix}(i,:) = Vinv_post{rix}(i,:) ./ norm(Vinv_post{rix}(i,:));
        
    end
   
    
    I = find(imag(diag(Dpre{rix})) == 0);
    Xpre{rix} = abs(Vinv_pre{rix}(I,:) * dIn_pre{rix}) / norm(dIn_pre{rix});
    D = diag(Dpre{rix});
    Ypre{rix} = D(I);
    
    I = find(imag(diag(Dpost{rix})) == 0);
    Xpost{rix} = abs(Vinv_post{rix}(I,:) * dIn_post{rix}) / norm(dIn_post{rix});
    D = diag(Dpost{rix});
    Ypost{rix} = D(I);
    
end

Xpre_tot = vertcat(Xpre{:});
Ypre_tot = vertcat(Ypre{:});
Xpost_tot = vertcat(Xpost{:});
Ypost_tot = vertcat(Ypost{:});

thetapre_tot = acos(Xpre_tot);
thetapost_tot = acos(Xpost_tot);


binstep = 0.01;
binwidth = 0.1;
Ybins = -1.25:binstep:0.05;



for i=1:length(Ybins)
    
    Xbin_pre(i) = mean(Xpre_tot( and(Ypre_tot >= Ybins(i) - binwidth, Ypre_tot < Ybins(i) + binwidth)));
    Xbin_post(i) = mean(Xpost_tot( and(Ypost_tot >= Ybins(i) - binwidth, Ypost_tot < Ybins(i) + binwidth)));
    
    resultantbin_pre(i) = mean(exp(1i * thetapre_tot( and(Ypre_tot >= Ybins(i) - binwidth, Ypre_tot < Ybins(i) + binwidth))));
    resultantbin_post(i) = mean(exp(1i * thetapost_tot( and(Ypost_tot >= Ybins(i) - binwidth, Ypost_tot < Ybins(i) + binwidth))));

    thetabin_pre(i) = angle(resultantbin_pre(i));
    thetabin_post(i) = angle(resultantbin_post(i));

    thetaSEMbin_pre(i) = sqrt((1 - abs(resultantbin_pre(i)))) ./ sqrt(sum( and(Ypre_tot >= Ybins(i) - binwidth, Ypre_tot < Ybins(i) + binwidth)));
    thetaSEMbin_post(i) = sqrt((1 - abs(resultantbin_post(i))))  ./ sqrt(sum(  and(Ypost_tot >= Ybins(i) - binwidth, Ypost_tot < Ybins(i) + binwidth)));

    
    Xbin_pre_sem(i) = std(Xpre_tot( and(Ypre_tot >= Ybins(i) - binwidth, Ypre_tot < Ybins(i) + binwidth))) ./ sqrt(sum( and(Ypre_tot >= Ybins(i) - binwidth, Ypre_tot < Ybins(i) + binwidth)));
    Xbin_post_sem(i) = std(Xpost_tot( and(Ypost_tot >= Ybins(i) - binwidth, Ypost_tot < Ybins(i) + binwidth))) ./ sqrt(sum( and(Ypost_tot >= Ybins(i) - binwidth, Ypost_tot < Ybins(i) + binwidth)));
    
end

hold on
errorbar(Ybins, Xbin_pre, Xbin_pre_sem, 'k')
errorbar(Ybins, Xbin_post, Xbin_post_sem, 'r')


% taupre_tot = -1./Ypre_tot;
% taupost_tot = -1./Ypost_tot;

% Ipre = find(Ypre_tot + 1 > 0);
% Ipost = find(Ypost_tot + 1 > 0);

Ipre = find(and(0 < Ypre_tot + 1, Ypre_tot + 1 < 0.99));
Ipost = find(and(0 < Ypost_tot + 1, Ypost_tot + 1 < 0.99));

taupre_tot = -125./log(Ypre_tot(Ipre)+1);
taupost_tot = -125./log(Ypost_tot(Ipost)+1);


binstep = 25;
taubins = 0:binstep:2750;
binwidth = 100;

clear XXbin_pre
clear XXbin_post
clear XXbin_pre_sem
clear XXbin_post_sem
clear thetabin_pre
clear thetabin_post
clear thetaSEMbin_pre
clear thetaSEMbin_post


tbins = 1
if tbins

for i=1:length(taubins)
        
    XXbin_pre(i) = mean(Xpre_tot( and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth)));
    XXbin_post(i) = mean(Xpost_tot( and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth)));
    
    resultantbin_pre(i) = mean(exp(1i * thetapre_tot( Ipre(and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth)))));
    resultantbin_post(i) = mean(exp(1i * thetapost_tot( Ipost(and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth)))));

    thetabin_pre(i) = angle(resultantbin_pre(i));
    thetabin_post(i) = angle(resultantbin_post(i));

    thetaSEMbin_pre(i) = sqrt((1 - abs(resultantbin_pre(i)))) ./ sqrt(sum( and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth)));
    thetaSEMbin_post(i) = sqrt((1 - abs(resultantbin_post(i))))  ./ sqrt(sum( and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth)));
    
    XXbin_pre_sem(i) = std(Xpre_tot( and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth))) ./ sqrt(sum( and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth)));
    XXbin_post_sem(i) = std(Xpost_tot( and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth))) ./ sqrt(sum( and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth)));
    
end

shadedErrorBar(taubins, thetabin_pre * 180/pi, thetaSEMbin_pre * 180/pi,'lineprops','k');
hold on
shadedErrorBar(taubins, thetabin_post * 180/pi, thetaSEMbin_post * 180/pi,'lineprops','b');
axis([0,1250,76,92])
set(gca, 'fontsize', 28)
xlabel('Time Constant (ms)')
ylabel('Angle from Stimulus-Input Differential')
plotrawdata=1;
if plotrawdata

    scatter(taupre_tot, 180 / pi * thetapre_tot(Ipre), 100, '+', 'markeredgecolor', 'k')
    scatter(taupost_tot, 180 / pi * thetapost_tot(Ipost), 100, '+', 'markeredgecolor', 'b')    
    
end

else
    
    thetabins = 82.5:0.15:89.5
    
    for i=1:length(thetabins)
        
    taubin_pre(i) = mean(taupre_tot( and(180 / pi * thetapre_tot(Ipre) >= thetabins(i) - 0.5, 180 / pi * thetapre_tot(Ipre) < thetabins(i) + 0.5)));
    taubin_post(i) = mean(taupost_tot( and(180 / pi * thetapost_tot(Ipost) >= thetabins(i) - 0.5, 180 / pi * thetapost_tot(Ipost) < thetabins(i) + 0.5)));
    
    taubin_pre_sem(i) = std(taupre_tot( and(180 / pi * thetapre_tot(Ipre) >= thetabins(i) - 0.5, 180 / pi * thetapre_tot(Ipre) < thetabins(i) + 0.5))) / sqrt(sum(and(180 / pi * thetapre_tot(Ipre) >= thetabins(i) - 0.5, 180 / pi * thetapre_tot(Ipre) < thetabins(i) + 0.5)));
    taubin_post_sem(i) = std(taupost_tot( and(180 / pi * thetapost_tot(Ipost) >= thetabins(i) - 0.5, 180 / pi * thetapost_tot(Ipost) < thetabins(i) + 0.5))) / sqrt(sum(and(180 / pi * thetapost_tot(Ipost) >= thetabins(i) - 0.5, 180 / pi * thetapost_tot(Ipost) < thetabins(i) + 0.5)));
    
    end

shadedErrorBar(thetabins, taubin_pre, taubin_pre_sem,'lineprops','k');
hold on
shadedErrorBar(thetabins, taubin_post, taubin_post_sem,'lineprops','b');
set(gca, 'fontsize', 18)
ylabel('Time Constant (ms)')
xlabel('Angle from Stimulus-Input Differential')
    
end


    
for i=1:6
Jpre = find(and(taupre_tot> (i-1) * 200, taupre_tot < i * 200));
Jpost = find(and(taupost_tot> (i-1) * 200, taupost_tot < i * 200));
p(i) = ranksum(thetapre_tot(Ipre(Jpre)), thetapost_tot(Ipost(Jpost)));
end

for i=1:6
Jpre = find(and(180 / pi * thetapre_tot(Ipre) > thetabins(i) - 0.5, 180 / pi * thetapre_tot(Ipre) > thetabins(i) + 0.5));
Jpost = find(and(180 / pi * thetapost_tot(Ipost) > thetabins(i) - 0.5, 180 / pi *  thetapost_tot(Ipost) > thetabins(i) + 0.5));
p(i) = ranksum(taupre_tot(Jpre), taupost_tot(Jpost));
end


for rix = 1:length(RXS)
for i=1:length(taubins)

    taupre{rix} = -1./Ypre{rix};
    taupost{rix} = -1./Ypost{rix};
  
    XXXbin_pre{rix}(i) = nanmean(Xpre{rix}( and(taupre{rix} >= taubins(i) - binwidth, taupre{rix} < taubins(i) + binwidth)));
    XXXbin_post{rix}(i) = nanmean(Xpost{rix}( and(taupost{rix} >= taubins(i) - binwidth, taupost{rix} < taubins(i) + binwidth)));
    
    XXXbin_pre_sem{rix}(i) = nanstd(Xpre{rix}( and(taupre{rix} >= taubins(i) - binwidth, taupre{rix} < taubins(i) + binwidth))) ./ sqrt(sum( and(taupre{rix} >= taubins(i) - binwidth, taupre{rix} < taubins(i) + binwidth)));
    XXXbin_post_sem{rix}(i) = nanstd(Xpost{rix}( and(taupost{rix} >= taubins(i) - binwidth, taupost{rix} < taubins(i) + binwidth))) ./ sqrt(sum( and(taupost{rix} >= taubins(i) - binwidth, taupost{rix} < taubins(i) + binwidth)));
    
end
end

for rix = 1:length(RXS)
subplot(2, 4,rix)
errorbar(taubins * 125, XXXbin_pre{rix}, XXXbin_pre_sem{rix}, 'k')
hold on
errorbar(taubins * 125, XXXbin_post{rix}, XXXbin_post_sem{rix}, 'r')
end





