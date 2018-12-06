close all

T = ADDTRG{1}.P.t;
X = ADDTRG{1}.PL{1}.PSTH;
XERR = ADDTRG{1}.PL{1}.PSTHERR;

TMOT = ADDTRG{1}.P.tmot;
V = ADDTRG{1}.P.PSTHMOT;
VERR = ADDTRG{1}.P.PSTHMOTERR;

subplot(221)

T0 = find(T==0);

% Vertical grating, relevant vs irrelevant

TRIG = 1;

hold on
title('Visual Relevant', 'fontsize', 18)
ylabel('Population Averaged Response (dF/F)', 'fontsize', 18)
xlabel('Time From Stimulus (seconds)', 'fontsize', 18)

hPV = plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) + squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'r');
plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) - squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'r')
hSOM = plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) - squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'g');
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) + squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'g')
hVIP = plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) + squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'b');
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) - squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'b')

hPYR = plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) + squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'k');
plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) - squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'k')

TRIG = 2;

hold on

plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) + squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'color', 'r', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) - squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'color', 'r', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) - squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'color', 'g', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) + squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'color', 'g', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) + squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'color', 'b', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) - squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'color', 'b', 'linestyle', '--')

plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) + squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'color', 'k', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) - squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'color', 'k', 'linestyle', '--')

box on
set(gca, 'fontsize',  18)
legend([hPV, hSOM, hVIP, hPYR], {'PV', 'SOM', 'VIP', 'PYR'})


subplot(223)

box on

hold on

TRIG = 1;

plot(TMOT, V(TRIG,:) + VERR(TRIG,:), 'linewidth', 3)
plot(TMOT, V(TRIG,:) - VERR(TRIG,:), 'linewidth', 3)

TRIG = 2;

plot(TMOT, V(TRIG,:) + VERR(TRIG,:), 'linewidth', 3, 'linestyle', '--')
plot(TMOT, V(TRIG,:) - VERR(TRIG,:), 'linewidth', 3, 'linestyle', '--')

axis([min(T), max(T), min(min(V)), max(max(V))])

title('Visual Relevant', 'fontsize', 18)
ylabel('Running speed', 'fontsize', 18)
xlabel('Time From Stimulus (seconds)', 'fontsize', 18)

% axis([-0.1, 0.25, -0.06, 0.06]) 

subplot(222)

% Angled grating, relevant vs irrelevant

TRIG = 5;

hold on
title('Visual Irrelevant', 'fontsize', 18)
ylabel('Population Averaged Response (dF/F)', 'fontsize', 18)
xlabel('Time From Stimulus (seconds)', 'fontsize', 18)


hPV = plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) + squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'r');
plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) - squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'r')
hSOM = plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) - squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'g');
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) + squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'g')
hVIP = plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) + squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'b');
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) - squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'b')

hPYR = plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) + squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'k');
plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) - squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'k')


TRIG = 6;

hold on

plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) + squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'color', 'r', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==3,:)) - mean(X(TRIG, CELLLAB==3, T0)) - squeeze(XERR(TRIG,CELLLAB==3,:)),1), 'color', 'r', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) - squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'color', 'g', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==4,:)) - mean(X(TRIG, CELLLAB==4, T0)) + squeeze(XERR(TRIG,CELLLAB==4,:)),1), 'color', 'g', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) + squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'color', 'b', 'linestyle', '--')
plot(T, mean(squeeze(X(TRIG,CELLLAB==5,:)) - mean(X(TRIG, CELLLAB==5, T0)) - squeeze(XERR(TRIG,CELLLAB==5,:)),1), 'color', 'b', 'linestyle', '--')

plot(T, mean(squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) + squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'color', 'k', 'linestyle', '--')
plot(T, (squeeze(X(TRIG,CELLLAB==1,:)) - mean(X(TRIG, CELLLAB==1, T0)) - squeeze(XERR(TRIG,CELLLAB==1,:)),1), 'color', 'k', 'linestyle', '--')

box on
set(gca, 'fontsize',  18)
legend([hPV, hSOM, hVIP, hPYR], {'PV', 'SOM', 'VIP', 'PYR'})
% axis([-0.1, 0.25, -0.06, 0.06]) 

subplot(224)

box on

hold on

TRIG = 5;

plot(TMOT, V(TRIG,:) + VERR(TRIG,:), 'linewidth', 3)
plot(TMOT, V(TRIG,:) - VERR(TRIG,:), 'linewidth', 3)

TRIG = 6;

plot(TMOT, V(TRIG,:) + VERR(TRIG,:), 'linewidth', 3, 'linestyle', '--')
plot(TMOT, V(TRIG,:) - VERR(TRIG,:), 'linewidth', 3, 'linestyle', '--')

axis([min(T), max(T), min(min(V)), max(max(V))])


title('Visual Irrelevant', 'fontsize', 18)
ylabel('Running speed', 'fontsize', 18)
xlabel('Time From Stimulus (seconds)', 'fontsize', 18)
set(gca, 'fontsize',  18)

