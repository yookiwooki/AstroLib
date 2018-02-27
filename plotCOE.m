%Generates Plot of Classical Orbital Elements over Time
%Not plotting fast variable (true or mean anomaly)

%Sean McArdle
%1/19/2017

function plotCOE(COEs,t,title)

%Plot COE's over time
fig1 = figure('Position', [100, 50, 1000, 950]);

subplot(5,1,1)
plot(t,COEs(:,1))
ylabel('a (km)')

subplot(5,1,2)
plot(t,COEs(:,2))
ylabel('e (unitless)')

subplot(5,1,3)
plot(t,COEs(:,3))
ylabel('i (radians)')

subplot(5,1,4)
plot(t,COEs(:,4))
ylabel('\omega (radians)')

subplot(5,1,5)
plot(t,COEs(:,5))
ylabel('\Omega (radians)')
xlabel('Time (sec)')

suptitle(title)

saveas(fig1, [pwd '\Pictures\' title '.png'])

end