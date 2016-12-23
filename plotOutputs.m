function [ ] = plotOutputs( phat, half, time, efficiency, fileName )
% Plot the confidence interval of the estimated probability of convex,
% iteration time, and iteration efficiency vs. iteration number.

maxIte = length(phat);
num = [1:1:maxIte];

% Plot the estimated probability of convex
h1 = figure;
hold on;
axis([0 maxIte 0 1]);
plot(num, phat + half, '--');
area(phat + half, 'LineStyle', ':', 'FaceColor', 'y');
plot(num, phat - half,'--');
area(phat - half, 'LineStyle', ':', 'FaceColor', 'w');
plot(num, phat,'b', 'LineWidth', 2);
xlabel('Number of iterations');
ylabel('Estimated probability of convex');
figureHandle = gcf;
%# make all text in the figure to size 15 and bold
set(gca,'FontSize',40,'fontWeight','normal')
set(findall(figureHandle,'type','text'),'fontSize',40,'fontWeight','normal')
hold off;
saveas(h1, [pwd, '\outputs\', fileName, '_phat.fig']);

% Plot the estimated probability of convex
h2 = figure;
hold on;
axis auto;
plot(num, time,'b');
xlabel('Number of iterations');
ylabel('Iteration time');
figureHandle = gcf;
set(gca,'FontSize',40,'fontWeight','normal')
set(findall(figureHandle,'type','text'),'fontSize',40,'fontWeight','normal')
hold off;
saveas(h2, [pwd, '\outputs\', fileName, '_time.fig']);

% Plot the estimated probability of convex
h3 = figure;
hold on;
axis auto;
plot(num, efficiency,'b');
xlabel('Number of iterations');
ylabel('Iteration efficiency');
figureHandle = gcf;
set(gca,'FontSize',25,'fontWeight','normal')
set(findall(figureHandle,'type','text'),'fontSize',25,'fontWeight','normal')
hold off;
saveas(h3, [pwd, '\outputs\', fileName, '_efficiency.fig']);

end

