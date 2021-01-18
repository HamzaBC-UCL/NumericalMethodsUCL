%% Question 1
clear
%% Part a
n = 2000; % set sample size

u1 = rand(1,n); % random numbers
u2 = rand(1,n);

p = sqrt(-2*log(u1)); % transform u1 into radius form

theta = 2*pi*u2; % transform u2 into angle form

% cartesian coordinates
x = p.*cos(theta);
y = p.*sin(theta);

% transforming further into normal random numbers with mean and sd
sigma = 0.1;
mu = 0.2;

x1 = x*sigma+mu;
y1 = y*sigma+mu;

%% test plot of numbers, not part of the question
subplot(2,2,1);
histogram(x1,50,'FaceColor','b');
colormap hot; axis square
title(sprintf('Box-Muller Samples Y\n Mean = %1.5f\n Variance = %1.5f',round(mean(x1),2),round(var(x1),2)))

subplot(2,2,2);
histogram(y1,50,'FaceColor','b');
colormap hot; axis square
title(sprintf('Box-Muller Samples X\n Mean = %1.5f\n Variance = %1.5f',round(mean(y1),2),round(var(y1),2)))

%% Part b
grid = linspace(-0.2,0.6,201); % create vector of evenly spaced numbers
ran = rand(1,n); % random numbers

pd = 1/(sqrt(2*pi)*sigma)*exp(-((grid-mu)/sigma).^2/2); % pdf for steps

pdU = 1/(sqrt(2*pi)*sigma)*exp(-((x1-mu)/sigma).^2/2); % pdf for box muller random numbers we generated, uses one of the random numbers instead of 1/2 from each. Either way, it gives the same answer

% I liked the way the subplot looked so i kept it here. Only Running this
% section should result in only this plot
subplot(2,1,2)
plot(grid,pd,'g','LineWidth',2)
hold on
plot(x1,ran.*pdU,'b.')
xlim([-0.2,0.6]) 
title(sprintf('Analytical pdf vs Box Muller\n mu = 0.2\n sigma = 0.1'))
legend('Analytical pdf', 'Box Muller')