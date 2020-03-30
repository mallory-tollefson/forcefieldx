%% This script plots the lambda values experienced by individual walkers. The script is currently designed to work with 16 walkers or fewer.
%
% walkers.txt should consist of two columns of numbers, where column one
% is the walker number (0-15) and column two is the lambda value (0.000-1.000).
%
% To prepare the walkers.txt file from a raw FFX log file, use the
% prepareWalkers.sh script in the forcefieldx/matlab directory:
% ./prepareWalkers.sh MC-OST.log

clear all;
close all;
figure(1);
walkers = load('walkers.txt');
ln = size(walkers);
n=ln(1);

walkerNumbers = walkers(:,1);
lambdaValues = walkers(:,2);
currentStep = zeros(1,16);
colors = [255,0,0;255,128,0;128,255,0;0,255,255;0,0,255;127,0,255;255,0,255;150,150,150;255,255,255;255,255,0;115,60,60;35,105,75;225,175,0;100,0,0;100,175,255;150,255,200];
legends = zeros(1,16);

hold on;
for i = 1:n   
    walkerNumber = walkerNumbers(i) + 1;
    lambdaValue = lambdaValues(i);
    step = currentStep(walkerNumber);
    color = colors(walkerNumber,:);
    if currentStep(walkerNumber) == 0
        legends(walkerNumber) = plot(step,lambdaValue,'.','MarkerSize',30, 'Color', color/255);
    else
        plot(step,lambdaValue,'.','MarkerSize',30, 'Color', color/255);
    end
    currentStep(walkerNumber) = currentStep(walkerNumber) + 1;
end
hold off;
axis square
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
title('Trajactory \lambda States','FontSize', 16);
xlabel('Time (picoseconds)');
ylabel('\lambda');

numWalkers=0;
for i = 1:16
    if legends(i)==0.0
        numWalkers = i;
    end
    if numWalkers ~= 0
        numWalkers = i-1;
        break
    end
end

legendString="";
walkerString="";
for i = 1:numWalkers
    legendString = strcat(legendString, "legends(");
    legendString = strcat(legendString, num2str(i));
    
    walkerString = strcat(walkerString,"'Walker ");
    walkerString = strcat(walkerString, num2str(i));
    
    if i==numWalkers
        legendString = strcat(legendString, ")");
        walkerString = strcat(walkerString, "'");
    else
        legendString = strcat(legendString, "), ");
        walkerString = strcat(walkerString, "',");
    end
end

finalLegend = strcat("legend([",legendString);
finalLegend = strcat(finalLegend,"],{");
finalLegend = strcat(finalLegend,walkerString);
finalLegend = strcat(finalLegend,"});");
eval(finalLegend)