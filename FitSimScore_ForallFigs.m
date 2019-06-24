%%% Analysis of DG network data %%%
% This Matlab code fits the data using least-squares methods:
% "Least squares" means that the overall solution minimizes the sum of the squares of the errors made in the results of every single data point
% http://en.wikipedia.org/wiki/Least_squares
% Enter the files to import.

% ModelDB file along with publication:
% Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
% http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract
% modified and augmented by
% Man Yi Yim / 2015
% Alexander Hanuschkin / 2011

close all;
set(0,'defaultaxesfontsize',10);

offset_corrected = 4;   % 0= not corrected; 1= corrected <AH> 2= corrected <MY> 3= two parameters (formula 1) 4= two parameters (formula 2) 5= three parameters
simscore1_corrected = 0;
fig1 = figure(1);

idname = '-pp10-gaba1-kir1-st0';
for figure_nr = 1:1
        A = importdata(strcat('sim_score',idname,'.txt'));

        % power law fit (http://en.wikipedia.org/wiki/Power_function)
        fun = @(x,xdata)xdata.^x(1);    % function to fit
        x0 = .1;                        % initial values for fitting parameters
        xdata = A(:,1);                 % x values data
        ydata = A(:,2);                 % y values data
        
        switch offset_corrected
            case 1
                ydata = ydata - (mean(ydata(find(xdata<0))));
                if (simscore1_corrected==1)
                    ydata = ydata * 1/ydata(find(A(:,1)==1));
                end
            case 2
                % power law fit (http://en.wikipedia.org/wiki/Power_function)
                % offset  
                offset = mean(ydata(find(xdata<0)));
                fun = @(x,xdata)offset+(1-offset)*xdata.^x(1);    % function to fit
            case 3
                fun = @(x,xdata)x(1)*(1-xdata)+xdata.^x(2);
                x0 = [0.1, 3];
            case 4
                fun = @(x,xdata)x(1) + (1-x(1))*xdata.^x(2);
                x0 = [0.1, 1];
            case 5
                %offset = mean(ydata(find(xdata<0)));
                fun = @(x,xdata)x(1)*(1-xdata).^x(2)+xdata.^x(3);
                x0 = [0.1, 1, 3];
        end
        
        
        x =  lsqcurvefit(fun,x0,xdata,ydata);   %
        line([0 1],[0 1],'LineStyle','--','Color','blue','LineWidth',3,'MarkerSize',1); hold on;
        line([0.6 0.6],[0 1],'LineStyle','-.','Color','green','LineWidth',3,'MarkerSize',1); hold on;
        h1 = plot([0:0.005:1],fun(real(x),[0:0.005:1]),'-','Color','red','LineWidth',3,'MarkerSize',1);hold on;
        ts = sprintf('fit x^{%.3f}',x);
        h2 = plot(xdata,ydata,'o','Color',[1 0.50 0.25],'LineWidth',2,'MarkerSize',10);hold on;
        xlabel('Input');
        xlim([-0.05 1.05]);
        ylabel('Output');
        ylim([-0.05 1.05]);
        title_str = sprintf('Similarity Score Fig. %d',figure_nr);
        title(title_str);

        ts = sprintf('%.3f',x);
        text(0.2,0.5,ts,'FontSize',16,'Color','red')
        ts = sprintf('%.3f',fun(real(x),0.6));
        text(0.2,0.8,ts,'FontSize',16,'Color','green')
        axis square
end

switch offset_corrected
    case 0
        text(-0.5,4.35,'(sim_{in})^c','FontSize',22,'Color','black');
    case 1
        text(-0.5,4.35,'(sim_{in})^c - mean subtraced','FontSize',22,'Color','black');
    case 2
        text(-1.3,4.35,'<sim_{out}(0)>+(1-<sim_{out}(0)>)(sim_{in})^c','FontSize',22,'Color','black');
    case 3
        text(-1.3,4.35,'a(1-sim_{in})+(sim_{in})^c','FontSize',22,'Color','black');
    case 4
        text(-1.3,4.35,'a+(1-a)(sim_{in})^c','FontSize',22,'Color','black');
    case 5
        text(-1.3,4.35,'a(1-sim_{in})^b+(sim_{in})^c','FontSize',22,'Color','black');
end