%% Apply CT Corrections to Patient Data and Compare Uncorr vs Corr. Results

% Plot uncorrected Cherenkov vs. Dose
%% 6X
m6 = -9.59;
b6 = 450.4;
addpath('/Users/f004mg7/Library/CloudStorage/GoogleDrive-savannah.m.decker.th@dartmouth.edu/My Drive/Research/Scripts')
load cmap1 
for i = [11]
    for k = 1:3
        figure; imagesc(all6X(i).Cherenkov); colormap(cmap1)
        h = imellipse(); m = createMask(h);
        c = mean(nonzeros(double(all6X(i).Cherenkov).*m), 'all');
        d = mean(nonzeros(all6X(i).Dose.*m), 'all');
        i*3+k-1
        unC6_total(i*3+k-1) = c;
        d6_total(i*3+k-1) = d;
        close; 
        figure; imagesc(all6X(i).CT); colormap(gray)
        h = imellipse(); m = createMask(h);
        HU = nanmean(nonzeros((all6X(i).CT.*m)), 'all');
        CF = (m6.*(-135) + b6) / (m6*HU + b6);
        HU6_total(i*3+k-1) = HU;
        CF6_total(i*3+k-1) = CF;
        %count = count+1
        close 
    end
 
end

figure; plot(d6_total, unC6_total, 'ro'); lsline

% plot corrected CHerenkov vs. Dose
corr6 = unC6_total.*CF6_total;
hold on; plot(d6_total, corr6, 'bo'); lsline


%% calculate correction factor based on HU
% m6 = -9.59;
% b6 = 450.4;
% for i = 1:length(all6X)
%     figure; imagesc(all6X(i).CT); colormap(gray)
%     h = imellipse(); m = createMask(h);
%     HU = nanmean((all6X(i).CT.*m), 'all');
%     CF = (m6.*(-135) + b6) / (m6*HU + b6);
%     close 
%     HU6(i) = HU;
%     CF6(i) = CF;
% end

%% plot corrected vs uncorrected + fits
figure;
% scatter(d6, unC6, "filled", 20, 'r', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.8)
% hold on
% scatter(d6, corr6,"filled",  20, 'b', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.8)
plot(d6, unC6, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
hold on 
plot(d6, corr6, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b')
xlabel('Dose (cGy)')
ylabel('Cherenkov Intensity (a.u.)')
title('Moffitt - 6X')

x = linspace(0, 300, 300);
y1 = unC6_fit.p1.*x+ unC6_fit.p2;
y2= corr6_fit.p1 .*x + corr6_fit.p2;

hold on 
plot(x, y1, 'r--', 'LineWidth', 2)
hold on 
plot(x, y2, 'b--', 'LineWidth', 2)
set(gca, 'FontSize', 24)

str1=['Uncorr. R^2 =',num2str(unC6_goodness.rsquare)];
%'y = ',num2str(e_fit.p1),'x+',num2str(e_fit.p2)]
t=annotation('textbox',[.15 .9 0 0],'string',str1,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'r','LineWidth', 1)
t.FontSize = 20;

str2=[ 'Corr. R^2 =' num2str(corr6_goodness.rsquare)]
%'y = ',num2str(e_fit.p1),'x+',num2str(e_fit.p2)]
t=annotation('textbox',[.15 .7 0 0],'string',str2,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'b','LineWidth', 1)
t.FontSize = 20;

xlim([0 200])
legend({'Uncorrected', 'Corrected'})

%% 
count = 1;
for i = [1,10,17,19,21,25, 1,3,5,7,10,12, 17,19,21,25] % 15: [1,3,5,7,12,17,19,21,25, 1,3,12, 17,19,21]
    all6X(count).CT = ctData3(i).CTnotape;
    count = count+1;
end

%% 15X
m15 = -12.31;
b15 = 654.86;

for i = [5]%1:length(all15X)
    for k = 1:3
        figure; imagesc(all15X(i).Cherenkov); colormap(cmap1)
        h = imellipse(); m = createMask(h);
        c = mean(nonzeros(double(all15X(i).Cherenkov).*m), 'all');
        d = mean(nonzeros(all15X(i).Dose.*m), 'all');
        i*3+k-1
        unC15_total(i*3+k-1) = c;
        d15_total(i*3+k-1) = d;
        close; 
        figure; imagesc(all15X(i).CT); colormap(gray)
        h = imellipse(); m = createMask(h);
        HU = nanmean(nonzeros((all15X(i).CT.*m)), 'all');
        CF = (m15.*(-135) + b15) / (m15*HU + b15);
        HU15_total(i*3+k-1) = HU;
        CF15_total(i*3+k-1) = CF;
        %count = count+1
        close 
    end
 
end

figure; plot(d15_total, unC15_total, 'ro'); lsline

% plot corrected CHerenkov vs. Dose
corr15_total = unC15_total.*CF15_total;
hold on; plot(d15_total, corr15_total, 'bo'); lsline
%% plot corrected vs uncorrected + fits
figure;
% scatter(d6, unC6, "filled", 20, 'r', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.8)
% hold on
% scatter(d6, corr6,"filled",  20, 'b', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.8)
plot(d15, unC15, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
hold on 
plot(d15, corr15, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b')
xlabel('Dose (cGy)')
ylabel('Cherenkov Intensity (a.u.)')
title('Moffitt - 15X')
corr15
x = linspace(0, 300, 300);
y1 = uncorr15_fit.p1.*x+ uncorr15_fit.p2;
y2= corr15_fit.p1 .*x + corr15_fit.p2;

hold on 
plot(x, y1, 'r--', 'LineWidth', 2)
hold on 
plot(x, y2, 'b--', 'LineWidth', 2)
set(gca, 'FontSize', 24)

str1=['Uncorr. R^2 =',num2str(uncorr15_goodness.rsquare)];
%'y = ',num2str(e_fit.p1),'x+',num2str(e_fit.p2)]
t=annotation('textbox',[.15 .9 0 0],'string',str1,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'r','LineWidth', 1)
t.FontSize = 20;

str2=[ 'Corr. R^2 =' num2str(corr15_goodness.rsquare)]
%'y = ',num2str(e_fit.p1),'x+',num2str(e_fit.p2)]
t=annotation('textbox',[.15 .7 0 0],'string',str2,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'b','LineWidth', 1)
t.FontSize = 20;

xlim([0 200])
legend({'Uncorrected', 'Corrected'})
%% show all uncorrected vs corrected data
figure; 
s1 = scatter(d6_total, unC6_total, 120, 'red', 'filled')
s1.MarkerFaceAlpha = 0.5
hold on
s2 = scatter(d15_total, unC15_total, 120, 'blue', 'filled')
s2.MarkerFaceAlpha = 0.5

x = linspace(0, 300, 300);
y1 = uncorr6t_fit.p1.*x+ uncorr6t_fit.p2;
y2= uncorr15t_fit.p1 .*x + uncorr15t_fit.p2;

hold on 
plot(x, y1, 'r--', 'LineWidth', 2)
hold on 
plot(x, y2, 'b--', 'LineWidth', 2)
set(gca, 'FontSize', 24)
xlim([0 200])
xlabel('Dose (cGy)')
ylabel('Cherenkov Intensity (a.u.)')
title('Uncorrected Data')
legend({'6X', '15X'})
%%
figure; 
s1 = scatter(d6_total, corr6_total, 120, 'red', 'filled')
s1.MarkerFaceAlpha = 0.5
hold on
s2 = scatter(d15_total, corr15_total, 120, 'blue', 'filled')
s2.MarkerFaceAlpha = 0.5

x = linspace(0, 300, 300);
y1 = corr6t_fit.p1.*x+ corr6t_fit.p2;
y2= corr15t_fit.p1 .*x + corr15t_fit.p2;

hold on 
plot(x, y1, 'r--', 'LineWidth', 2)
hold on 

plot(x, y2, 'b--', 'LineWidth', 2)
set(gca, 'FontSize', 24)
xlim([0 200])
xlabel('Dose (cGy)')
ylabel('Cherenkov Intensity (a.u.)')
title('CT-Corrected Data')
legend({'6X', '15X'})

%% apply CT correction to Moffitt skin color data
