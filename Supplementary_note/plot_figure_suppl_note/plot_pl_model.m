function plot_pl_model
%% This code plots the figure in the Supplementary note

%%


fileID = fopen('gene_name.txt','r');
gene = fscanf(fileID,'%s');
fclose(fileID);


data = dlmread(['./fit_s_u_',gene,'.txt']);


data = sortrows(data,1);



t0 = 0.5;

data(data(:,1)<t0,:) = [];



data(1:6,:) = [];

n = size(data,1);

measured=data(:,3:-1:2);

m = max(measured);

measured = measured./repmat(m,n,1);





st = reshape(measured(:,2),n/10,10);
ut = reshape(measured(:,1),n/10,10);



mst = mean(st);
mut = mean(ut);


mst = mst/max(mst);
mut = mut/max(mut);



tt = reshape(data(:,1),n/10,10);
mtt = mean(tt);




err_u = 0.01445;
err_s = 0.01445;



best = dlmread('smim1_mod.txt');

phase1 = dlmread('bounds1.txt');
phase2 = dlmread('bounds2.txt');

time = best(:,1);




X = [time;flipud(time)]';
Y1 = [phase1(:,3);flipud(phase1(:,2))]';
Y2 = [phase2(:,3);flipud(phase2(:,2))]';


%% Plot

figure(1)
clf

subplot(2,1,1)

hold on


errorbar(mtt,mut,err_u' * ones(length(mtt),1),'ok','MarkerFaceColor','k')

plot(time,phase1(:,1),'r','Linewidth',2)

fill(X,Y1,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
plot(time,phase1(:,1),'r','Linewidth',2)

errorbar(mtt,mut,err_u' * ones(length(mtt),1),'ok','MarkerFaceColor','k')



ylabel('unspliced')

legend('data','model','95% c.i.','location','northwest')

title(gene)



subplot(2,1,2)


hold on

fill(X,Y2,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
plot(time,phase2(:,1),'r','Linewidth',2)

errorbar(mtt,mst,err_s' * ones(length(mtt),1),'ok','MarkerFaceColor','k')

ylabel('spliced')


xlabel('latent time')




% set(gcf, 'PaperUnits', 'centimeters');
% exportfig(gcf,'fit_smim1.eps','FontMode', 'fixed','Fontsize',20,'color', 'cmyk','width',13,'height',20,'Renderer','painters','Lockaxes',0);%

end