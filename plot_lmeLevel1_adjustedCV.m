function [Tvalue,Pvalue]=plot_lmeLevel1_adjustedCV(lme,groupingVname,X_name,Y_name,contextFtitle,outPath,prefix,adjustedCV_ind,paper_ind,visible_ind,xlim_array,ylim_array)
% X_name have to be exactly the name in lme
% NaN in the table of fitting lme doesn't trouble this script
if ~exist(outPath,'dir')
    mkdir(outPath);
end

% obtain X_data
coefName=lme.CoefficientNames(2:end);% get rid of intercept
coef=lme.Coefficients.Estimate(2:end);% get rid of intercept
designM=designMatrix(lme);designM=designM(:,2:end);% get rid of intercept
ind_X=find(strcmp(coefName,X_name));
X_data=designM(:,ind_X);
% obtain the Y_adjusted
Y_data=response(lme);
if adjustedCV_ind==1
    ind_CV=find(~strcmp(coefName,X_name));
    CV_coef=coef(ind_CV);
    CV_designM=designM(:,ind_CV);
    CV=(CV_coef'*CV_designM')';
    c=nanmean(CV_designM,1);
    c=repmat(c,size(CV_designM,1),1);
    c=(CV_coef'*c')';
    Y_used=Y_data-CV+c;
elseif adjustedCV_ind==0
    Y_used=Y_data;
    c=0;
end
% obtain the slope and intercept using to plot refline
intercept=lme.Coefficients.Estimate(1)+mean(c);
slope=coef(ind_X);
% obtain T and P
Tvalue=lme.Coefficients.tStat(2:end);
Pvalue=lme.Coefficients.pValue(2:end);
Tvalue=Tvalue(ind_X);
Pvalue=Pvalue(ind_X);
% plot
figure('visible',visible_ind);%adjust it yourself
eval(['groupingV_data=lme.Variables.' groupingVname ';']);
[content_G]=unique(groupingV_data);
subWithNoodles=0;
for i=1:length(content_G)
    ind = find(ismember(groupingV_data,content_G{i}));
    if length(ind)>1
        table_tmp=table;
        table_tmp.X=X_data(ind);
        table_tmp.Y=Y_used(ind);
        table_tmp(arrayfun(@isnan,table_tmp.Y),:)=[];
        table_tmp(arrayfun(@isnan,table_tmp.X),:)=[];
        if size(table_tmp,1)>1
            plot(table_tmp.X,table_tmp.Y,'LineWidth',1.5,'color',[0.80392,0.36078,0.36078]);%line of light red
            subWithNoodles=subWithNoodles+1;
            hold on;
        end
    end
end

if ~isempty(xlim_array) && ~isempty(ylim_array)
    xlim(xlim_array);
    ylim(ylim_array);
end

plot(X_data,Y_used,'ko','markersize',4,'MarkerFaceColor',[0.69804,0.13333,0.13333]);%dot of red
hold on;
lme_line=refline(slope,intercept);
set(lme_line,'color','black');
set(lme_line,'LineWidth',5);
%%%%%%%%%%%%% adjust it yourself %%%%%%%%%%%%%%
%set(gca,'box','OFF');
 set(gcf,'Units','centimeters','Position',[0 0 12 20]);
 set(gca,'Position',[.13 .05 .8 .90]);% set axes position
%%%%%%%%%%%%% adjust it yourself %%%%%%%%%%%%%%
xlabel(X_name);
if adjustedCV_ind==0
    ylabel(Y_name);
elseif adjustedCV_ind==1
    ylabel([Y_name ' AdjustedForCV']);
end
%%%%
tmp=length(contextFtitle);
contextFtitle{tmp+1}=['T:' num2str(Tvalue) ' P:' num2str(Pvalue) ' observa:' num2str(lme.NumObservations) ' noodles:' num2str(subWithNoodles)];
title(contextFtitle);
%%%%
if paper_ind==1
    print(gcf,'-depsc2','-r600',[outPath '/' prefix '.eps']);
else
   saveas(gcf,[outPath '/' prefix],'jpg');
end
clf;
close(figure(gcf));