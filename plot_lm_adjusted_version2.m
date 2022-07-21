function [GLM_T,GLM_P,partial_R,partial_P]=plot_lm_adjusted_version2(model,X_name_real,X_name_model,divDot_name,contextFtitle,outpath,prefix,xlim_array,ylim_array,adjustedCV_ind,paper_ind)

if ~exist(outpath,'dir')
    mkdir(outpath);
end

% obtain X_data
tableX=model.Variables;
designM=CYJ_myDesignMatrix(tableX,model.PredictorNames,0);
coefName=model.CoefficientNames(2:end);% get rid of intercept
coef=model.Coefficients.Estimate(2:end);% get rid of intercept
ind_X=find(strcmp(coefName,X_name_model));
X_data=designM(:,ind_X);
% obtain the Y_adjusted
eval(['Y_data=tableX.' model.ResponseName]);
if adjustedCV_ind==1
    ind_CV=find(~strcmp(coefName,X_name_model));
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
figure('visible','on');
if ~isempty(divDot_name)
    eval(['divDot_var=model.Variables.' divDot_name ';']);
    [content_G]=unique(divDot_var);
    for i=1:length(content_G)
        ind = find(ismember(divDot_var,content_G{i}));
        plot(X_data(ind),Y_used(ind),'ko','markersize',5,'MarkerFaceColor',[rand rand rand]);
        hold on
    end
    hold off
    legend(content_G','Location','Best');
else
    plot(X_data,Y_used,'ko','markersize',5,'MarkerFaceColor','r');
end
coef=model.Coefficients;
slope=table2array(coef(X_name_model,'Estimate'));
intercept=table2array(coef('(Intercept)','Estimate'))+mean(c);
lm_line=refline(slope,intercept);
set(lm_line,'color','black');
set(lm_line,'LineWidth',2);
GLM_T=table2array(coef(X_name_model,'tStat'));
GLM_P=table2array(coef(X_name_model,'pValue'));
% xMin=min(X_data)-0.00785;
% xMax=max(X_data)+0.00785;
% yMin=min(Y_used)-1.9;
% yMax=max(Y_used)+1.9;
% xlim([xMin xMax]);
% ylim([yMin,yMax]);
% ax = gca;
% ax.XAxis.Exponent = -2;
% set(gcf,'Units','centimeters','Position',[0 0 20 20]);
% %set(gca,'ytick',-2:4:18);
% %set(gca,'xtick',-0.006:0.004:0.02);
%%%
tmp=length(contextFtitle);
Y_name=model.ResponseName;
ind_CV=find(~strcmp(model.PredictorNames,X_name_model));
CV_names=model.PredictorNames(ind_CV);
[partial_R,partial_P]=myPartialCorr(tableX,X_name_real,Y_name,CV_names);
contextFtitle{tmp+1}=['partialCorr:' num2str(partial_R) ' P:' num2str(partial_P)];
contextFtitle{tmp+2}=['Sub:' num2str(model.NumObservations) '  T:' num2str(GLM_T) ' P:' num2str(GLM_P)];
title(contextFtitle);
if ~isempty(xlim_array) && ~isempty(ylim_array)
    xlim(xlim_array);
    ylim(ylim_array);
end
axis square
xlabel(X_name_real);
if adjustedCV_ind==0
    ylabel(strrep(Y_name,'_','\_'));
elseif adjustedCV_ind==1
    ylabel([strrep(Y_name,'_','\_') ' AdjustedForCV']);
end
if paper_ind==1
    print(gcf,'-depsc2','-r600',[outpath '/' prefix '.eps']);
else
    saveas(gcf,[outpath '/' prefix],'jpg');
end
%clf
close(figure(gcf));

%%
function [XMatrix]=CYJ_myDesignMatrix(tableX,X_names,intercept_ind)
XMatrix=[];
for i=1:length(X_names)
    eval(['str_ind=iscellstr(tableX.' X_names{i} '(1));']);
    if str_ind==1
        eval(['XMatrix=[XMatrix dummyVar_refGroup(tableX.' X_names{i} ')];']);
    else
        eval(['tmp2=tableX.' X_names{i} ';'])
        XMatrix=[XMatrix tmp2];
    end
end
if intercept_ind==1
    XMatrix=[ones(size(XMatrix,1),1),XMatrix];
end
%%
function [codedGroup]=dummyVar_refGroup(group)
uniqueG={group{1}};%using the first item as reference
for i=1:length(group)-1
    
    if isempty(find(strcmp(group{i+1},uniqueG)))
        uniqueG=[uniqueG;group{i+1}];
    end
end
remainLevels=uniqueG(2:end);
codedGroup=zeros(length(group),length(remainLevels));
for i=1:length(remainLevels)
   ind=find(ismember(group,remainLevels{i})==1);
   codedGroup(ind,i)=1;
end