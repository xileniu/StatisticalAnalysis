function CYJ_lvOneOut_valida_corr(tableX,x_name,y_name,CV_names,permuNum,outpath,prefix,xlim_array,ylim_array,paper_ind)
if ~exist(outpath,'dir')
    mkdir(outpath);
end

eval(['y_data=tableX.' y_name ';']);
XMatrix=myDesignMatrix(tableX,[x_name,CV_names],0);% need not intercept
for i=1:size(tableX,1)
    XMatrix_tmp=XMatrix;
    theOne_X=XMatrix_tmp(i,:);
    XMatrix_tmp(i,:)=[];
    y_tmp=y_data;
    y_tmp(i,:)=[];
    % real
    lm=fitlm(XMatrix_tmp,y_tmp);
    [predictedY(i,1),~]=predict(lm,theOne_X);
    % permutation
    predictedY_permutation(i,:)=LvOneOut_valida_corr_permutation(XMatrix_tmp,y_tmp,theOne_X,permuNum);
end

eval(['realY=tableX.' y_name ';']);
[realR,~]=corr(predictedY,realY);
permu_R=corr(predictedY_permutation,realY);
pValue=(length(find(abs(permu_R)>=abs(realR)))+1)/(permuNum+1);
tableX_run=table;tableX_run.predictedY=predictedY;tableX_run.realY=realY;
model=fitlm(tableX_run,'realY~predictedY');
plot_lm_adjusted_version2(model,'predictedY','predictedY','',{strrep(prefix,'_','\_'),['permutation p: ' num2str(pValue,'%.3f')]},outpath,prefix,xlim_array,ylim_array,0,paper_ind);


function [predictedY_permu]=LvOneOut_valida_corr_permutation(XMatrix_permu,y_data_permu,theOne_X,permuNum)
for i=1:permuNum
    rowrank = randperm(length(y_data_permu));
    y_tmp_permu=y_data_permu(rowrank);
    lm=fitlm(XMatrix_permu,y_tmp_permu);
    [predictedY_permu(i,1),~]=predict(lm,theOne_X);
end