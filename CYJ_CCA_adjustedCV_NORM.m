function CYJ_CCA_adjustedCV_NORM(tableX,brain_names,beh_names,p_ind,Res_ind,CV_names_brain,CV_names_beh,outpath,prefix,per_number)
% perform CCA between variables of the brain and variables of the behavior
% tableX: all data stored
% brain_names,beh_names: string cells£¬respective names of brain Vars and behavioral Vars£¬have to be row vector
% p_ind£ºstring£¬'pF'/'pChisq'/'permutation'
% Res_ind£ºwhether to take residuals of the covariates for the original Vars
% CV_names_brain,CV_names_beh: string cells,name of covariates
% per_number: number of permutation if p_ind is set to 'permutation'

table_brain=table;
table_beh=table;

if ~exist(outpath,'dir')
    mkdir(outpath);
end

if Res_ind==1
    % make sure brain beh CV no NaN
    [nonNaN_tableX,~]=specifiedVar_notNan(tableX,unique([beh_names,brain_names,CV_names_brain,CV_names_beh]));
    % take the residuals
    if ~isempty(CV_names_brain)
        for i=1:size(brain_names,2)
            formula=[brain_names{i} '~' strjoin(CV_names_brain,'+')];
            model=fitlm(nonNaN_tableX,formula);
            nonZ_brain_Res_matrix(:,i)=model.Residuals.Raw;
        end
    else
        for i=1:size(brain_names,2)
            eval(['nonZ_brain_Res_matrix(:,i)=nonNaN_tableX.' brain_names{i} ';']);
        end
    end
    if ~isempty(CV_names_beh)
        for i=1:size(beh_names,2)
            formula=[beh_names{i} '~' strjoin(CV_names_beh,'+')];
            model=fitlm(nonNaN_tableX,formula);
            nonZ_beh_Res_matrix(:,i)=model.Residuals.Raw;
        end
    else
        for i=1:size(beh_names,2)
            eval(['nonZ_beh_Res_matrix(:,i)=nonNaN_tableX.' beh_names{i} ';']);
        end
    end
else
    [nonNaN_tableX,~]=specifiedVar_notNan(tableX,[beh_names,brain_names]);
    for i=1:size(brain_names,2)
        eval(['nonZ_brain_Res_matrix(:,i)=nonNaN_tableX.' brain_names{i} ';']);
    end
    
    for i=1:size(beh_names,2)
        eval(['nonZ_beh_Res_matrix(:,i)=nonNaN_tableX.' beh_names{i} ';']);
    end
end

% normalize
brain_Res_matrix=zscore(nonZ_brain_Res_matrix);
beh_Res_matrix=zscore(nonZ_beh_Res_matrix);

% CCA
if strcmp(p_ind,'pF') || strcmp(p_ind,'pChisq')
    [brain_coef,beh_coef,CCC,U,V,stats]=canoncorr(brain_Res_matrix,beh_Res_matrix);
    save([outpath '/cca_coef'],'brain_coef','beh_coef');
elseif strcmp(p_ind,'permutation')
    [CCC,P_CCC,brain_loading,beh_loading,brain_loading_p,beh_loading_p,U,V,CCC_nullD]=CCA_permutation(brain_Res_matrix,beh_Res_matrix,per_number);
    hist(CCC_nullD);
    saveas(gcf,[outpath '/' prefix '_nullD'],'jpg');
    close(figure(gcf));
end



% plot
for i=1:length(CCC)
    if strcmp(p_ind,'pF') || strcmp(p_ind,'pChisq')
        %calculate loading
        [brain_loading(:,i),brain_loading_p(:,i)]=corr(brain_Res_matrix,U(:,i));
        [beh_loading(:,i),beh_loading_p(:,i)]=corr(beh_Res_matrix,V(:,i));
        eval(['p_value=stats.' p_ind '(i);']);
    elseif strcmp(p_ind,'permutation')
        p_value=P_CCC(i);
    end
    
    tmp_U=U(:,i);tmp_V=V(:,i);
    plot(tmp_U,tmp_V,'bo');
    xlabel('brain Adjusted');
    ylabel('beh Adjusted');
    
    contextFtitle{1}=strrep(prefix,'_','\_');
    contextFtitle{2}=['The ' num2str(i) ' element of r: ' num2str(CCC(i)) ' ' p_ind ': ' num2str(p_value,'%.4f')];
    contextFtitle{3}=['total var:' num2str(length(brain_names)+length(beh_names)) '  subject:' num2str(size(brain_Res_matrix,1))];
    title(contextFtitle);
    prefix2=[prefix '_No' num2str(i) 'elementOfR'];
    saveas(gcf,[outpath '/' prefix2],'jpg');
    close(figure(gcf));
    % save loading
    %if p_value<0.05
    table_brain.brain_VARname=brain_names';
    brain_loading_p1=brain_loading_p(:,i);
    %ind=find(brain_loading_p1>0.05);
    %brain_loading_p1(ind)=NaN;
    brain_loading1=brain_loading(:,i);
    %brain_loading1(ind)=NaN;
    %
    table_brain.brain_loading(:,i)=cellstr(num2str(brain_loading1,'%.3f'));
    table_brain.brain_loading_p(:,i)=cellstr(num2str(brain_loading_p1,'%.3f'));
    %
    table_beh.beh_VARname=beh_names';
    beh_loading_p1=beh_loading_p(:,i);
    %ind=find(beh_loading_p1>0.05);
    %beh_loading_p1(ind)=NaN;
    beh_loading1=beh_loading(:,i);
    %beh_loading1(ind)=NaN;
    %
    table_beh.beh_loading(:,i)=cellstr(num2str(beh_loading1,'%.3f'));
    table_beh.beh_loading_p(:,i)=cellstr(num2str(beh_loading_p1,'%.3f'));
    %     else
    %         table_brain.brain_loading(:,i)=cellfun(@(x) 'NaN',cell(size(brain_loading,1),1),'UniformOutput',false);
    %         table_brain.brain_loading_p(:,i)=cellfun(@(x) 'NaN',cell(size(brain_loading,1),1),'UniformOutput',false);
    %         table_beh.beh_loading(:,i)=cellfun(@(x) 'NaN',cell(size(beh_loading,1),1),'UniformOutput',false);
    %         table_beh.beh_loading_p(:,i)=cellfun(@(x) 'NaN',cell(size(beh_loading,1),1),'UniformOutput',false);
    %     end
end

writetable(table_brain,[outpath '/' prefix '_table_brain.csv']);
writetable(table_beh,[outpath '/' prefix '_table_beh.csv']);

%%
function [nonNaN_tableX,nonNaNs_ind]=CYJ_specifiedVar_notNan(tableX,var_names)
% deal with nan and empty
sum_ind=zeros(size(tableX,1),1);
for i=1:length(var_names)
    eval(['a=iscellstr(tableX.' var_names{i} ');']);
    if a
        eval(['ind=cellfun(@isempty,tableX.' var_names{i} ');']);
    else
        eval(['ind=isnan(tableX.' var_names{i} ');']);
    end
    sum_ind=sum_ind+ind;
end
nonNaNs_ind=~sum_ind>0;
nonNaN_tableX=tableX(nonNaNs_ind,:);

%%
function [real_CCC,P_CCC,real_Matrix1_loading,real_Matrix2_loading,P_Matrix1_loadings,P_Matrix2_loadings,real_U,real_V,CCC_nullD]=CCA_permutation(Matrix1,Matrix2,per_number)
% the first one is the one who going through permutation

% permutation CCA
for i=1:per_number
    rowrank = randperm(size(Matrix1, 1));
    PERed_Matrix=Matrix1(rowrank,:);
    [~,~,PERed_CCC(:,i),U,V,~]=canoncorr(PERed_Matrix,Matrix2);

    for j=1:size(PERed_CCC,1)
        %calculate loading
        PERed_Matrix1_loading(:,j)=corr(Matrix1,U(:,j));
        PERed_Matrix2_loading(:,j)=corr(Matrix2,V(:,j));
    end
    
    PERed_Matrix1_loadings(:,:,i)=PERed_Matrix1_loading;
    PERed_Matrix2_loadings(:,:,i)=PERed_Matrix2_loading;
end

% real CCA
[~,~,real_CCC,real_U,real_V,~]=canoncorr(Matrix1,Matrix2);
for i=1:length(real_CCC)
    real_Matrix1_loading(:,i)=corr(Matrix1,real_U(:,i));
    real_Matrix2_loading(:,i)=corr(Matrix2,real_V(:,i));
end

% test for CCC
CCC_nullD=PERed_CCC(1,:);% read value is not inside
cell_real_CCC=num2cell(real_CCC);
get_P=@(x) (length(find(abs(CCC_nullD)>=abs(x)))+1)/(per_number+1);
P_CCC=cellfun(get_P,cell_real_CCC);

% test for loadings
PERed_Matrix1_loadings(:,:,per_number+1)=real_Matrix1_loading;
PERed_Matrix2_loadings(:,:,per_number+1)=real_Matrix2_loading;

tmp_Matirx1=reshape(PERed_Matrix1_loadings,size(PERed_Matrix1_loadings,1)*size(PERed_Matrix1_loadings,2),size(PERed_Matrix1_loadings,3));
tmp_Matirx2=reshape(PERed_Matrix2_loadings,size(PERed_Matrix2_loadings,1)*size(PERed_Matrix2_loadings,2),size(PERed_Matrix2_loadings,3));

cell_PERed_Matrix1_loadings=mat2cell(tmp_Matirx1,ones(size(tmp_Matirx1,1),1),per_number+1);
cell_PERed_Matrix2_loadings=mat2cell(tmp_Matirx2,ones(size(tmp_Matirx2,1),1),per_number+1);

get_P=@(x) (length(find(abs(x)>=abs(x(end)))))/(per_number+1);
tmp_Matirx1=cellfun(get_P,cell_PERed_Matrix1_loadings);
tmp_Matirx2=cellfun(get_P,cell_PERed_Matrix2_loadings);

P_Matrix1_loadings=reshape(tmp_Matirx1,size(PERed_Matrix1_loadings,1),size(PERed_Matrix1_loadings,2));
P_Matrix2_loadings=reshape(tmp_Matirx2,size(PERed_Matrix2_loadings,1),size(PERed_Matrix2_loadings,2));


