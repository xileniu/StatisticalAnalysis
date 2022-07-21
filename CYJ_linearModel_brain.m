function [DFE,T_0001,T_bonferroni_005,T_005]=CYJ_linearModel_brain(numberOfworkers,fileName_cell,tableX,formula, outpath, Mask_fileName,linearModel_flag,parpool_flag)
%fileName_cell:image of subjects, can not have empty
%tableX:a table contain all variables with variableNames
%formula:Y~X+C
%Mask_fileName:mask
%linearModel_flag:LME for linear mixed model,GLM for general linear model
if ~isempty(find(cellfun(@isempty,fileName_cell), 1))
    disp('There are empty cells in your fileName_cell, please check');
    return
end
if ~exist(outpath,'dir')
    mkdir(outpath);
end
if exist([outpath '/diary.txt'],'file')
    cmd=['rm -rf ' outpath '/diary.txt'];
    system(cmd);
end
diary([outpath '/diary.txt']);
diary on;

num_subjects=length(fileName_cell);
Nii1=load_untouch_nii(fileName_cell{1});
header=Nii1.hdr;
if header.dime.dim(1)==3
    dimension=header.dime.dim(2:4);
else
    error('the images have to be 3-dimensional');
end
slice_num=dimension(3);
disp(cat(2, 'There are totally ',num2str(slice_num), ' slices'));
disp('         ');

% to get df and para
tableX.volumeData=rand(num_subjects,1);
if strcmp(linearModel_flag,'GLM')==1
    model=fitlm(tableX,formula);
elseif strcmp(linearModel_flag,'LME')==1
    model=fitlme(tableX,formula);
end
DFE=model.DFE;
para=model.CoefficientNames;
% initialize the matrix
MSE_image=zeros(dimension(1),dimension(2),dimension(3));
AdjustedR_image=zeros(dimension(1),dimension(2),dimension(3));
Residuals4D_image=zeros(dimension(1),dimension(2),dimension(3),num_subjects);
T_image=zeros(dimension(1),dimension(2),dimension(3),length(para));
P_image=zeros(dimension(1),dimension(2),dimension(3),length(para));
if strcmp(parpool_flag,'matlabpool')
    matlabpool('local',numberOfworkers);
elseif strcmp(parpool_flag,'parpool')
    myparpool=parpool('local',numberOfworkers);
end
%
parfor i=1:slice_num
    data=[];
    mask_nii = load_untouch_nii(Mask_fileName, [], [], [], [], [], [i]);%load this slice of mask
    mask_data=mask_nii.img;
    index = find(mask_data);
    fileName_cell_parfor=fileName_cell;
    if ~isempty(index)
        for j=1:num_subjects
            tmp=load_untouch_nii(fileName_cell_parfor{j}, [], [], [], [], [], [i]);%load this slice of subject
            data(j,:)=tmp.img(index);% row is allvoxel of one subject,column is one voxel across all subjects
        end
    else
        disp(cat(2, num2str(i), 'th data is empty'));
        continue; % leave this function
    end
    disp(cat(2, 'processing slice ', num2str(i), 'th data, this slice has ',num2str(length(data(1,:))),' voxel'));
    
    tableX_parfor=tableX;
    T_slice_para=zeros(length(data(1,:)),length(para));
    P_slice_para=zeros(length(data(1,:)),length(para));
    MSE_slice=zeros(length(data(1,:)),1);
    AdjustedR_slice=zeros(length(data(1,:)),1);
    Residuals4D_slice=zeros(length(data(1,:)),num_subjects);
    if strcmp(linearModel_flag,'GLM')==1
        for voxel=1:length(data(1,:))
            tableX_parfor.volumeData = data(:,voxel);
            lm = fitlm(tableX_parfor,formula);
            MSE_slice(voxel)=lm.MSE;
            AdjustedR_slice(voxel)=lm.Rsquared.Adjusted;
            Residuals4D_slice(voxel,:)=lm.Residuals.Raw;
            for p = 1:length(lm.CoefficientNames)
                T_slice_para(voxel,p)= lm.Coefficients{p,3};%tstat
                P_slice_para(voxel,p)= lm.Coefficients{p,4};%p value
            end
            disp(['slice_' num2str(i) ' voxel_' num2str(voxel) ' is done']);
        end
    elseif strcmp(linearModel_flag,'LME')==1
        for voxel=1:length(data(1,:))
            tableX_parfor.volumeData = data(:,voxel);
            lme = fitlme(tableX_parfor,formula);
            MSE_slice(voxel)=lme.MSE;
            AdjustedR_slice(voxel)=lme.Rsquared.Adjusted;
            Residuals4D_slice(voxel,:)=residuals(lme);
            for p = 1:length(lme.CoefficientNames)
                T_slice_para(voxel,p)= lme.Coefficients{p,4};%tstat
                P_slice_para(voxel,p)= lme.Coefficients{p,6};%p value
            end
            disp(['slice_' num2str(i) ' voxel_' num2str(voxel) ' is done']);
        end
    end
    
    % MSE
    dimension_parfor=dimension;
    tmp = zeros(dimension_parfor(1), dimension_parfor(2));
    tmp(index)=MSE_slice;
    MSE_image(:,:,i)=tmp;
    % AdjustedR
    tmp = zeros(dimension_parfor(1), dimension_parfor(2));
    tmp(index)=AdjustedR_slice;
    AdjustedR_image(:,:,i)=tmp;
    % Residual4D
    tmp = zeros(dimension_parfor(1), dimension_parfor(2),num_subjects);
    tmp_4D=reshape(tmp,dimension_parfor(1)*dimension_parfor(2),num_subjects);
    tmp_4D(index,:)=Residuals4D_slice;
    tmp_4DD=reshape(tmp_4D,dimension_parfor(1), dimension_parfor(2),num_subjects);
    Residuals4D_image(:,:,i,:)=tmp_4DD;
    % T
    tmp = zeros(dimension_parfor(1), dimension_parfor(2),length(para));
    tmp_T=reshape(tmp,dimension_parfor(1)*dimension_parfor(2),length(para));
    tmp_T(index,:)=T_slice_para;
    tmp_TT=reshape(tmp_T,dimension_parfor(1), dimension_parfor(2),length(para));
    T_image(:,:,i,:)=tmp_TT
    % P
    tmp = zeros(dimension_parfor(1), dimension_parfor(2),length(para));
    tmp_P=reshape(tmp,dimension_parfor(1)*dimension_parfor(2),length(para));
    tmp_P(index,:)=P_slice_para;
    tmp_PP=reshape(tmp_P,dimension_parfor(1), dimension_parfor(2),length(para));
    P_image(:,:,i,:)=tmp_PP
end
if strcmp(parpool_flag,'matlabpool')
    matlabpool close
elseif strcmp(parpool_flag,'parpool')
    delete(myparpool);
end

% save T nifti file
T_nii=Nii1;
T_nii.hdr.hk.sizeof_hdr= 348;			% must be 348!
T_nii.hdr.hk.data_type= '';
T_nii.hdr.hk.db_name= '';
T_nii.hdr.hk.extents= 0;
T_nii.hdr.hk.session_error= 0;
T_nii.hdr.hk.regular= 'r';
T_nii.hdr.hk.dim_info= 0;
T_nii.hdr.dime.datatype=16;
T_nii.hdr.dime.bitpix=32;
T_nii.fileprefix=outpath;
T_nii.hdr.dime.scl_slope = 1;
for p=1:length(para)
    T_nii.img=T_image(:,:,:,p);
    save_untouch_nii(T_nii, cat(2,outpath,'/T_',para{p},'.nii'));
end
% save P nifti file
P_nii=Nii1;
P_nii.hdr.hk.sizeof_hdr= 348;			% must be 348!
P_nii.hdr.hk.data_type= '';
P_nii.hdr.hk.db_name= '';
P_nii.hdr.hk.extents= 0;
P_nii.hdr.hk.session_error= 0;
P_nii.hdr.hk.regular= 'r';
P_nii.hdr.hk.dim_info= 0;
P_nii.hdr.dime.datatype=16;
P_nii.hdr.dime.bitpix=32;
P_nii.fileprefix=outpath;
P_nii.hdr.dime.scl_slope = 1;
for p=1:length(para)
    P_nii.img=P_image(:,:,:,p);
    save_untouch_nii(P_nii, cat(2,outpath,'/P_',para{p},'.nii'));
end
% save MSE_image nifti file
MSE_image_nii=Nii1;
MSE_image_nii.hdr.hk.sizeof_hdr= 348;			% must be 348!
MSE_image_nii.hdr.hk.data_type= '';
MSE_image_nii.hdr.hk.db_name= '';
MSE_image_nii.hdr.hk.extents= 0;
MSE_image_nii.hdr.hk.session_error= 0;
MSE_image_nii.hdr.hk.regular= 'r';
MSE_image_nii.hdr.hk.dim_info= 0;
MSE_image_nii.hdr.dime.datatype=16;
MSE_image_nii.hdr.dime.bitpix=32;
MSE_image_nii.fileprefix=outpath;
MSE_image_nii.hdr.dime.scl_slope = 1;
MSE_image_nii.img=MSE_image;
save_untouch_nii(MSE_image_nii,cat(2,outpath, '/MSE.nii'));

% save AdjustedR_image nifti file
AdjustedR_image_nii=Nii1;
AdjustedR_image_nii.hdr.hk.sizeof_hdr= 348;			% must be 348!
AdjustedR_image_nii.hdr.hk.data_type= '';
AdjustedR_image_nii.hdr.hk.db_name= '';
AdjustedR_image_nii.hdr.hk.extents= 0;
AdjustedR_image_nii.hdr.hk.session_error= 0;
AdjustedR_image_nii.hdr.hk.regular= 'r';
AdjustedR_image_nii.hdr.hk.dim_info= 0;
AdjustedR_image_nii.hdr.dime.datatype=16;
AdjustedR_image_nii.hdr.dime.bitpix=32;
AdjustedR_image_nii.fileprefix=outpath;
AdjustedR_image_nii.hdr.dime.scl_slope = 1;
AdjustedR_image_nii.img=AdjustedR_image;
save_untouch_nii(AdjustedR_image_nii,cat(2,outpath, '/AdjustedR.nii'));

% save Residuals4D_image nifti file
Residuals4D_image_nii=Nii1;
Residuals4D_image_nii.hdr.dime.dim=[4,dimension(1),dimension(2),dimension(3),num_subjects,1,1,1];
Residuals4D_image_nii.hdr.hk.sizeof_hdr= 348;			% must be 348!
Residuals4D_image_nii.hdr.hk.data_type= '';
Residuals4D_image_nii.hdr.hk.db_name= '';
Residuals4D_image_nii.hdr.hk.extents= 0;
Residuals4D_image_nii.hdr.hk.session_error= 0;
Residuals4D_image_nii.hdr.hk.regular= 'r';
Residuals4D_image_nii.hdr.hk.dim_info= 0;
Residuals4D_image_nii.hdr.dime.datatype=16;
Residuals4D_image_nii.hdr.dime.bitpix=32;
Residuals4D_image_nii.fileprefix=outpath;
Residuals4D_image_nii.hdr.dime.scl_slope = 1;
Residuals4D_image_nii.img=Residuals4D_image;
save_untouch_nii(Residuals4D_image_nii,cat(2,outpath, '/Residuals4D.nii'));

%%%%% Bonferroni
mask_voxel = numOFnonzero_voxels(Mask_fileName);
Bonferroni=0.025./mask_voxel;
% Tthreshold
T_0001=tinv(0.9995,DFE);
T_bonferroni_005=tinv(1-Bonferroni,DFE);
T_005=tinv(0.975,DFE);

T_threshold_fileName=cat(2,outpath, '/Tthreshold');
fid=fopen(T_threshold_fileName, 'w');
fprintf(fid, cat(2, 'Degree of freedom(df): ', num2str(DFE))); fprintf(fid, '\n');fprintf(fid, '\n');
fprintf(fid, cat(2, 'Bonferroni(0.05): ', num2str(tinv(1-Bonferroni,DFE)))); fprintf(fid, '\n');
fprintf(fid, cat(2, '2-tailed(p=0.001): ', num2str(tinv(0.9995,DFE)))); fprintf(fid, '\n');
fprintf(fid, cat(2, '2-tailed(p=0.01): ', num2str(tinv(0.995,DFE)))); fprintf(fid, '\n');
fprintf(fid, cat(2, '2-tailed(p=0.05): ', num2str(tinv(0.975,DFE)))); fprintf(fid, '\n');fprintf(fid, '\n');
fprintf(fid, cat(2, '1-tailed(p=0.001): ', num2str(tinv(0.999,DFE)))); fprintf(fid, '\n');
fprintf(fid, cat(2, '1-tailed(p=0.01): ', num2str(tinv(0.99,DFE)))); fprintf(fid, '\n');
fprintf(fid, cat(2, '1-tailed(p=0.05): ', num2str(tinv(0.95,DFE)))); fprintf(fid, '\n');
fclose(fid);

diary off;