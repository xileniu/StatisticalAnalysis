function afni_save_png(source_path,underlay_name,overlay_name,colorBar,xyz_num,picNum,range_num,threshold_num,outpath,prefix)
% source_path: paths for the images of underlay and overlay
% underlay_name: filename of the underlay. e.g. underlay_name='Sym_Template_6_MNI.nii'
% overlay_name: filename of the overlay
% colorBar: name of the colorbar in afni
% xyz_num: coordinates. e.g. xyz_num='0 18 6';
% picNum: specify rows, columns, and intervals of the slices
% range_num: max value of the colorbar. e.g. range_num='5'
% threshold_num: threshold value. e.g. threshold_num='3400 1'/'0000 0'

cd=['cd ' source_path ';'];
window=['afni -com ''OPEN_WINDOW A geom=+2+79'' -com ''OPEN_PANEL A.Define_Overlay''  -com ''OPEN_WINDOW A.axialimage geom=346x207+5+490 ifrac=0.8 mont=' num2str(picNum(1)) 'x' num2str(picNum(2)) ':' num2str(picNum(3)) ':0:none opacity=9''  '];
underlay=[' -com ''SWITCH_UNDERLAY A.' underlay_name ''' '];
overlay=[' -com ''SWITCH_OVERLAY A.' overlay_name ''' '];
xhairs=' -com ''SET_XHAIRS A.OFF'' ';
no_understand=[' -com ''SET_FUNC_VISIBLE A.+'' -com ''SET_PBAR_ALL A.-99 1.000000 ' colorBar ''' -com ''SET_FUNC_RESAM A.NN.NN'' '];
range=[' -com ''SET_FUNC_RANGE A.' range_num ''' '];
threshold = ['-com ''SET_THRESHOLD A.' threshold_num ''' '];
xyz=[' -com ''SET_DICOM_XYZ A ' xyz_num ''' '];
savejpg=[' -com ''SAVE_PNG A.axialimage ' outpath '/' prefix ''   ''' '];
quit=' -com ''QUIT''';
cmd = [cd window underlay overlay xhairs no_understand range threshold xyz  savejpg quit];
system(cmd);
% make sure the image has been generated
while exist([outpath '/' prefix '.png'],'dir')==0
    if exist([outpath '/' prefix '.png'],'dir')==1
        break
    end
end