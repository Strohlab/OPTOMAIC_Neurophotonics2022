function [DC_Bilder, err_log]=DC_read_Bilder(prot_fid,analyse_animalfolder);
%%
s=['- \n']; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_read_Bilder V20220519';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%% - Line Scanning fMRI -
%% Read of images
%% DC_Bilder consist of e.g. 200 images with a size of 64x128

folder_program=cd (['../data/' analyse_animalfolder '/dicom']) 
if exist('Bilder.mat')
    load(['Bilder.mat'])
    file_generation=0;
else %DC_ this part is not tested
    file_generation=1;
    filelist=dir('MRI*');
    Info=dicominfo(filelist(1).name);
    Bilder=zeros(Info.Height,Info.Width,numel(filelist));
    laenge=double(double(numel(filelist))*double(Info.Height));
    
    for num_im=1:numel(filelist)
        num_im
        num_im=double(num_im);
        Bild = dicomread(filelist(num_im).name);
        Bilder(:,:,num_im)=Bild;
        save('Bilder','Bilder');
    end
    %_DC
end
%%
if (num2str(size(Bilder))~=0),
    err_log=0;
    s='### Images were read. \n'; fprintf(prot_fid,s);disp(s);
else
    err_log=1;
    s='### Error: Read of Images! \n'; fprintf(prot_fid,s);disp(s);
end
%% shift data to DC_Bilder for use in main script
DC_Bilder=Bilder;

%% end of function
cd(folder_program)
if (file_generation==1),
    s=['### Variable Bilder had to be generated. \n']; fprintf(prot_fid,s);disp(s)
else
    s=['### Variable Bilder loaded. \n']; fprintf(prot_fid,s);disp(s)
end
s=['### Image data size = ' num2str(size(DC_Bilder)) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)

