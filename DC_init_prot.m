function [datum,prot_fid,prot_filename,analyse_animalfolder]=DC_init_prot;

analyse_function_name_version='FUNCTION init_prot V20220519';
analyse_function_author='Dirk Cleppien';
%% Function to initialize protocol file
prot_fid=0;

%% infos to create protocol filename
%datum=datestr(now,'yymmdd_HH-MM-SS');
datum=input('Input of the actual date [yyyymmdd]:','s');
analyse_initials=input('Input of initials of analyzing person:','s');

%% choose of analyzing animal
dname=uigetdir(['../'],'For Documentation: Please select the animal folder to analyze (with a /dicom subfolder containing measured dicom files):');

%% change into analyzing folder
folder_program=cd (dname);

%% creating a new protocol file
d=findstr('\',dname);
analyse_animalfolder=dname(d(end)+1:end)
prot_filename=['OPTOMAIC_Analyse_' analyse_animalfolder '_' analyse_initials '_' datum]
prot_fid=fopen(prot_filename,'w');

%% change back in matlab folder
cd(folder_program);

%% general information of analysis written in protocol file
fprintf(prot_fid,'OPTOMAIC \n');
fprintf(prot_fid,'- Analysis of a single experiment - \n');
fprintf(prot_fid,'General protocol information \n');
fprintf(prot_fid,' \n');
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s);
s=['### Date: ' datum ' \n']; fprintf(prot_fid,s);disp(s);
s=['### Initials of analyzing person: ' analyse_initials ' \n']; fprintf(prot_fid,s);disp(s);
s=['### Analyzing folder: ../data/' analyse_animalfolder ' \n']; fprintf(prot_fid,s);disp(s);
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s);

%% end of file




