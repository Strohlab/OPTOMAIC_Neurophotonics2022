function [DC_data_long,mean_DC_spatial_profile, DC_ProtImages, err_log]=DC_reshape_Bilder(DC_Bilder,DC_sync_StimMRI,DC_stim_offset,cutoff_stst,Rep_stim,datum,DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);
%% 
s=['- \n']; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_reshape_Bilder V20220519';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%% Global variable to use
s=['### Global Var: Rep_Stim = ' num2str(Rep_stim) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### Global Var: cutoff_stst = ' num2str(cutoff_stst) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### Global Var: DC_sync_StimMRI = ' num2str(DC_sync_StimMRI) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### Global Var: DC_stim_offset = ' num2str(DC_stim_offset) ' \n']; fprintf(prot_fid,s);disp(s)

%% just to depict an voxel example
depicted_pixel=40;
s=['### local Var: depicted_pixel = ' num2str(depicted_pixel) ' \n']; fprintf(prot_fid,s);disp(s)

%%
err_log=1;

%% reshape of images from 3d array to 2d array
dim=size(DC_Bilder);s=['### Var dim: size = ' num2str((dim)) ' \n']; fprintf(prot_fid,s);disp(s);
                    s=['### Var dim: 1. Voxel per Line; 2. Lines per Image; 3. Number of Images: (Lines per Image) x (Number of Images) = (time course) \n']; fprintf(prot_fid,s);disp(s);

DC_data_long=reshape(DC_Bilder,dim(1),dim(2)*dim(3));s=['### Var DC_data_long: size = ' num2str(size(DC_data_long)) ' \n']; fprintf(prot_fid,s);disp(s);

%% analysis figure
fig=figure('units','normalized','outerposition',[0 0 1 1]);
fig_row=2;
fig_col=3;
fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=['Subroutine Reshape Bilder'];
fig.FileName=['Subroutine Reshape Bilder'];

%% image DC_Data_long
subplot(2,3,1)
    imagesc(DC_data_long(:,:));
    title('Raw DC-data-long')
s=['### data size = ' num2str(size(DC_data_long)) ' \n']; fprintf(prot_fid,s);disp(s)
s=['### 1. dim voxel position in the head; 2. dim evolution over time \n']; fprintf(prot_fid,s);disp(s)

%% cut off signal propagation into steady state
s=['### cut off signal propagation into signal steady state \n']; fprintf(prot_fid,s);disp(s);
DC_data_long_stst=DC_data_long;s=['### Var DC_data_long_stst: size = ' num2str(size(DC_data_long_stst)) ' \n']; fprintf(prot_fid,s);disp(s);

%% cut off propagation into signal stst
DC_data_long_stst(:,1:cutoff_stst)=[];s=['### Var DC_data_long_stst: size = ' num2str(size(DC_data_long_stst)) ' \n']; fprintf(prot_fid,s);disp(s);
subplot(2,3,4)
    imagesc(DC_data_long_stst(:,:))
    title('Signal without Steady State settling')

subplot(2,3,2)
    x_data=1:400;
    plot(x_data,DC_data_long(depicted_pixel,1:400));
    hold on
    plot(x_data(cutoff_stst+1:end),DC_data_long_stst(depicted_pixel,1:200)+500)
    hold off
    title({'Raw vs cutted signal';['Offset = ' num2str(cutoff_stst)]})
    
%% shift for use in main script
DC_data_long=DC_data_long_stst;s=['### Var DC_data_long: size = ' num2str(size(DC_data_long)) ' \n']; fprintf(prot_fid,s);disp(s);

%% spatial signal intensity 
size(DC_data_long)
mean_DC_spatial_profile=mean(DC_data_long(:,10:100),2);
size(mean_DC_spatial_profile)

%% depiction
x_spatial=1:size(mean_DC_spatial_profile,1);
subplot(2,3,3)
    plot(mean_DC_spatial_profile)
    DC_noise=mean_DC_spatial_profile;
    s_DC_noise=size(DC_noise,1)
    DC_noise(5:s_DC_noise-5)=0;
    thres_Signal=zeros(s_DC_noise,1)
    thres_Signal(:)=mean(DC_noise(DC_noise>0))+2*std(DC_noise(DC_noise>0))
    hold on
    plot(DC_noise)
    plot(thres_Signal)
    d=mean_DC_spatial_profile.*(mean_DC_spatial_profile>thres_Signal)
    box=d>0;
    plot(d)
    plot(box.*mean_DC_spatial_profile)
    xcount=1:64;
    pixel=(xcount(d>0))
    xlim([1 size(mean_DC_spatial_profile,1)])
    [py,px]=max(mean_DC_spatial_profile);
    title({['Depiction of signal profile over voxel'],...
        ['Pixel ' num2str(pixel(1:4)) ' - ' num2str(pixel(end-3:end))]...
        ['Voxel with max. Signal = ' num2str(px)]});
    hold off

spatial_SI=x_spatial(box(:)==1);
spatial_SI_min=min(spatial_SI);;s=['### Min. spatial SI voxel = ' num2str(spatial_SI_min) ' \n']; fprintf(prot_fid,s);disp(s);
spatial_SI_max=max(spatial_SI);;s=['### Max. spatial SI voxel = ' num2str(spatial_SI_max) ' \n']; fprintf(prot_fid,s);disp(s);

%% sync between MRI and stimulus pulse
s=['### Synchronization between MRI and stimulation pulse: \n']; fprintf(prot_fid,s);disp(s);
box(:)=1;
DC_data_long_sync=DC_data_long;s=['### Var DC_data_long_sync: size = ' num2str(size(DC_data_long_sync)) ' \n']; fprintf(prot_fid,s);disp(s);

% cut off for sync between MRI and stimulus pulse - baseline before
d_sync=DC_sync_StimMRI-DC_stim_offset;s=['### First d_sync= ' num2str(d_sync) ' are removed \n']; fprintf(prot_fid,s);disp(s);
% arbitrary offset before stimulus has to be included:
DC_data_long_sync(:,1:d_sync)=[];s=['### Var DC_data_long: size = ' num2str(size(DC_data_long_sync)) ' \n']; fprintf(prot_fid,s);disp(s);

%%
subplot(2,3,5)
    x_data=1:400;
    plot(x_data,DC_data_long(depicted_pixel,1:400));
    hold on
    plot(x_data(d_sync+1:end),DC_data_long_sync(depicted_pixel,1:400-d_sync)+800)
    hold off
    title({'Cutted vs synchronized signal';['Sync offset = ' num2str(d_sync)];'first stim +20'})
DC_data_long=DC_data_long_sync;s=['### Var DC_data_long: size = ' num2str(size(DC_data_long)) ' \n']; fprintf(prot_fid,s);disp(s);

subplot(2,3,6)
     imagesc(DC_data_long(:,:))
     title({'Synchronized signal';['Size = ' num2str(size(DC_data_long))]})

%% save analysis figure
image_suffix='Subroutine_DC_Reshape_Bilder';
old=cd(['../data/' analyse_animalfolder '/']);
DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)

%% end of function
err_log=0;
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)
