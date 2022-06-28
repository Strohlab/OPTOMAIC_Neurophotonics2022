function [DC_data_collect_arr,DC_ProtImages, err_log]=DC_collect_data_stim_MRI(DC_data_long_arr,Rep_stim,datum,DC_ProtImages,err_log,prot_fid,analyse_animalfolder,animal_ident);
%%
s=['- \n']; fprintf(prot_fid,s);disp(s)
analyse_function_name_version='FUNCTION DC_collect_data_stim V20200402';
analyse_function_author='Dirk Cleppien';
s=['### (' analyse_function_name_version ' - ' analyse_function_author ')  \n']; fprintf(prot_fid,s);disp(s)

%% Global variables to use
s=['### Global Var: Rep_stim = ' num2str(Rep_stim) ' \n']; fprintf(prot_fid,s);disp(s)

%% 
% Single 2D data set depending on filter
DC_var_size_data_long=size(DC_data_long_arr(:,:,1));s=['### Var: DC_var_size_data_long = ' num2str(DC_var_size_data_long) ' \n']; fprintf(prot_fid,s);disp(s)
% Number of stimulation interval:
DC_var_length=fix(DC_var_size_data_long(2)/Rep_stim);s=['### Var: DC_var_length = ' num2str(DC_var_length) ' \n']; fprintf(prot_fid,s);disp(s)
% Spatial resolution
DC_var_dim1=size(DC_data_long_arr(:,:,1),1);s=['### Var: DC_var_dim1 = ' num2str(DC_var_dim1) ' \n']; fprintf(prot_fid,s);disp(s)
% Number of stimulations
DC_var_length=fix(DC_var_size_data_long(2)/Rep_stim);s=['### Var: DC_var_length = ' num2str(DC_var_length) ' \n']; fprintf(prot_fid,s);disp(s)

%% Output data arr (4D)
DC_data_collect_arr=zeros(DC_var_dim1,Rep_stim,DC_var_length,size(DC_data_long_arr,3));
s=['### Var DC_data_collect_arr - size = ' num2str(size(DC_data_collect_arr)) '  \n']; fprintf(prot_fid,s);disp(s)

%% Analysis figure initialisation
fig=figure('units','normalized','outerposition',[0 0 1 1]);
fig_row=size(DC_data_long_arr,3);
fig_col=2;
fig.ToolBar='none';
fig.NumberTitle='off';
fig.Name=['Subroutine Collect data stim MRI'];
fig.FileName=['Subroutine Collect data stim MRI'];

%% loop over un-/ filtered data
for z_filter=1:size(DC_data_long_arr,3),

    err_log=1;
    s=['### For-loop Var: z_filter = ' num2str(z_filter) ' \n']; fprintf(prot_fid,s);disp(s)
    
    %% filtered 2D data set
    DC_data_long=DC_data_long_arr(:,:,z_filter);
    DC_data_collect=DC_data_long;s=['### Var DC_data_long - size = ' num2str(size(DC_data_long)) '  \n']; fprintf(prot_fid,s);disp(s)
    
    %% reshape of DC_data long to stimulus, result is DC_data_collect
    dummy=DC_data_long(1:DC_var_dim1,1:Rep_stim*DC_var_length);
    size(dummy)
    DC_data_collect=reshape(dummy,DC_var_dim1,Rep_stim,DC_var_length);
    DC_var_size_data_collect=size(DC_data_collect);s=['### Size DC_data_collect = ' num2str(DC_var_size_data_collect) ' \n']; fprintf(prot_fid,s);disp(s)
    
    %% depiction of reshaping
    depicted_pixel=40;
    subplot(fig_row,fig_col,z_filter*fig_col-1)
        m_dummy=mean(dummy(depicted_pixel,1:600));
        plot(dummy(depicted_pixel,1:600)/m_dummy)
        hold on
        visualize_offset=0.5;
        plot(1:200,DC_data_collect(depicted_pixel,:,1)/m_dummy+visualize_offset,'g')
        plot(201:400,DC_data_collect(depicted_pixel,:,2)/m_dummy+visualize_offset,'r')
        plot(401:600,DC_data_collect(depicted_pixel,:,3)/m_dummy+visualize_offset,'k')
        hold off
        title({['z-filter: ' num2str(z_filter') '; Reshaping 2D to 3D: DC data collect(' num2str(size(DC_data_collect)) ')' ]; 'Offset for visualization'})
    s=['### data size = ' num2str(size(DC_data_collect)) ' \n']; fprintf(prot_fid,s);disp(s)
    s=['### 1. dim voxel position in the head; 2. dim evolution over time; 3. dim number of stimulus \n']; fprintf(prot_fid,s);disp(s)
    
    %% reshaped data into 4D arr for use in main script
    DC_data_collect_arr(:,:,:,z_filter)=DC_data_collect;
    
    %%
    mean_Bild=mean(DC_data_collect_arr(:,:,:,z_filter),3);
    subplot(fig_row,fig_col,z_filter*fig_col)
        imagesc(mean_Bild(:,:))
        colorbar
        title('Mean signal over FOV')
end % z_filter

%% save figure
image_suffix='Subroutine_DC_collect_data_stim_MRI';
old=cd(['../data/' analyse_animalfolder '/']);
DC_hg_name=[ animal_ident '_' datum '_' image_suffix '.emf' ];
saveas(fig,DC_hg_name,'emf')
cd(old)

%% end of function
s=['### Var DC_data_collect_arr - size = ' num2str(size(DC_data_collect_arr)) '  \n']; fprintf(prot_fid,s);disp(s)
err_log=0;
s=['### (' analyse_function_name_version ') - end \n']; fprintf(prot_fid,s);disp(s)