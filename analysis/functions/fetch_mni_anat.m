%Fetches mni anatomy so can plot electrodes and other in MNI space.

function mni=fetch_mni_anat(root_dir)
    mni.anat.fs_outputdir=[root_dir '/mni/mni-icbm152-t1-tal-nlin-asym-09b-hires_freesurfer-output'];
    mni.anat.fs_t1=ft_read_mri([mni.anat.fs_outputdir '/mri/orig.mgz']);
    mni.anat.pial=ft_read_headshape({[mni.anat.fs_outputdir '/surf/lh.pial.T1'],[mni.anat.fs_outputdir '/surf/rh.pial.T1']});
    mni.anat.l_pial=ft_read_headshape({[mni.anat.fs_outputdir '/surf/lh.pial.T1']});
    mni.anat.r_pial=ft_read_headshape({[mni.anat.fs_outputdir '/surf/rh.pial.T1']});
    
    %transform to ras space:
    transform=mni.anat.fs_t1.transform*inv(mni.anat.fs_t1.hdr.tkrvox2ras);
    temp_pos=mni.anat.pial.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.pial.pos=transform*temp_pos;
    mni.anat.pial.pos=mni.anat.pial.pos';
    mni.anat.pial.pos(:,4)=[];
    
    %transform to ras space:
    temp_pos=mni.anat.l_pial.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.l_pial.pos=transform*temp_pos;
    mni.anat.l_pial.pos=mni.anat.l_pial.pos';
    mni.anat.l_pial.pos(:,4)=[];
    
    %transform to ras space:
    temp_pos=mni.anat.r_pial.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.r_pial.pos=transform*temp_pos;
    mni.anat.r_pial.pos=mni.anat.r_pial.pos';
    mni.anat.r_pial.pos(:,4)=[];
    
    [mni.anat.mesh_lh,mni.anat.mesh_rh]=create_anat_model(mni.anat.fs_outputdir,{'Left-Hippocampus','Right-Hippocampus'});
    
    %transform to ras space:
    temp_pos=mni.anat.mesh_lh.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.mesh_lh.pos=transform*temp_pos;
    mni.anat.mesh_lh.pos=mni.anat.mesh_lh.pos';
    mni.anat.mesh_lh.pos(:,4)=[];
    
    %transform to ras space:
    temp_pos=mni.anat.mesh_rh.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.mesh_rh.pos=transform*temp_pos;
    mni.anat.mesh_rh.pos=mni.anat.mesh_rh.pos';
    mni.anat.mesh_rh.pos(:,4)=[];
    
    [mni.anat.mesh_la,mni.anat.mesh_ra]=create_anat_model(mni.anat.fs_outputdir,{'Left-Amygdala','Right-Amygdala'});
    
    %transform to ras space:
    temp_pos=mni.anat.mesh_la.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.mesh_la.pos=transform*temp_pos;
    mni.anat.mesh_la.pos=mni.anat.mesh_la.pos';
    mni.anat.mesh_la.pos(:,4)=[];
    
    %transform to ras space:
    temp_pos=mni.anat.mesh_ra.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    mni.anat.mesh_ra.pos=transform*temp_pos;
    mni.anat.mesh_ra.pos=mni.anat.mesh_ra.pos';
    mni.anat.mesh_ra.pos(:,4)=[];
    mni.anat.electrodes_info=importdata([root_dir '/mni/mni-icbm152-t1-tal-nlin-asym-09b-hires_electrode-models/sub-mni-icbm152-t1-tal-nlin-asym-09b-hires_electrodes_info.mat'],'electrodes_info');
    
end