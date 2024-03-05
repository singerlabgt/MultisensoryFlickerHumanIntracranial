
function [varargout]=create_anat_model(fs_outputdir,anat_names_array)
    
    %get atlas and set some parameters for mask:
    atlas = ft_read_atlas([fs_outputdir '/mri/aparc+aseg.mgz']);
    atlas.coordsys = 'ras';
    cfg1            = [];
    cfg1.inputcoord = 'ras';
    cfg1.atlas      = atlas;
    
    %set some parameters to make the mesh:
    seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    cfg2             = [];
    cfg2.method      = 'iso2mesh';
    cfg2.radbound    = 2;
    cfg2.maxsurf     = 0;
    cfg2.tissue      = 'brain';
    cfg2.numvertices = 1000;
    cfg2.smooth      = 3;
    cfg2.spmversion  = 'spm12';
    
    %build structure:
    for i=1:length(anat_names_array)
        %make mask for this structure
        cfg1.roi        = anat_names_array(i);
        mask = ft_volumelookup(cfg1, atlas);
        
        %make mesh for this structure:
        seg.brain = mask;
        mesh{i} = ft_prepare_mesh(cfg2, seg);
        
        %adjust to ras coordinate system:
        mesh{i}.pos(:,1)=mesh{i}.pos(:,1)-atlas.hdr.c_r;
        mesh{i}.pos(:,2)=mesh{i}.pos(:,2)-atlas.hdr.c_a;
        mesh{i}.pos(:,3)=mesh{i}.pos(:,3)-atlas.hdr.c_s;
    end
    
    varargout=mesh;
end