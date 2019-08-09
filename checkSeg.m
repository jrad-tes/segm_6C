function checkSeg(mri,seg)
% plots a given fieldtrip-mri datastructure, masked with the binary seg
    mri.mask_h = seg;
    cfg = [];
    cfg.method = 'ortho'; %'slice';
    cfg.interactive = 'yes';
    cfg.funcolorlim = 'maxabs';
    cfg.funparameter = 'mask_h';
    cfg.colorbar = 'yes';
    cfg.funcolormap = 'hot';
    cfg.funcolorlim = 'zeromax';
    ft_sourceplot(cfg,mri);
end 