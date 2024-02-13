#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 00:51:46 2023

@author: s205272
"""
import numpy as np 

def normalize_mi_ma(x, mi, ma, clip=False, eps=1e-20, dtype=np.float32):
    if dtype is not None:
        x   = x.astype(dtype,copy=False)
        mi  = dtype(mi) if np.isscalar(mi) else mi.astype(dtype,copy=False)
        ma  = dtype(ma) if np.isscalar(ma) else ma.astype(dtype,copy=False)
        eps = dtype(eps)
    try:
        import numexpr
        x = numexpr.evaluate("(x - mi) / ( ma - mi + eps )")
    except ImportError:
        x =                   (x - mi) / ( ma - mi + eps )
    if clip:
        x = np.clip(x,0,1)
    return x

def normalize(x, pmin=2, pmax=99.8, axis=None, clip=False, eps=1e-20, dtype=np.float32):
    """Percentile-based image normalization."""

    mi = np.percentile(x,pmin,axis=axis,keepdims=True)
    ma = np.percentile(x,pmax,axis=axis,keepdims=True)
    return normalize_mi_ma(x, mi, ma, clip=clip, eps=eps, dtype=dtype)



def mkdir(directory):
    """ check if directory exists and create it through Python if it does not yet.

    Parameters
    ----------
    directory : str
        the directory path (absolute or relative) you wish to create (or check that it exists)

    Returns
    -------
        void function, no return
    """
    import os 
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    return []



if __name__=="__main__":
    
    import pylab as plt 
    import numpy as np 
    import skimage.io as skio 
    import glob
    import os 
    import scipy.ndimage as ndimage 
    from skimage.metrics import peak_signal_noise_ratio
    import skimage.filters as skfilters
    import scipy.io as spio
    
    # """
    # 1st set 
    # """
    # rootfolder = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMv2/Bo-Jui/Confocal Projection/ZebrafishBrain/Calcium/GFP/230718/Processing_55-62'
    # savefolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium'
    # imfiles = np.sort(glob.glob(os.path.join(rootfolder, 'Reg_*Sampled_t51-300.tif')))
    
    # """
    # 2nd set 
    # """
    # rootfolder = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMv2/Bo-Jui/Confocal Projection/ZebrafishBrain/Calcium/GFP/230718/Processing_39-46'
    # savefolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_39-46'
    # mkdir(savefolder)
    # imfiles = np.sort(glob.glob(os.path.join(rootfolder, 'Reg_*Sampled_t51-300.tif')))
    
    # """
    # 3rd set 
    # """
    # rootfolder = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMv2/Bo-Jui/Confocal Projection/ZebrafishBrain/Calcium/GFP/230718/Processing_47-54'
    # savefolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_47-54'
    # mkdir(savefolder)
    # imfiles = np.sort(glob.glob(os.path.join(rootfolder, 'Reg_*Sampled_t51-300.tif')))
    
    
    # """
    # 4th set 
    # """
    # rootfolder = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMv2/Bo-Jui/Confocal Projection/ZebrafishBrain/Calcium/GFP/230718/Processing_63-70'
    # savefolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_63-70'
    # mkdir(savefolder)
    # imfiles = np.sort(glob.glob(os.path.join(rootfolder, 'Reg_*Sampled_t51-300.tif')))
    
    
    """
    5th set 
    """
    rootfolder = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMv2/Bo-Jui/Confocal Projection/ZebrafishBrain/Calcium/GFP/230718/Processing_22-29'
    savefolder = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_22-29'
    mkdir(savefolder)
    imfiles = np.sort(glob.glob(os.path.join(rootfolder, 'Reg_*Sampled_t51-300.tif')))
    
    
    
    # 1.1, 2.1, 4.2, 8.5, 21.2, 42.4, 84.8, 217.2 (full) um.
    # Cell55: confocal depth= 1.1um
    # Cell56: confocal depth= 2.1um
    # Cell57: confocal depth= 4.2um
    # Cell58: confocal depth= 8.5um
    # Cell59: confocal depth= 21.2um
    # Cell60: confocal depth= 42.4um
    # Cell61: confocal depth= 84.8um
    # Cell62 : 217.2um
    
    confocal_depth = np.hstack([1.1, 2.1, 4.2, 8.5, 21.2, 42.4, 84.8, 217.2])
    mean_signals_time = [] # as a control 
    mean_signals_time_detrend = []
    
    ref_image = []
    psnr_metrics = []
    sobel_metrics = [] 
    
    for imfile_ii in np.arange(len(imfiles)):
        
        basename = os.path.split(imfiles[imfile_ii])[-1]
        im = skio.imread(imfiles[imfile_ii])
        
        im_mean = np.hstack([np.mean(imm) for imm in im])
        mean_signals_time.append(im_mean)
        
        
        # std norm per time slice. 
        im_norm = np.array([(sli-sli.mean())/ sli.std() for sli in im]) # bleach correction .
        # im_norm = np.array([sli/sli.mean() for sli in im]) # bleach correction .
        # im_norm = im.copy() 
        
        # renorm
        # im_norm = (im_norm - im_norm.min()) / (im_norm.max()-im_norm.min()) # this takes out. 
        # percent norm
        im_norm = np.array([normalize(imm, pmin=2, pmax=99.8, clip=True, eps=1e-20, dtype=np.float32) for imm in im_norm])
        
        im_norm = np.uint8(255*im_norm)
        
        im_norm_sobel = np.array([np.mean(np.abs(skfilters.sobel(sli))) for sli in im_norm])
        sobel_metrics.append(im_norm_sobel)
        
        
        """
        Save this version. 
        """
        skio.imsave(os.path.join(savefolder,
                                 'norm_'+basename), 
                    im_norm)
        
        # background subtraction 
        if imfile_ii == 0:
            # ref_image = np.max(im, axis=0)
            ref_image = np.max(im_norm, axis=0)
            # ref_image = im_norm[0]
            
        # now compute. 
        psnr_im = np.hstack(peak_signal_noise_ratio(ref_image, im_norm_tt) for im_norm_tt in im_norm)
        # psnr_im = np.hstack(peak_signal_noise_ratio(ref_image, im_tt) for im_tt in im)
        psnr_metrics.append(psnr_im)
        
        # very clear offset 
        mean_signals_time_detrend.append((im_mean-np.mean(im_mean))/ np.std(im_mean))
        
    mean_signals_time = np.vstack(mean_signals_time)
    mean_signals_time_detrend = np.vstack(mean_signals_time_detrend)
    
    
    psnr_metrics = np.vstack(psnr_metrics)
    sobel_metrics = np.vstack(sobel_metrics)
    
    
    depth_points = np.vstack([np.ones(psnr_metrics.shape[1])*dd for dd in confocal_depth])
    
    
    # do a plot of PSNR. 
    plt.figure(figsize=(10,10))
    plt.boxplot(psnr_metrics.T)
    plt.ylabel('PSNR', fontsize=20)
    plt.xlabel('confocal depth [um]', fontsize=20)
    plt.xticks(np.arange(len(confocal_depth))+1, confocal_depth, fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    
    # do a plot of sobel. 
    plt.figure(figsize=(10,10))
    plt.boxplot(sobel_metrics.T)
    plt.ylabel('Sobel Magnitude', fontsize=20)
    plt.xlabel('confocal depth [um]', fontsize=20)
    plt.xticks(np.arange(len(confocal_depth))+1, confocal_depth, fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    
    mean_psnrs = np.nanmean(psnr_metrics, axis=1)
    mean_sobels = np.nanmean(sobel_metrics, axis=1)
    
    # plot these points on a scatter plot. 
    plt.figure(figsize=(10,10))
    plt.scatter(depth_points.ravel(), psnr_metrics.ravel(), s=50, c='k')
    plt.plot(depth_points[:,0], mean_psnrs, 'k-', lw=3)
    plt.ylabel('PSNR', fontsize=20, fontname='Liberation Sans')
    plt.xlabel('confocal depth [um]', fontsize=20, fontname='Liberation Sans')
    plt.xticks(fontsize=18, fontname='Liberation Sans')
    plt.yticks(fontsize=18, fontname='Liberation Sans')
    plt.ylim([5,35])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    plt.savefig(os.path.join(savefolder,
                             'mean_psnr_SNR_vs_confocal_depth.pdf'), 
                 dpi=300, bbox_inches='tight')
    plt.show()
    
    
    
    # plot these points on a scatter plot. 
    plt.figure(figsize=(10,10))
    plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
    plt.plot(depth_points[:,0], mean_sobels, 'k-', lw=3)
    plt.ylabel('Sobel Magnitude', fontsize=20, fontname='Liberation Sans')
    plt.xlabel('confocal depth [um]', fontsize=20, fontname='Liberation Sans')
    plt.xticks(fontsize=18, fontname='Liberation Sans')
    plt.yticks(fontsize=18, fontname='Liberation Sans')
    # plt.ylim([0.005,0.025])
    plt.ylim([0.005,0.035])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    plt.savefig(os.path.join(savefolder,
                             'mean_sobel_SNR_vs_confocal_depth.pdf'), 
                 dpi=300, bbox_inches='tight')
    plt.show()
    
    
    """
    Save the metrics so we can replot 
    """
    savedict = {'depth_points': depth_points,
                'sobel_metrics': sobel_metrics,
                'psnr_metrics': psnr_metrics,
                'mean_psnrs': mean_psnrs,
                'mean_sobels':mean_sobels,
                'mean_signals_time': mean_signals_time,
                'mean_signals_time_detrend': mean_signals_time_detrend,
                'confocal_depth': confocal_depth}
    
    spio.savemat(os.path.join(savefolder, 
                              'SNR_statistics.mat'), 
                 savedict)
    
    
    
    
    