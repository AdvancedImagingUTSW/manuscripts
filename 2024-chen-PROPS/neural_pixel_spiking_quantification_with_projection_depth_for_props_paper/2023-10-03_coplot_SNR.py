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
    
    
    all_folders = ['/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium', 
                   '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_39-46',
                   '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_47-54',
                   '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_63-70']
                   # '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_22-29']
    
    # names = ['ctrl', 
    #          '39-46',
    #          '47-54',
    #          '63-70',
    #          '22-29']
    
    
    ##Could you change "ctrl", "39-46", .etc to fish1 group1, fish1 group2, ...and I think there's one is fish2 group1.
    # names = ['Fish 1 - ROI 1', 
    #          'Fish 1 - ROI 2',
    #          'Fish 1 - ROI 3',
    #          'Fish 1 - ROI 4',
    #          'Fish 2 - ROI 1']
    
    names = ['ROI 1', 
             'ROI 2',
             'ROI 3',
             'ROI 4']#,
             # 'Fish 2 - ROI 1']
    
    confocal_depth = np.hstack([1.1, 2.1, 4.2, 8.5, 21.2, 42.4, 84.8, 217.2])
    
    # """
    # Save the metrics so we can replot 
    # """
    # savedict = {'depth_points': depth_points,
    #             'sobel_metrics': sobel_metrics,
    #             'psnr_metrics': psnr_metrics,
    #             'mean_psnrs': mean_psnrs,
    #             'mean_sobels':mean_sobels,
    #             'mean_signals_time': mean_signals_time,
    #             'mean_signals_time_detrend': mean_signals_time_detrend,
    #             'confocal_depth': confocal_depth}
    
    # spio.savemat(os.path.join(savefolder, 
    #                           'SNR_statistics.mat'), 
    #              savedict)
    
    all_psnrs = []
    all_sobels = []
    
    for folder in all_folders:
        
        stats = spio.loadmat(os.path.join(folder, 
                                  'SNR_statistics.mat'))
        
        psnr = stats['mean_psnrs']
        sobel = stats['mean_sobels']
        
        all_psnrs.append(psnr)
        all_sobels.append(sobel)
        
    all_psnrs = np.vstack(all_psnrs)
    all_sobels = np.vstack(all_sobels)
    
    # plot these points on a scatter plot. 
    plt.figure(figsize=(7,8))
    
    for iii in np.arange(len(all_psnrs)):
    # plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
        plt.plot(confocal_depth, all_psnrs[iii], 'o-', lw=3, label=names[iii], ms=10)
    plt.xscale('log')
    plt.ylabel('PSNR', fontsize=28, fontname='Liberation Sans')
    plt.xlabel(r'Projection depth [um]', fontsize=28, fontname='Liberation Sans')
    plt.xticks(fontsize=24, fontname='Liberation Sans')
    plt.yticks(fontsize=24, fontname='Liberation Sans')
    plt.legend(fontsize="18")
    # plt.ylim([0.005,0.025])
    plt.ylim([0.005,0.035])
    plt.ylim([10,35])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    # plt.legend(bbox_to_anchor=(1.15, 1.0), loc='upper left', fontsize="18")
    plt.tight_layout()
    plt.savefig(os.path.join('.',
                              'mean_SNR_vs_confocal_depth.pdf'), 
                  dpi=300, bbox_inches='tight')
    plt.show()
    
    
    # plot these points on a scatter plot. 
    plt.figure(figsize=(7,8))
    
    for iii in np.arange(len(all_psnrs)):
    # plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
        plt.plot(confocal_depth, all_psnrs[iii], 'o-', lw=3, label=names[iii], ms=10)
    # plt.xscale('log')
    plt.ylabel('PSNR', fontsize=28, fontname='Liberation Sans')
    plt.xlabel(r'Projection depth [um]', fontsize=28, fontname='Liberation Sans')
    plt.xticks(fontsize=24, fontname='Liberation Sans')
    plt.yticks(fontsize=24, fontname='Liberation Sans')
    plt.legend(fontsize="18")
    # plt.ylim([0.005,0.025])
    plt.ylim([0.005,0.035])
    plt.ylim([10,35])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    # plt.legend(bbox_to_anchor=(1.15, 1.0), loc='upper left', fontsize="18")
    plt.tight_layout()
    plt.savefig(os.path.join('.',
                              'mean_SNR_vs_confocal_depth_linear.pdf'), 
                  dpi=300, bbox_inches='tight')
    plt.show()
    
    
    