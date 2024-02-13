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


# potentially extend this to handle anistropic! 
def smooth_vol(vol_binary, ds=4, smooth=5):
    
    from skimage.filters import gaussian
    from scipy.ndimage import gaussian_filter
    import skimage.transform as sktform
    import numpy as np 
    
    small = sktform.resize(vol_binary, np.array(vol_binary.shape)//ds, preserve_range=True, mode='reflect')
    small = gaussian_filter(small, sigma=smooth)
    
    return sktform.resize(small, np.array(vol_binary.shape), preserve_range=True)


def xcorr(x, y=None, norm=True, eps=1e-12, mode='full'):
	r""" Computes the discrete crosscorrelation of two read 1D signals as defined by

	.. math::
		c_k = \sum_n x_{n+k} \cdot y_n

	If norm=True, :math:`\hat{x}=\frac{x-\mu}{\sigma}` and :math:`\hat{y}=\frac{y-\mu}{\sigma}` is normalised and the zero-normalised autocorrelation is computed,  
	
	.. math::
		c_k = \frac{1}{T}\sum_n \hat{x}_{n+k} \cdot \hat{y}_n

	where :math:`T` is the length of the signal :math:`x`

	Parameters
	----------
	x : 1D numpy array
		The input 1D signal
	y : 1D numpy array
		The optional second 1D signal. If None, the autocorrelation of x is computed.   
	norm : bool
		If true, the normalized autocorrelation is computed such that all values are in the range [-1.,1.]
	eps :  scalar
		small constant to prevent zero division when norm=True

	Returns
	-------
	result : 1D numpy array
		the 1-sided autocorrelation if y=None, or the full cross-correlation otherwise

	Notes
	-----
	The definition of correlation above is not unique and sometimes correlation
	may be defined differently. Another common definition is:

	.. math:: 
		c'_k = \sum_n x_{n} \cdot {y_{n+k}}
	
	which is related to :math:`c_k` by :math:`c'_k = c_{-k}`.
	"""
	import numpy as np 

	if norm: 
		a = (x - np.nanmean(x)) / (np.nanstd(x) * len(x) + eps)
		if y is not None:
			b = (y - np.nanmean(y)) / (np.nanstd(y) + eps)
		else:
			b = a.copy()
	else:
		a = x.copy()
		if y is not None:
			b = y.copy()
		b = x.copy()
	result = np.correlate(a, b, mode=mode) # this is not normalized!. 

	if y is None: 
		# return the 1-sided autocorrelation. 
		result = result[result.size // 2:]

	return result

def xcorr_timeseries_set_1d(timeseries_array1, timeseries_array2=None, norm=True, eps=1e-12, stack_final=False):
	r""" Computes the discrete crosscorrelation of two read 1D signals as defined by

	.. math::
		c_k = \sum_n x_{n+k} \cdot y_n

	If norm=True, :math:`\hat{x}=\frac{x-\mu}{\sigma}` and :math:`\hat{y}=\frac{y-\mu}{\sigma}` is normalised and the zero-normalised autocorrelation is computed,  
	
	.. math::
		c_k = \frac{1}{T}\sum_n \hat{x}_{n+k} \cdot \hat{y}_n

	where :math:`T` is the length of the signal :math:`x`

	given two arrays or lists of 1D signals. The signals in the individual arrays need not have the same temporal length. If they do, by setting stack_final=True, the result can be returned as a numpy array else will be returned as a list  

	Parameters
	----------
	timeseries_array1 : array_like of 1D signals 
		A input list of 1D signals 
	timeseries_array2 : array_like of 1D signals 
		An optional second 1D signal set. If None, the autocorrelation of each timeseries in timeseries_array1 is computed.   
	norm : bool
		If true, the normalized autocorrelation is computed such that all values are in the range [-1.,1.]
	eps :  scalar
		small constant to prevent zero division when norm=True
	stack_final : bool
		if timeseries_array1 or timeseries_array2 are numpy arrays with individual signals within of equal temoporal length, setting this flag to True, will return a numpy array else returns a list of cross-correlation curves

	Returns
	-------
	xcorr_out : array_list of 1D numpy array
		the 1-sided autocorrelation if y=None, or the full cross-correlation otherwise

	Notes
	-----
	The definition of correlation above is not unique and sometimes correlation
	may be defined differently. Another common definition is:

	.. math:: 
		c'_k = \sum_n x_{n} \cdot {y_{n+k}}
	
	which is related to :math:`c_k` by :math:`c'_k = c_{-k}`.
	"""
	import numpy as np 
	compute_xcorr=True 
	if timeseries_array2 is None:
		timeseries_array2 = timeseries_array1.copy() # create a copy.
		compute_xcorr=False

	xcorr_out = []
	for ii in np.arange(len(timeseries_array1)):
		timeseries_ii_1 = timeseries_array1[ii].copy()
		timeseries_ii_2 = timeseries_array2[ii].copy()
		if compute_xcorr:
			xcorr_timeseries_ii = xcorr(timeseries_ii_1, y=timeseries_ii_2, norm=norm, eps=eps)
		else:
			xcorr_timeseries_ii = xcorr(timeseries_ii_1, y=None, norm=norm, eps=eps)
		xcorr_out.append(xcorr_timeseries_ii)

	if stack_final:
		xcorr_out = np.vstack(xcorr_out)
		return xcorr_out
	else:
		return xcorr_out


def xcorr_timeseries_set_1d(timeseries_array1, timeseries_array2=None, norm=True, eps=1e-12, stack_final=False):
	r""" Computes the discrete crosscorrelation of two read 1D signals as defined by

	.. math::
		c_k = \sum_n x_{n+k} \cdot y_n

	If norm=True, :math:`\hat{x}=\frac{x-\mu}{\sigma}` and :math:`\hat{y}=\frac{y-\mu}{\sigma}` is normalised and the zero-normalised autocorrelation is computed,  
	
	.. math::
		c_k = \frac{1}{T}\sum_n \hat{x}_{n+k} \cdot \hat{y}_n

	where :math:`T` is the length of the signal :math:`x`

	given two arrays or lists of 1D signals. The signals in the individual arrays need not have the same temporal length. If they do, by setting stack_final=True, the result can be returned as a numpy array else will be returned as a list  

	Parameters
	----------
	timeseries_array1 : array_like of 1D signals 
		A input list of 1D signals 
	timeseries_array2 : array_like of 1D signals 
		An optional second 1D signal set. If None, the autocorrelation of each timeseries in timeseries_array1 is computed.   
	norm : bool
		If true, the normalized autocorrelation is computed such that all values are in the range [-1.,1.]
	eps :  scalar
		small constant to prevent zero division when norm=True
	stack_final : bool
		if timeseries_array1 or timeseries_array2 are numpy arrays with individual signals within of equal temoporal length, setting this flag to True, will return a numpy array else returns a list of cross-correlation curves

	Returns
	-------
	xcorr_out : array_list of 1D numpy array
		the 1-sided autocorrelation if y=None, or the full cross-correlation otherwise

	Notes
	-----
	The definition of correlation above is not unique and sometimes correlation
	may be defined differently. Another common definition is:

	.. math:: 
		c'_k = \sum_n x_{n} \cdot {y_{n+k}}
	
	which is related to :math:`c_k` by :math:`c'_k = c_{-k}`.
	"""
	import numpy as np 
	compute_xcorr=True 
	if timeseries_array2 is None:
		timeseries_array2 = timeseries_array1.copy() # create a copy.
		compute_xcorr=False

	xcorr_out = []
	for ii in np.arange(len(timeseries_array1)):
		timeseries_ii_1 = timeseries_array1[ii].copy()
		timeseries_ii_2 = timeseries_array2[ii].copy()
		if compute_xcorr:
			xcorr_timeseries_ii = xcorr(timeseries_ii_1, y=timeseries_ii_2, norm=norm, eps=eps)
		else:
			xcorr_timeseries_ii = xcorr(timeseries_ii_1, y=None, norm=norm, eps=eps)
		xcorr_out.append(xcorr_timeseries_ii)

	if stack_final:
		xcorr_out = np.vstack(xcorr_out)
		return xcorr_out
	else:
		return xcorr_out


def _ma_average(sig, winsize=3, mode='reflect', avg_func=np.nanmean):
    
    sig_ = np.pad(sig, [winsize//2, winsize//2], mode=mode)
    sig_out = []
    for ii in np.arange(len(sig)):
        data = sig_[ii:ii+winsize]
        sig_out.append(avg_func(data))
    return np.hstack(sig_out)


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
    from tqdm import tqdm 
    from scipy.signal import find_peaks
    import skimage.morphology as skmorph
    import skimage.segmentation as sksegmentation 
    import skimage.measure as skmeasure 
    import scipy.io as spio
    
    all_folders = ['/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium', 
                   '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_39-46',
                   '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_47-54',
                    '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_63-70']#,
                   # '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/fzhou/Bo-Jui/OPMv2/Zebrafish_Brain_Calcium/230718/Processing_22-29']
    
    # names = ['ctrl', 
    #          '39-46',
    #          '47-54',
    #          '63-70',
    #          '22-29']
    # names = ['Fish 1 - ROI 1', 
    #          'Fish 1 - ROI 2',
    #          'Fish 1 - ROI 3',
    #          'Fish 1 - ROI 4',
    #          'Fish 2 - ROI 1']
    names = ['ROI 1', 
             'ROI 2',
             'ROI 3',
             'ROI 4']
             # 'Fish 2 - ROI 1']
    
    confocal_depth = np.hstack([1.1, 2.1, 4.2, 8.5, 21.2, 42.4, 84.8, 217.2])
    
    
    # import scipy.io as spio
    # spio.savemat(os.path.join(savefolder,
    #                           'depth_calcium_firing_domains.mat'),
    #              {'area_ratio': np.hstack(total_area)/np.hstack(total_tissue_area),
    #               'confocal_depth': confocal_depth,
    #               'total_area': total_area,
    #               'total_tissue': total_tissue_area})
    
    all_areas = []
    # all_sobels = []
    
    for folder in all_folders:
        
        stats = spio.loadmat(os.path.join(folder, 
                                  'depth_calcium_firing_domains.mat'))
        
        area = np.hstack(stats['total_area'].ravel()) / np.hstack(stats['total_tissue'].ravel()) 
        all_areas.append(area)
        
    all_areas = np.vstack(all_areas)
    
    mean_area = np.nanmean(all_areas, axis=0)
    # all_sobels = np.vstack(all_sobels)
    
    # plot these points on a scatter plot. 
    plt.figure(figsize=(10,10))
    
    for iii in np.arange(len(all_areas)):
    # plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
        plt.plot(confocal_depth, all_areas[iii], 'o-', lw=3, ms=15, label=names[iii])
        
    plt.plot(confocal_depth, mean_area, 'k-o', lw=10, ms=20, label='Mean', alpha=1)
    plt.xscale("log")
    plt.ylabel('Signal Area / Tissue Area', fontsize=36, fontname='Liberation Sans')
    plt.xlabel('Projection depth [um]', fontsize=36, fontname='Liberation Sans')
    plt.xticks(fontsize=32, fontname='Liberation Sans')
    plt.yticks(fontsize=32, fontname='Liberation Sans')
    plt.legend(fontsize="24")
    # plt.ylim([0.005,0.025])
    # plt.ylim([0.005,0.035])
    plt.ylim([0,0.25])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    plt.tight_layout()
    plt.savefig(os.path.join('.',
                              'mean_signal-to-tissue_area_vs_confocal_depth.pdf'), 
                  dpi=300, bbox_inches='tight')
    plt.show()

    
    
    
    
    median_area = np.nanmedian(all_areas, axis=0)
    
    plt.figure(figsize=(10,10))
    
    for iii in np.arange(len(all_areas)):
    # plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
        plt.plot(confocal_depth, all_areas[iii], 'o-', lw=3, ms=15, label=names[iii])
        
    plt.plot(confocal_depth, median_area, 'k-o', lw=10, ms=20, label='Median', alpha=1)
    plt.xscale("log")
    plt.ylabel('Signal Area / Tissue Area', fontsize=36, fontname='Liberation Sans')
    plt.xlabel('Projection depth [um]', fontsize=36, fontname='Liberation Sans')
    plt.xticks(fontsize=32, fontname='Liberation Sans')
    plt.yticks(fontsize=32, fontname='Liberation Sans')
    plt.legend(fontsize="24")
    # plt.ylim([0.005,0.025])
    # plt.ylim([0.005,0.035])
    # plt.ylim([10,35])
    plt.ylim([0,0.25])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    plt.tight_layout()
    plt.savefig(os.path.join('.',
                              'median_signal-to-tissue_area_vs_confocal_depth.pdf'), 
                  dpi=300, bbox_inches='tight')
    plt.show()
    
    
    """
    linear scales
    """
    # plot these points on a scatter plot. 
    plt.figure(figsize=(10,10))
    
    for iii in np.arange(len(all_areas)):
    # plt.scatter(depth_points.ravel(), sobel_metrics.ravel(), s=10, c='k')
        plt.plot(confocal_depth, all_areas[iii], 'o-', lw=3, ms=15, label=names[iii])
        
    plt.plot(confocal_depth, mean_area, 'k-o', lw=10, ms=20, label='Mean', alpha=1)
    # plt.xscale("log")
    plt.ylabel('Signal Area / Tissue Area', fontsize=36, fontname='Liberation Sans')
    plt.xlabel('Projection depth [um]', fontsize=36, fontname='Liberation Sans')
    plt.xticks(fontsize=32, fontname='Liberation Sans')
    plt.yticks(fontsize=32, fontname='Liberation Sans')
    plt.legend(fontsize="24")
    # plt.ylim([0.005,0.025])
    # plt.ylim([0.005,0.035])
    plt.ylim([0,0.25])
    plt.xlim([-5,255])
    plt.tick_params(length=10, right=True)
    plt.tight_layout()
    plt.savefig(os.path.join('.',
                              'mean_signal-to-tissue_area_vs_confocal_depth_linear.pdf'), 
                  dpi=300, bbox_inches='tight')
    plt.show()

    

    
    
    # mean_signals_time = [] # as a control 
    # mean_signals_time_detrend = []
    
    # ref_image = []
    # psnr_metrics = []
    # sobel_metrics = [] 
    
    # poolsize = 2
    
    
    # total_area =[ ]
    # total_tissue_area = []
    # thresholds = []
    
    
    # print(len(imfiles))
    
    # for imfile_ii in np.arange(len(imfiles))[:]:
        
    #     basename = os.path.split(imfiles[imfile_ii])[-1]
    #     im = skio.imread(imfiles[imfile_ii])
        
    #     # im_mean = np.hstack([np.mean(imm) for imm in im])
    #     # mean_signals_time.append(im_mean)
    #     im_norm = im.copy()
        
    #     # # std norm per time slice. 
    #     # im_norm = np.array([(sli-sli.mean())/ sli.std() for sli in im]) # bleach correction .
    #     # # im_norm = np.array([sli/sli.mean() for sli in im]) # bleach correction .
    #     # im_norm = np.array([normalize(ndimage.gaussian_filter(ndimage.zoom(imm, zoom=[poolsize,poolsize], order=1, mode='reflect'),sigma=1), pmin=2, pmax=99.8, clip=True, eps=1e-20, dtype=np.float32) for imm in im_norm])
    #     # im_norm = np.array([normalize(ndimage.zoom(imm, zoom=[poolsize,poolsize], order=1, mode='reflect'),sigma=1), pmin=2, pmax=99.8, clip=True, eps=1e-20, dtype=np.float32) for imm in im_norm])
    #     im_norm = normalize(im_norm, clip=True)
    #     # # # grab the normalized11
    #     # # im_norm = np.uint8(255*im_norm)
        
    #     """
    #     1. background subtact in every timepoint. 
    #     """
    #     bg_norm = np.array([smooth_vol(sli, ds=16, smooth=5) for sli in im_norm]) # was 8 
    #     # bg_norm = np.array([smooth_vol(sli, ds=8, smooth=5) for sli in im_norm]) # was 8 
        
        
    #     tissue = np.max(im_norm, axis=0)
    #     tissue_binary = tissue> skfilters.threshold_multiotsu(tissue,3)[0]
    #     tissue_binary = ndimage.binary_fill_holes(tissue_binary)
        
    #     total_tissue_area.append(np.sum(tissue_binary))
        
        
        
    #     # check background normed. 
    #     im_norm_norm = im_norm - bg_norm
    #     # im_norm_norm = normalize(im_norm_norm)
    #     im_norm_norm = normalize(im_norm_norm, pmin=0, pmax=100, clip=True, eps=1e-20, dtype=np.float32)
        
    #     skio.imsave(os.path.join(savefolder,
    #                              'bg_subtract_'+basename), np.uint8(255*im_norm_norm))
        
    #     plt.figure(figsize=(10,10))
    #     plt.plot(im_norm[0][im_norm.shape[0]//2]); 
    #     plt.plot(bg_norm[0][bg_norm.shape[0]//2]); 
    #     plt.plot(im_norm_norm[0][im_norm_norm.shape[0]//2])
    #     plt.show()
        
        
    #     plt.figure(figsize=(10,10))
    #     plt.plot(im_norm[0][:,im_norm.shape[1]//2]); 
    #     plt.plot(bg_norm[0][:,bg_norm.shape[1]//2]); 
    #     plt.plot(im_norm_norm[0][:,im_norm_norm.shape[1]//2])
    #     plt.show()
        
        
    #     """
    #     2. detect oscillatory signal in every pixel domain!. 
    #     """
    #     pixel_timeseries = im_norm_norm.transpose(1,2,0).reshape(-1, len(im_norm_norm))
        
        
    #     pixel_timeseries_autocorr_raw = xcorr_timeseries_set_1d(pixel_timeseries[:], 
    #                                                         timeseries_array2=None, 
    #                                                         norm=True, 
    #                                                         eps=1e-12, 
    #                                                         stack_final=True)
    #     pixel_timeseries_autocorr_raw = pixel_timeseries_autocorr_raw / pixel_timeseries_autocorr_raw[:,0][:,None]
        
        
    #     plt.figure(figsize=(10,10))
    #     plt.plot(np.median(pixel_timeseries_autocorr_raw, axis=0))
    #     plt.show()
        
    #     autocorr_peaks = find_peaks(np.median(pixel_timeseries_autocorr_raw, axis=0), distance=10)
        
    #     print(autocorr_peaks)
    #     print('------')
    #     plt.figure(figsize=(10,10))
    #     plt.plot(np.median(pixel_timeseries_autocorr_raw, axis=0))
    #     plt.plot(autocorr_peaks[0], np.median(pixel_timeseries_autocorr_raw, axis=0)[autocorr_peaks[0]], 'go')
    #     plt.show()
        
        
    #     # ma average. # this is taking a while.... ## might parallelize this 
    #     # winsize = len(im_norm) // 20
    #     # winsize = 25
    #     winsize = 11 # this works.  # why ? we should get the average. 
    #     # winsize = 9 ### this is justified!  ish 
        
    #     #### seems we get 16 ....
        
        
    #     # winsize = autocorr_peaks[0][0]//4 + 1
    #     mode='reflect'
        
    #     pixel_timeseries_pad = np.pad(pixel_timeseries, pad_width=[[0,0], [winsize//2,winsize//2]], mode=mode)
    #     pixel_timeseries_ma = []
        
    #     for tt in tqdm(np.arange(pixel_timeseries.shape[1])):
            
    #         data = pixel_timeseries_pad[:,tt:tt+winsize].copy()
    #         data = np.nanmedian(data, axis=-1)
    #         # data = np.nanmean(data, axis=-1)
    #         pixel_timeseries_ma.append(data)
    #     pixel_timeseries_ma = np.array(pixel_timeseries_ma)
    #     pixel_timeseries_ma = pixel_timeseries_ma.T
 
        
    #     # for every pixel derive the autocorrelation! 
    #     pixel_timeseries_autocorr = xcorr_timeseries_set_1d(pixel_timeseries_ma[:], 
    #                                                         timeseries_array2=None, 
    #                                                         norm=True, 
    #                                                         eps=1e-12, 
    #                                                         stack_final=True)
        
    #     pixel_timeseries_autocorr = pixel_timeseries_autocorr / pixel_timeseries_autocorr[:,0][:,None]
        
        
    #     # find peaks
    #     all_inds = []
    #     all_periods = []
    #     all_peak_heights = []
    #     all_n_peaks = []
        
    #     for iii in tqdm(np.arange(len(pixel_timeseries_autocorr))):
    #         # inds, stats = find_peaks(pixel_timeseries_autocorr[283*186], height=0, distance=15)
    #         inds, stats = find_peaks(pixel_timeseries_autocorr[iii], height=0, distance=15)
            
    #         if len(inds)>0:
    #             all_inds.append(inds[0])
    #             all_periods.append(np.nanmean(np.diff(inds)))
    #             all_n_peaks.append(len(inds))
    #             all_peak_heights.append(np.nanmean(stats['peak_heights']))
    #         else:
    #             # why is this 0 ?
    #             all_inds.append(np.nan)
    #             all_periods.append(np.nan)
    #             all_n_peaks.append(0)
    #             all_peak_heights.append(np.nan)
                
                
    #     all_inds = np.hstack(all_inds)
    #     all_periods = np.hstack(all_periods)
    #     all_peak_heights = np.hstack(all_peak_heights)
    #     all_n_peaks = np.hstack(all_n_peaks)
        
        
    #     all_inds = all_inds.reshape(im_norm_norm[0].shape)
    #     all_periods = all_periods.reshape(im_norm_norm[0].shape)
    #     all_peak_heights = all_peak_heights.reshape(im_norm_norm[0].shape)
    #     all_n_peaks = all_n_peaks.reshape(im_norm_norm[0].shape)
        
        
        
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(all_inds)
    #     plt.show()
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(all_periods)
    #     plt.show()
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(all_peak_heights)
    #     plt.show()
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(all_n_peaks==1)
    #     plt.show()
            
            
    #     # threshold = skfilters.threshold_otsu(all_inds[np.logical_not(np.isnan(all_inds))])
    #     threshold = skfilters.threshold_multiotsu(all_inds[np.logical_not(np.isnan(all_inds))])[0] # to capture low.. stuff
    #     thresholds.append(threshold)
        
        
    #     # if imfile_ii == 0: 
    #     threshold_final = threshold # set the single threshold ? 1
        
        
        
    #     # binary_detection = np.logical_and(all_inds>threshold_final, all_inds < len(im)//2)
    #     binary_detection = all_inds>threshold_final
    #     binary_detection = binary_detection * tissue_binary
    #     # binary_detection = skmorph.binary_erosion(binary_detection, skmorph.disk(1))
    #     binary_detection = skmorph.remove_small_objects(binary_detection,min_size=25)
    #     binary_detection = skmorph.binary_closing(binary_detection, skmorph.disk(1)) # does make a diff
    #     # binary_detection = ndimage.binary_fill_holes(binary_detection)
    #     # binary_detection = skmorph.binary_dilation(binary_detection, skmorph.disk(1))
        
        
    #     labelled_detection = skmeasure.label(binary_detection)
        
    #     # remove large objects.
    #     regprops = skmeasure.regionprops(labelled_detection)
    #     areas = np.hstack([re.area for re in regprops])
        
    #     print(np.max(areas))
        
    #     remove_areas = np.setdiff1d(np.unique(labelled_detection), 0)[np.hstack(areas) > 10000]
        
    #     for rr in remove_areas:
    #         labelled_detection[labelled_detection==rr] = 0 
        
    #     binary_detection = labelled_detection>0
        
        
    #     total_area.append(np.sum(binary_detection))
        
        
    #     # I guess we replace with 0's in the first peak detections. and threshold out all lower (too fast)
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(binary_detection)
    #     plt.show()
        
        
    #     marked = sksegmentation.mark_boundaries(im_norm[0],
    #                                             binary_detection*1,
    #                                             color=(0,1,0))
        
    #     marked2 = sksegmentation.mark_boundaries(im_norm[0],
    #                                             labelled_detection,
    #                                             color=(0,1,0))
        
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(marked)
    #     plt.show()
        
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(marked2)
    #     plt.show()
        
        
    #     # save out. 
    #     savefolder_detections = os.path.join(savefolder, 
    #                                          basename)
    #     mkdir(savefolder_detections)
        
        
    #     all_marked2 = []
    #     for tt in np.arange(len(im_norm)):
            
    #         marked2 = sksegmentation.mark_boundaries(im_norm[tt],
    #                                                 labelled_detection,
    #                                                 color=(0,1,0))
    #         all_marked2.append(marked2)
    #     all_marked2 = np.array(all_marked2)
            
    #     skio.imsave(os.path.join(savefolder_detections,
    #                              'bg_subtract_'+basename), np.uint8(255*all_marked2))
            
        
        
        
        
    #     # take only the 1 half by definition of autocorr. (the first peak should be  zero. )
        
    #     # we detect the number of peaks
    #     # we detect the average period. (by the spacing between peaks!)
        
        
    #     # use joint as filter!. 
        
    # plt.figure(figsize=(10,10))
    # # plt.scatter(confocal_depth, sobel_metrics.ravel(), s=10, c='k')
    # plt.plot(confocal_depth, np.hstack(total_area)/np.hstack(total_tissue_area), 'k-o', lw=3)
    # plt.ylabel('Signal Area / Tissue Area', fontsize=20, fontname='Liberation Sans')
    # plt.xlabel('confocal depth [um]', fontsize=20, fontname='Liberation Sans')
    # plt.xticks(fontsize=18, fontname='Liberation Sans')
    # plt.yticks(fontsize=18, fontname='Liberation Sans')
    # # plt.ylim([0.005,0.025])
    # plt.xlim([-5,255])
    # plt.tick_params(length=10, right=True)
    # plt.savefig(os.path.join(savefolder,
    #                          'mean_Ca_signal_vs_confocal_depth.pdf'), 
    #              dpi=300, bbox_inches='tight')
    # plt.show()
        
        
        