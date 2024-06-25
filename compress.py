# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 12:08:49 2022

@author: Sveekruth

Description: Script to find compressibility of all images in the Natural Images Dataset.
"""

# Import libraries:
import os, io, time
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set(context='talk')
import PIL

# Directory paths:
pdir = "E:/IISc/BAI/PhD/V1ImagesGammaProject"
data_dir = f"{pdir}/data"
img_dir = f"{data_dir}/images"
resource_dir = "E:/IISc/BAI/PhD/NIA_22/Resources"

# Helper functions:
reorient_matlab = lambda M: M.transpose(1, 2, 0) # Reorients a 3D array as rows x columns x depth (MATLAB convention, also used by PIL and Matplotlib)
reorient_numpy = lambda M: M.transpose(2, 0, 1) # Reorients a 3D array as depth x rows x columns (NumPy convention, also used by PyTorch)
    
def patch(f, rf_center, half_length=2, deg=True, \
          monitor_specs={'height':11.8, 'width':20.9, 'x_res':1280, 'y_res':720}, \
          view_dist=50/2.54):
    """
    Returns an RF centered (in degrees) square patch from an image file. Monitor 
    dimensions in inches and pixels. Patch dimensions in degrees or pixels as 
    specified (follows half-space convention, i.e., the full square will be of 
    side length 2*patch_deg or 2*patch_px).
    """  
    
    x_deg = np.arctan((monitor_specs['width']/2)/view_dist)*180/np.pi # Half the space, in degrees
    y_deg = np.arctan((monitor_specs['height']/2)/view_dist)*180/np.pi 

    x_axis_deg = np.linspace(-x_deg, x_deg, monitor_specs['x_res'])
    y_axis_deg = np.linspace(-y_deg, y_deg, monitor_specs['y_res'])

    azi_deg, ele_deg = rf_center
    
    # Convert to pixel coordinates
    azi_px = abs(x_axis_deg - azi_deg).argmin()
    ele_px = abs(y_axis_deg[::-1] - ele_deg).argmin()
    
    if deg: # Extract patch of 2 x specified degrees width and height
        x_ppd = (monitor_specs['x_res']/2)/x_deg
        y_ppd = (monitor_specs['y_res']/2)/y_deg
        if type(half_length) == np.ndarray:
            hlx = int(round(half_length[0]*x_ppd, 0))
            hly = int(round(half_length[1]*y_ppd, 0))
        else:
            hlx = int(round(half_length*x_ppd, 0))
            hly = int(round(half_length*y_ppd, 0))
        
    else: # Extract patch of 2 x specified pixels width and height
        if type(half_length) == np.ndarray:
            hlx = half_length[0]
            hly = half_length[1]
        else:
            hlx = half_length
            hly = hlx
            
    if f.mode == 'RGB': # Has color channels. Assumes NumPy 3D convention.
        f_array = reorient_numpy(np.asarray(f)) 
        g_array = f_array[:, (ele_px - hly):(ele_px + hly), (azi_px - hlx):(azi_px + hlx)]
        g = PIL.Image.fromarray(reorient_matlab(g_array)) # Patch image file
    else: # Assumes 2D grayscale or 'L'
        f_array = np.asarray(f)
        g_array = f_array[(ele_px - hly):(ele_px + hly), (azi_px - hlx):(azi_px + hlx)]
        g = PIL.Image.fromarray(g_array)
        
    return g

   
def SSIM(img1, img2, L=2**8-1, k1=1e-2, k2=3e-2):
    """
    Takes two image arrays im1 (typically original) and im2 (typically 
    compressed) and computes the Structural Similarity Index Metric [-1, 1] 
    between them. L is the Dynamic Range, a constant = 2^bpc - 1, where bpc is
    8 by default. k1 and k2 are small constants.
    Refer: https://en.wikipedia.org/wiki/Structural_similarity
    """
    
    mu1, mu2 = img1.mean(), img2.mean()
    var1, cov, var2 = np.cov(np.concatenate([img1.ravel()[np.newaxis, :], \
                           img2.ravel()[np.newaxis, :]], axis=0))\
                            [np.triu_indices(2)]
                            
    c1, c2 = (k1*L)**2, (k2*L)**2  # Two variables to stabilize the division with weak denominator
    
    return ((2*mu1*mu2 + c1)*(2*cov + c2))/((mu1**2 + mu2**2 + c1)*(var1 + var2 + c2))


def bpp(f, thr_q=50):
    """
    An alternative approach to bpp computation. Keeps JPEG quality hyper-
    parameter 'q' fixed. Returns the bpp, SSIM, and compressed image.
    """
       
    if f.mode == 'RGB':
        f_array = reorient_numpy(np.asarray(f))
    else: # Assumes grayscale or 'L'
        f_array = np.asarray(f)
    
    buf = io.BytesIO()
    f.save(buf, format='JPEG', quality=thr_q)
        
    bits = buf.getbuffer().nbytes*8
    N = f_array.shape[-2]*f_array.shape[-1] # Pixel count
    bpp = bits/N 
    
    if f.mode == 'RGB':
        buf = reorient_numpy(np.asarray(PIL.Image.open(buf)))
    else:
        buf = np.asarray(PIL.Image.open(buf))
    ssim = SSIM(f_array, buf)
    
    buf = PIL.Image.fromarray(reorient_matlab(buf), mode='RGB')
        
    return bpp, ssim, buf


def alt_bpp(f, thr_ssim=0.99, q=50, q_upper=100, q_lower=1):
    """
    Takes a PIL Image object and computes its compressed bits per pixel (bpp).
    Uses JPEG compression and a binary search algorithm to ensure that the
    compressed file has an SSIM > 0.99 (default).
    """
    
    if f.mode == 'RGB':
        f_array = reorient_numpy(np.asarray(f))
    else: # Assumes grayscale or 'L'
        f_array = np.asarray(f)
        
    # Working with a RAM buffer to increase speed and avoid read/write to flash memory (REFER: https://stackoverflow.com/questions/30771652/how-to-perform-jpeg-compression-in-python-without-writing-reading):
    # REFER: https://jdhao.github.io/2019/07/06/python_opencv_pil_image_to_bytes/
    buf = io.BytesIO() # Memory (bytes) buffer. Read more on bytearray() too.
    f.save(buf, format='JPEG', quality=q) # The buffer essentially allows us to store data in a variable instead of a file.
    while (q_upper - q_lower) != 1:
        if f.mode == 'RGB':
            buf = reorient_numpy(np.asarray(PIL.Image.open(buf)))
        else:
            buf = np.asarray(PIL.Image.open(buf))
        if SSIM(f_array, buf) > thr_ssim:
            q_upper = q
            q = int(round((q + q_lower)/2, 0))
        else:
            q_lower = q
            q = int(round((q + q_upper)/2, 0))
        buf = io.BytesIO()
        f.save(buf, format='JPEG', quality=q)
    buf = io.BytesIO() # Memory (bytes) buffer. Read more on bytearray() too.
    f.save(buf, format='JPEG', quality=q_upper)
    
    bits = buf.getbuffer().nbytes*8
    N = f_array.shape[-2]*f_array.shape[-1] # Pixel count
    bpp = bits/N
    
    buf = PIL.Image.open(buf)
    # if f.mode == 'RGB':
    #     buf = reorient_numpy(np.asarray(PIL.Image.open(buf)))
    # else:
    #     buf = np.asarray(PIL.Image.open(buf))
                
    return bpp, q_upper, buf

"""
Each image is 1280x720px. Each pixel has 3 channels, each containing a uint8 
(0-255) numeric value which carries a byte of info. This amounts to 2764800 
bytes, approx. 2.8MB. No image should be larger than this upper bound. Our 
source files were in TIF format, then converted to PNG, both which perform 
lossless compression.

NOTE: These TIFs have a channel depth of 24 bits per pixel (bpp) (8bits/channel 
for 3 channels) while PNGs are 32 bpp (8bpc with an extra alpha transparency
channel). However, on inspection it turns out the transparency values for all
images used here are unity (or 255 throughout), making it redundant.

However, JPEG involves the discrete cosine transform (DCT) to discard high 
frequency components in the image and results in lossy compression. This
means there will be some reconstruction error while trying to retrieve the
original image.
"""
# Loading RF Data for patch centering:
RF_df = pd.read_csv(f"{resource_dir}/rfStats.csv")
RF_df.index = [f"Elec{j + 1}" for j in RF_df.index]
RF_df.loc['Elec90', :] = np.NaN
RF_df.dropna(how='all', axis=0, inplace=True) # Removes all bad electrodes
RF_df.drop('Elec28', axis=0, inplace=True) # 28 is not a high RMS electrode

# Preparing DF to store compressibility (-bpp):
index = pd.Index([f"Image{i}" for i in range(1, 33)]) # Max 32 images in each directory
columns = pd.Index([f"Elec{i}" for i in range(1, 91)]) # 90 electrodes in total in the (81 grid + 9 ECoG)
df_C = pd.DataFrame(np.NaN*np.empty([len(index), len(columns)]), index=index, \
                    columns=columns)  
# stack_C = np.empty(np.concatenate([df_C.shape, [3, 86, 92]], axis=0), \
#                                   dtype='uint8')
# stack_C = np.concatenate([stack_C, stack_C], axis=0)

# Iterating -bpp over all RF centered patches and images as a measure of compressibility:
if not os.path.isfile(f"{resource_dir}/compress.xlsx"):
    start = time.time()
    sets = ['TL', 'AF'] 
    for c, s in enumerate(sets):
        for i in index:
            if os.path.isfile(f"{img_dir}/Images{s}/{i}.png"):
                f = PIL.Image.open(f"{img_dir}/Images{s}/{i}.png").convert('RGB')
                for j in RF_df.index:
                    df_C.loc[i, j] = -bpp(patch(f, RF_df.loc[j, :'meanEle']))[0]
                    # df_C = -df_C # Convention
                    # stack_C[32*c + int(i[5:]) - 1, int(j[4:]) - 1]
        if c == 0:
            with pd.ExcelWriter(f"{resource_dir}/compress.xlsx", engine="xlsxwriter") as writer:  
                df_C.to_excel(writer, f"{s}", index_label='Index')
        else:
            with pd.ExcelWriter(f"{resource_dir}/compress.xlsx", engine="openpyxl", mode='a') as writer:  
                df_C.to_excel(writer, f"{s}", sheet_name=None,\
                              index_label='Index')
        df_C.loc[:, :] = np.NaN
    stop = time.time()
    print(f"Time Elapsed: {stop - start}s")
else:
    df_C = pd.read_excel(f"{resource_dir}/compress.xlsx", sheet_name=None, index_col='Index')
        
# Looking at what Sudhanshu has done for Predictability:
if not os.path.isfile(f"{resource_dir}/predict.xlsx"):
    with np.load("predictability_96_RF.npz") as f: 
        # print(f.files)
        P = f['predictability']
    sets = ['TL', 'AF'] 
    for c, s in enumerate(sets):
        if c == 0:
            with pd.ExcelWriter(f"{resource_dir}/predict.xlsx", engine="xlsxwriter") as writer:  
                pd.DataFrame(P[32:64, :], index=index, columns=RF_df.index).\
                             to_excel(writer, f"{s}", index_label='Index')
        else:
            with pd.ExcelWriter(f"{resource_dir}/predict.xlsx", engine="openpyxl", mode='a') as writer:  
                pd.DataFrame(P[:32, :], index=index, columns=RF_df.index).\
                             to_excel(writer, f"{s}", index_label='Index')
else:
    df_P = pd.read_excel(f"{resource_dir}/predict.xlsx", sheet_name=None, index_col='Index')
    
with np.load(f"{resource_dir}/predictions_96_RF.npz") as f:
    stack_P = f['stack_P']

fig1, ax1 = plt.subplots(1, 2)
ax1[0].bar(RF_df.index, df_P['TL'].corrwith(df_C['TL'].dropna(axis=1, how='all')))
ax1[1].bar(RF_df.index, df_P['AF'].corrwith(df_C['AF'].dropna(axis=1, how='all')))
# df_P = pd.read_excel(f"{resource_dir}/predict.xlsx", sheet_name='TL', \
#                      index_col='Index') 
# df_C.dropna(axis=1, how='all', inplace=True)                    
# foo = df_P.corrwith(df_C, axis=0)

def displayPC(elec, dataset, subset):
    """
    Displays a 3x17 figure with image patches corresponding to valid electrode
    RF, their P and C values and reconstructions. Correlation plots are also
    shown.
    """
    assert elec in [int(i[4:]) for i in RF_df.index] # List of valid electrodes
    dataset_dict = {"ImagesTL":np.s_[:32, :], "ImagesAF":np.s_[32:64, :],\
                    "ImagesH":np.s_[64:, :]}
    assert dataset in dataset_dict.keys()
    if dataset == "ImageH":
        assert subset == 1 # ImagesH has no 2nd subset
    stack_P_dataset = stack_P[dataset_dict[dataset]]
    # with sns.set_style("dark") as style:
    fig, ax = plt.subplots(3, 16)
    for c, i in enumerate(range(1 + 16*(subset - 1), 17 + 16*(subset - 1))):
        f = PIL.Image.open(f"{img_dir}/{dataset}/Image{i}.png").convert('RGB')
        g = patch(f, RF_df.loc[f"Elec{elec}", :"meanEle"])
        ax[0, c].imshow(g)
        ax[1, c].imshow(stack_P_dataset[i - 1, \
                        np.where(RF_df.index == f"Elec{elec}")[0].item()])    
        ax[2, c].imshow(bpp(g)[-1])
    fig.subplots_adjust(wspace=0, hspace=0)
    return fig

"""
reject = [5, 12, 28, 29] # Electrodes to be rejected
LFP_electrodes = list(set(range(1, 82)) - set(reject)) # 77 in total
spike_electrodes = [int(s.split('_')[0][4:]) for s in os.listdir(f"{data_dir}/"
                    f"alpaH/Microelectrode/240817/GRF_002/segmentedData/Spikes")\
                    if 'SID1' in s]
spike_electrodes.sort()
"""

"""
if not os.path.isfile(f"{resource_dir}/TL_compress.csv"):
    # os.chdir(f"{img_dir}/ImagesTL") # Colored Texture and Landscape images
    
    # Images x Electrodes DataFrame with -bpp values (to be used)
    index = pd.Index([f"Image{i}" for i in range(1, 33)]) # Max 32 images in each directory
    columns = pd.Index([f"Elec{i}" for i in range(1, 82)]) # 81 electrodes in total in the grid
    df_C = pd.DataFrame(np.NaN*np.empty([32, 81]), index=index, columns=columns)  
    
    # Iterating -bpp over all VALID RF centered patches and images as a measure of compressibility
    for i in index:
        # if os.path.isfile(f"{i}.png"):
        if os.path.isfile(f"{img_dir}/ImagesTL/{i}.png"):
            # f = PIL.Image.open(f"{i}.png").convert('RGB')
            f = PIL.Image.open(f"{img_dir}/ImagesTL/{i}.png").convert('RGB')
            for j in LFP_electrodes:
                df_C.loc[i, f"Elec{j}"] = \
                    -bpp(patch(f, RF_df.loc[j, :'meanEle'], \
                     half_length=RF_df.loc[j, ['rfSizeAzi', 'rfSizeEle']].values, deg=True))
    # df.to_csv(f"{resource_dir}/TL_compress.csv")
    # df_compress.to_csv(f"{resource_dir}/TL_compress.csv")
    with pd.ExcelWriter(f"{resource_dir}/TL.xlsx", engine="openpyxl", mode='a') as writer:  
        df_C.to_excel(writer, "C", index_label='Index')
    # os.chdir(f"{resource_dir}")
else:
    # os.chdir(resource_dir)
    # df_compress = pd.read_csv(f"{resource_dir}/TL_compress.csv", index_col='Unnamed: 0')
    df_C = pd.read_excel(f"{resource_dir}/TL.xlsx", sheet_name='C', index_col='Index')    

# df_gamma = pd.read_csv(f"{resource_dir}/TL_gamma.csv", index_col='Row')
df_actual = pd.read_excel(f"{resource_dir}/TL.xlsx", sheet_name='Actual', index_col='Index')
# corrs = np.corrcoef(df_compress.dropna(axis=1).T, df_gamma.dropna(axis=1).T).\
#     diagonal(len(LFP_electrodes))

# sns.set_palette("Spectral", n_colors=32)

fig1, ax1 = plt.subplots(1, 2, sharey=True)
fig1.suptitle('Compressibility and Predictability')
for i in df_C.index:
    ax1[0].plot(df_C.loc[i, :], df_actual.loc[i, :], 'o', alpha=0.7, label=i)
ax1[0].legend()
ax1[0].set_title('Gamma Response vs. Compressibility')
ax1[0].set_ylabel('Gamma Response (30-80Hz, ST/BL)')
ax1[0].set_xlabel('-bpp (bits per pixel)')
ax1[0].set_ylim([0, None])
ax1[0].set_xlim([None, round(df_C.max().max())])

with np.load(f"{resource_dir}/pred_77.npz") as f: # Sudhanshu has currently sent 77 and 31 respectively
    # print(f.files)
    P = f['predictibility']
df_P =  pd.DataFrame(P, df_C.index, [f'Elec{e}' for e in \
                                                    LFP_electrodes]) 
for i in df_pred_1x1.index:
    ax1[1].plot(df_pred_1x1.loc[i, :], df_gamma.dropna(axis=1).loc[i, :], 'o', alpha=0.7, label=i)
ax1[1].set_title('Gamma Response vs. Predictability')
ax1[1].set_ylabel('Gamma (30-80Hz) deltaFFT (max - mean)')
ax1[1].set_xlabel('Predictability')
ax1[1].set_xlim([0, 1.1]) 
    
corrs = np.corrcoef(df_compress_1x1.dropna(axis=1).T, df_gamma.dropna(axis=1).T).\
    diagonal(len(LFP_electrodes))
    
fig2 = plt.figure()
ax2 = plt.axes()
ax2 = sns.violinplot(data=corrs)
ax1.set_title('Gamma Response vs. Compressibility')
ax1.set_ylabel('Gamma (30-80Hz) deltaFFT (max - mean)')
ax1.set_xlabel('-bpp (bits per pixel)')
ax1.set_ylim([0, None])
ax1.set_xlim([None, 0])
"""

"""
# Relationship between compressibility and predictability (courtesy Sudhanshu):
with np.load('pred_77.npz') as f: # Sudhanshu has currently sent 77 and 31 respectively
    # print(f.files)
    pred = f['predictibility']
df_pred77 =  pd.DataFrame(pred, df_compress.index, [f'Elec{e}' for e in \
                                                    LFP_electrodes])
"""

"""
# Comprehensive Dataframe with file sizes:
col_l1 = [d for d in os.listdir() if ((d[:6] == 'Images') & \
                                     (os.path.isdir(d)))] # Ensures image directory is chosen
col_l2 = ['TIF', 'PNG'] # Within every image directory, there are subdirectories to contain any non PNG images.
columns = pd.MultiIndex.from_product([col_l1, col_l2], names=['lvl1', 'lvl2'])    
index = pd.Index([f"Image{i}" for i in range(1, 33)], name='S. No.') # Max 32 images in each directory
df = pd.DataFrame(np.NaN*np.ones([32, len(col_l1)*len(col_l2)]), index=index, columns=columns)  

for d in col_l1:
    for t in col_l2:
        if t != 'PNG':
            if os.path.isdir(f"{img_dir}/{d}/{t}"):
                os.chdi'r(f"{img_dir}/{d}/{t}")
            else:
                continue # https://towardsdatascience.com/how-to-use-break-continue-and-pass-in-python-programming-9cd841763b3c
        else:
            os.chdir(f"{img_dir}/{d}")
        for i in index:
            if os.path.isfile(f"{i}.{t.lower()}"):
                df.loc[i, pd.IndexSlice[d, t]] = os.stat(f"{i}.{t.lower()}").st_size # https://stackoverflow.com/questions/25189575/pandas-dataframe-select-columns-in-multiindex
        os.chdir(img_dir)
"""

"""
The code so far works as expected. For ImagesH, there is only one category,
so after Image16, we have NaNs. For ImagesAS and ImagesTS, TIF folders are
absent.

NOTE: For ImageH TIFs, the sizes are ~5MB. Turns out they are 1600x1067px.
"""  

"""
# Emperical assessment of JPEG Compression Quality:
do_assessment = False
if do_assessment:    
    
# Creating Images x Electrodes Dataframes for all categories, storing patch compressibility:

    os.chdir('ImagesAF')
    f = PIL.Image.open('Image1.png')
    foo = reorient_numpy(np.asarray(f))[:3, :, :] # Example image
    bar = patch(foo, [2, -3], half_length=80, deg=False) # Patch from example image

    fig1, ax1 = plt.subplots(1, 2)
    ax1[0].imshow(reorient_matlab(foo)) # extent=[-x_deg, x_deg, -y_deg, y_deg])
    ax1[1].imshow(reorient_matlab(bar))

    f2 = PIL.Image.open('Image22.png')
    foo2 = reorient_numpy(np.asarray(f2))[:3, :, :] # Example image
    bar2 = patch(foo2, [2, -3], half_length=80, deg=False)  # Patch from example image

    fig2, ax2 = plt.subplots(1, 2)
    ax2[0].imshow(reorient_matlab(foo2)) # extent=[-x_deg, x_deg, -y_deg, y_deg])
    ax2[1].imshow(reorient_matlab(bar2))
    
    g = PIL.Image.fromarray(reorient_matlab(bar), mode='RGB')
    g.save('test1.jpg', 'JPEG', quality=100, subsampling=0)
    g2 = PIL.Image.fromarray(reorient_matlab(bar2), mode='RGB')
    g2.save('test2.jpg', 'JPEG', quality=100, subsampling=0)

    q_index = pd.Index([q for q in range(100, 0, -1)], name='q')
    baz_df = pd.DataFrame(np.NaN*np.empty([100, 4]), index=q_index, columns=\
                          ['SSIM_1', 'Compression_1', 'SSIM_2', 'Compression_2'])
    for q in q_index:
        if q != 100:
            g.save('test1_compressed.jpg', 'JPEG', quality=q)
            g2.save('test2_compressed.jpg', 'JPEG', quality=q)
        else:
            g.save('test1_compressed.jpg', 'JPEG', quality=q, subsampling=0)
            g2.save('test2_compressed.jpg', 'JPEG', quality=q, subsampling=0)
        baz = reorient_numpy(np.asarray(PIL.Image.open('test1_compressed.jpg')))
        baz2 = reorient_numpy(np.asarray(PIL.Image.open('test2_compressed.jpg'))) 
        baz_df.loc[q, :] = {'q':q, 'SSIM_1':SSIM(bar, baz), 'Compression_1':\
                            os.stat('test1_compressed.jpg').st_size, 'SSIM_2':\
                                SSIM(bar2, baz2), 'Compression_2':os.stat\
                                     ('test2_compressed.jpg').st_size}
    
    baz_df.loc[:, ['Compression_1', 'Compression_2']] /= \
        baz_df.loc[:, ['Compression_1', 'Compression_2']].max() # w.r.t. an uncompressed JPEG baseline
    
    fig3 = baz_df.plot()
    fig3.set_xlim(101, 0)
    
    test = baz_df.loc[:, 'SSIM_1'] > 0.99 # REPLACE with a user defined SSIM threshold
    q1 = np.arange(100, 1, -1)[test.values[:-1] ^ test.values[1:]].item() # To identify the lowest q value that is super-threshold
    test = baz_df.loc[:, 'SSIM_2'] > 0.99 # REPLACE with a user defined SSIM threshold
    q2 = np.arange(100, 1, -1)[test.values[:-1] ^ test.values[1:]].item() # To identify the lowest q value that is super-threshold
    del test

    # Summary of the above procedure:
    g.save('test1_compressed.jpg', 'JPEG', quality=q1)
    baz = reorient_numpy(np.asarray(PIL.Image.open('test1_compressed.jpg')))
    g2.save('test2_compressed.jpg', 'JPEG', quality=q2)
    baz2 = reorient_numpy(np.asarray(PIL.Image.open('test2_compressed.jpg')))
    fig4, ax4 = plt.subplots(2, 2, sharex=True, sharey=True)
    ax4[0, 0].imshow(reorient_matlab(bar))
    ax4[0, 1].imshow(reorient_matlab(baz))
    ax4[1, 0].imshow(reorient_matlab(bar2))
    ax4[1, 1].imshow(reorient_matlab(baz2))
"""