# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 12:08:49 2022

@author: Sveekruth

Description: Script to compute the compressibility of receptive field (RF) masked images.

Each image is 1280x720px. For color images, each pixel has 3 channels, each containing a uint8 (0-255) numeric value which carries a byte of info. This amounts to 2764800 bytes, approx. 2.8MB. No image should be larger than this upper bound. Our source files were in TIF format (uncompressed). These TIFs have a channel depth of 24 bits per pixel (bpp) (8bits/channel for 3 channels).

JPEG involves the discrete cosine transform (DCT) to discard high  frequency components in the image and results in lossy compression. This means there will be some reconstruction error while trying to retrieve the original image.
"""

# Import libraries:
import os, io, time
import numpy as np, scipy.io, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set(context='talk', style='dark')
import PIL

# Directory paths:
pdir = "E:" # "E:/IISc/BAI/PhD/V1ImagesGammaProject"
data_dir = f"{pdir}/data"
img_dir = f"{data_dir}/images/New_sets"
rf_dir = f"{data_dir}/rfData"
save_dir = f"{pdir}/IISc/BAI/PhD/V1ImagesGammaProject/programs/savedData/compressibility"

# Prerequisite vars:
subject = "kesariH"; # alpaH, kesariH, dona
elecs = np.arange(96) + 1 if subject == "dona" else np.arange(90) + 1
rms_elecs = scipy.io.loadmat(f"{rf_dir}/{subject}/{subject}MicroelectrodeRFData.mat")["highRMSElectrodes"].squeeze()

# Helper functions:
reorient_matlab = lambda M: M.transpose(1, 2, 0) # Reorients a 3D array as rows x columns x depth (MATLAB convention, also used by PIL and Matplotlib)

reorient_numpy = lambda M: M.transpose(2, 0, 1) # Reorients a 3D array as depth x rows x columns (NumPy convention, also used by PyTorch)

def monitor_XY_deg(view_dist=50, method=3, lab='Ray'):
    """
    This function returns the monitor x and y coordinate meshgrid in degrees, given monitor details
    and view distance in cm.

    Input arguments:
    view_dist - Viewing distance from screen center (cm)
    method - Method used to map pixels to degrees visual angle (dva)
    lab - Lab name, to fetch monitor specifications (in px and cm)

    Output arguments:
    X - X axis meshgrid matrix (x_res x y_res)
    Y - Same as above, for Y axis

    Typical use cases:
    i. [X, Y] = monitor_XY_deg; if analyzing data from Ray Lab, assuming none of the defaults need
    altering
    ii. [X, Y] = monitor_XY_deg(view_dist=23, lab='Arun'); if analyzing data from Arun Lab

    I have not included monitor_specs as an optional argument as creating it is usually more tedious
    than simply hardcoding it.
    """

    # DEFAULTS:
    if lab == 'Ray':
        monitor_specs = {'x_res':1280, 'y_res':720, 'width':53.086, 'height':29.972}
    elif lab == 'Arun':
        monitor_specs = {'x_res':1366, 'y_res':768, 'width':34.4232, 'height':19.3536}

    if method == 4: # Trigonometric method (most accurate, mathematically precise, no approximations
        x_axis_deg = np.arctan(np.linspace(-monitor_specs['width']/2, monitor_specs['width']/2,
                                    monitor_specs['x_res'])/view_dist)*(180/np.pi)
        y_axis_deg = np.arctan(np.linspace(-monitor_specs['height']/2, monitor_specs['height']/2,
                                    monitor_specs['y_res'])/view_dist)*(180/np.pi)

    x_deg = np.arctan(monitor_specs['width']/(2*view_dist))*(180/np.pi) # Half-width (degrees visual angle or dva)
    y_deg = np.arctan(monitor_specs['height']/(2*view_dist))*(180/np.pi) # Half-height (dva)

    if method == 1: # Full FOV method (includes both extremes at the expense of inexact pixel coordinates, but doesn't require dpp measurement):
        x_axis_deg = np.linspace(-x_deg, x_deg, monitor_specs['x_res'])
        y_axis_deg = np.linspace(-y_deg, y_deg, monitor_specs['y_res'])

    x_dpp = (2*x_deg)/monitor_specs['x_res'] # Degrees per pixel
    y_dpp = (2*y_deg)/monitor_specs['y_res']

    if method == 2: # Origin method (Supratim's default, includes the origin (0, 0) within machine precision but at the expense of full FOV at the positive ends of the axes):
        x_axis_deg = np.arange(-x_deg, x_deg, x_dpp)
        y_axis_deg = np.arange(-y_deg, y_deg, y_dpp)
    elif method == 3: # Pixel center method (more accurate, ensures every pixel number is mapped to its center):
        x_axis_deg = np.linspace(-(x_deg - x_dpp/2), (x_deg - x_dpp/2), monitor_specs['x_res'])
        y_axis_deg = np.linspace(-(y_deg - y_dpp/2), (y_deg - y_dpp/2), monitor_specs['y_res'])

    [X, Y] = np.meshgrid(x_axis_deg, y_axis_deg);

    return X, Y


    circ = lambda X, Y, h, k, r: ((X - h)**2 + (Y - k)**2) <= r**2 # Lambda function for circular boundaries


def patch(f, X, Y, rf_center=[0, 0], half_length=2, deg=True):
    """
    Returns an RF centered (in degrees) square patch from an image file. X and Y coordinate meshgrids
    to be provided. Patch dimensions in degrees or pixels as specified (follows half-space convention,
    i.e., the full square will be of side length 2*half_length).
    """
    x_axis_deg = x_axis_deg = X[0, :]
    y_axis_deg = Y[:, 0]

    x_dpp = x_axis_deg[1] - x_axis_deg[0]
    y_dpp = y_axis_deg[1] - y_axis_deg[0]

    azi_deg, ele_deg = rf_center

    # Convert to pixel coordinates
    azi_px = abs(x_axis_deg - azi_deg).argmin()
    ele_px = abs(y_axis_deg[::-1] - ele_deg).argmin()

    if deg: # Extract patch of 2 x specified degrees width and height
        x_ppd = 1/x_dpp
        y_ppd = 1/y_dpp
        try:
            len(half_length)
            hlx = int(round(half_length[0]*x_ppd, 0))
            hly = int(round(half_length[1]*y_ppd, 0))
        except:
            hlx = int(round(half_length*x_ppd, 0))
            hly = int(round(half_length*y_ppd, 0))

    else: # Extract patch of 2 x specified pixels width and height
        try:
            len(half_length)
            hlx = half_length[0]
            hly = half_length[1]
        except:
            hlx = half_length
            hly = hlx

    if f.mode == 'RGB': # Has color channels
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
    var1, cov, var2 = np.cov(np.concatenate([img1.ravel()[np.newaxis, :], img2.ravel()[np.newaxis, :]],
                                            axis=0))[np.triu_indices(2)]

    c1, c2 = (k1*L)**2, (k2*L)**2  # Two variables to stabilize the division with weak denominator

    return ((2*mu1*mu2 + c1)*(2*cov + c2))/((mu1**2 + mu2**2 + c1)*(var1 + var2 + c2))


def bpp(f, thr_q=50):
    """
    Computes the bits per pixel (bpp) of an image after applying JPEG compression (with a fixed quality
    hyperparameter 'q'). Returns the bpp, SSIM, and compressed image as a Numpy array.
    """

    with io.BytesIO() as buf:
        f.save(buf, format='JPEG', quality=thr_q)

        bits = buf.getbuffer().nbytes*8
        N = f.size[0]*f.size[1] # Pixel count
        bipp = bits/N

        g = PIL.Image.open(buf)
        if f.mode == 'RGB':
            f_array = reorient_numpy(np.asarray(f))
            f_compressed = reorient_numpy(np.asarray(g))
        else:
            f_array = np.asarray(f)  # Assumes grayscale or 'L'
            f_compressed = np.asarray(g)
        ssim = SSIM(f_array, f_compressed)

        return bipp, ssim, g


# Recreating RF table as DataFrame (from .mat):
mean_azi = scipy.io.loadmat(f"{rf_dir}/{subject}/{subject}MicroelectrodeRFData.mat")["rfStats"]["meanAzi"].squeeze()
mean_azi = np.array([np.ndarray.item(mean_azi[j]) if mean_azi[j].size else np.NaN for j in
                     range(mean_azi.size)])
mean_azi = np.concatenate([mean_azi, np.NaN*np.empty([len(elecs) - len(mean_azi)])])
mean_ele = scipy.io.loadmat(f"{rf_dir}/{subject}/{subject}MicroelectrodeRFData.mat")["rfStats"]["meanEle"].squeeze()
mean_ele = np.array([np.ndarray.item(mean_ele[j]) if mean_ele[j].size else np.NaN for j in
                     range(mean_ele.size)])
mean_ele = np.concatenate([mean_ele, np.NaN*np.empty([len(elecs) - len(mean_ele)])])

RF_data = pd.DataFrame(np.NaN*np.empty([elecs.size, 3]), index=[f"Elec{j}" for j in elecs],
                       columns=["mean_azi", "mean_ele", "high_RMS"])
RF_data.loc[:, "mean_azi"] = mean_azi
RF_data.loc[:, "mean_ele"] = mean_ele
RF_data.loc[:, "high_RMS"] = [j in rms_elecs for j in elecs]
del elecs, mean_azi, mean_ele, rms_elecs

# Logical mask:
X, Y = monitor_XY_deg(view_dist=50, method=3, lab='Ray') # Follows matrix Y axis convention
# Y = Y[::-1, :] # UD flip to follow image Y axis convention

# Compressibility:
categories = ["A", "F", "H", "T", "L"]
img_types = ["color", "grayscale"]
for category in categories:
    num_images = int(len([f for f in os.listdir(f"{img_dir}/{category}") if f[:5] == "Image"])/2)
    for img_type in img_types:
        img_init_idx = 1 + num_images*(img_type == "grayscale")
        C = pd.DataFrame(np.NaN*np.empty([num_images, RF_data.shape[0]]),
                         index=[f"Image{i}" for i in (np.arange(num_images) + img_init_idx)],
                         columns=RF_data.index)
        for i in C.index:
            f = PIL.Image.open(f"{img_dir}/{category}/{i}.tif").convert('RGB')
            for j in C.columns:
                if RF_data.loc[j, "high_RMS"]:
                    g = patch(f, X, Y, rf_center=RF_data.loc[j, ["mean_azi", "mean_ele"]])
                    C.loc[i, j] = -bpp(g)[0]
        if os.path.isfile(f"{save_dir}/{subject}.xlsx"):
            with pd.ExcelWriter(f"{save_dir}/{subject}.xlsx", engine="openpyxl", mode='a') as writer:
                C.to_excel(excel_writer=writer, sheet_name=f"{category}_{img_type}",
                           index_label='Index')
        else:
            with pd.ExcelWriter(f"{save_dir}/{subject}.xlsx", engine="xlsxwriter") as writer:
                C.to_excel(excel_writer=writer, sheet_name=f"{category}_{img_type}",
                           index_label='Index')