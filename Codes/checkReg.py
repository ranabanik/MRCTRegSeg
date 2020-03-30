__author__ = 'rana.banik@vanderbilt.edu'

import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from skimage.exposure import rescale_intensity

parser = argparse.ArgumentParser()
parser.add_argument('Image1', type=str, help='provide file path to the image1, e.g.: c:/home/data/train/120.nii.gz')
parser.add_argument('Image2', type=str, help='provide file path to the image2, e.g.: c:/home/data/train/121.nii.gz')
# parser.add_argument('Destination', type=str, help='provide file save destination, e.g.: c:/home/data/train/')
args = parser.parse_args()

# dest = r'{}'.format(args.Destination)
img1 = nib.load(args.Image2).get_data()
# img1 = (img1 - img1.std())/img1.mean()
# print(img1.shape)
# img1 = nib.Nifti1Image(img1, affine=np.eye(4))
print(img1.shape, img1.dtype, img1.min(), img1.max())
img2 = nib.load(args.Image1).get_data()  # (256, 256, 170)
img2 = rescale_intensity(img2, in_range=(img2.min(), img2.max()), out_range=(img1.min(), img1.max()))
# img2 = (img2 - img2.std())/img2.mean()
print(img2.shape, img2.dtype)
assert img1.shape == img2.shape
# print(dest)
sag, cor, axi = img1.shape    # [cor, axi, sag]

axial = np.zeros(shape=(sag, cor))  # (256, 256)
coronal = np.empty(shape=(sag, axi), dtype=img1.dtype)   # (256, 170)
sagittal = np.empty(shape=(cor, axi), dtype=img1.dtype)    # (256, 170)

a_step = np.int(axi/5)  # todo: 5 chosen arbitrarily
s_step  = c_step = np.int(cor/8)  # todo: same dimensions
# print(a_step, s_step, c_step)  #  34 32 32
ind = 0
for c in range(0, cor, c_step): # a 0  c_step
    ind += 1
    for s in range(0, axi, s_step):
        ind += 1
        if (ind % 2) == 0:
            # axial[c:c+c_step, s:s+s_step] = img1[c:c+c_step, s:s+s_step,  np.int(dim[2]/2)].astype('int16')
            sagittal[c:c+c_step, s:s+s_step] = img1[np.int(sag/2)+1, c:c+c_step, s:s+s_step].astype('int16')
        else:
            sagittal[c:c+c_step, s:s+s_step] = img2[np.int(sag/2)+1, c:c+c_step, s:s+s_step].astype('int16')
        # break
    # break
plt.imshow(sagittal, cmap='gray')
plt.axis('off')
# plt.imshow(img2[:, :, 100], alpha=0.3)
plt.show()


