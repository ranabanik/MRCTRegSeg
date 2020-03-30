import nibabel as nib
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

foldDir = r'/media/banikr2/DATA/emModels/MatlabTestRun/Temp'
foldPath = glob(os.path.join(foldDir, '*'))
# print(foldPath) # gets all the folders.

ind = 0
CTint = []#np.array([])
MRint = []#np.array([])
subjects = []
for fp in foldPath:
    ind += 1
    # print(fp[-3:])   # 262
    imgPath = glob(os.path.join(fp, '*.nii.gz'))
    if len(imgPath) != 3: # 0 - CT, 1 - MR, 2 - MR_norm
        print("File missing!!, {}".format(fp[:-3]))
        break
    subjects.append(fp[-3:])
    CTobj = nib.load(imgPath[0])
    CT = CTobj.get_data()
    CTint.append([CT.min(), CT.max()]) #np.append(CTint, (CT.min(), CT.max()), axis=1)
    MRobj = nib.load(imgPath[1])
    MR = MRobj.get_data()
    MRnormobj = nib.load(imgPath[2])
    MRnorm = MRnormobj.get_data()
    MRint.append([MR.min(), MR.max(), MRnorm.min(), MRnorm.max()])
print(MRint, np.shape(MRint))
fig, ax = plt.subplots(figsize=(10, 10))
ax = plt.figure().gca()
# ax.plot(np.array(CTint)[:, 0], linewidth=1.0, label="CT min")
# ax.plot(np.array(CTint)[:, 1], linewidth=1.0, label="CT max")
ax.plot(np.array(subjects), np.array(MRint)[:, 1]/10, linewidth=2.0, label="MR max")
ax.plot(np.array(subjects), np.array(MRint)[:, 0], linewidth=2.0, label="MR min")
ax.plot(np.array(subjects), np.array(MRint)[:, 3], linewidth=2.0, label="MR norm max")
ax.plot(np.array(subjects), np.array(MRint)[:, 2], linewidth=2.0, label="MR norm min")
ax.grid()
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xticks(fontsize=20, rotation=90)
plt.yticks(fontsize=20)
leg = ax.legend(fontsize=20, loc='best')
plt.show()
