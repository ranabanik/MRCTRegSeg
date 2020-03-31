import nibabel as nib
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pickle
import h5py

# foldDir = r'/media/banikr2/DATA/emModels/MatlabTestRun/Temp'
foldDir = r'/media/banikr2/DATA/emModels/CT-MR_batch2'
foldPath = glob(os.path.join(foldDir, '*'))
# print(foldPath) # gets all the folders.
figPath = r'/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Figures'
if __name__ != '__main__':
    # print(foldPath[155:157])
    ind = 0
    CTint = []  #np.array([])
    MRint = []  #np.array([])
    subjects = []
    for fp in foldPath:
        print(ind+1, fp[-3:])   # todo: change to -3
        # if __name__ != '__main__':
        imgPath = glob(os.path.join(fp, '*.nii.gz'))
        # print(imgPath[0])
        if len(imgPath) != 3:   # 0 - CT, 1 - MR, 2 - MR_norm
            print("File missing!!, {}".format(fp[-3:]))
            continue
        ind += 1
        # if ind == 5:
        #     break
        subjects.append(fp[-3:])
        CTobj = nib.load(imgPath[0])
        CT = CTobj.get_data()
        CTint.append([CT.min(), CT.max()])  # np.append(CTint, (CT.min(), CT.max()), axis=1)
        MRobj = nib.load(imgPath[1])
        MR = MRobj.get_data()
        MRnormobj = nib.load(imgPath[2])
        MRnorm = MRnormobj.get_data()
        MRint.append([MR.min(), MR.max(), MRnorm.min(), MRnorm.max()])

    # print(MRint, np.shape(MRint))
    with h5py.File(os.path.join(figPath, 'MRCTIntRange.h5'), 'w') as pfile:
        pfile['subID'] = np.array(subjects, dtype=h5py.string_dtype(encoding='utf-8'))
        pfile['CT'] = CTint
        pfile['MR'] = MRint
# ax = plt.figure().gca()


f = h5py.File(os.path.join(figPath, 'MRCTIntRange.h5'), 'r')
print(f.keys())
print(f['subID'])
CTint = f['CT']
MRint = f['MR']
subjects = f['subID']
print(np.array(subjects))

fig, ax = plt.subplots(2, 1, figsize=(32, 10), dpi=150)
axes = ax.ravel()
axes[0].plot(np.array(subjects), np.array(CTint)[:, 0], linewidth=2.0, label="CT norm min", linestyle='--')
axes[0].plot(np.array(subjects), np.array(CTint)[:, 1], linewidth=2.0, label="CT norm max", linestyle='--')
axes[0].grid()
axes[0].legend(fontsize=15, loc='best')
axes[0].xaxis.set_ticks(np.array(subjects))
axes[0].xaxis.set_ticklabels(np.array(subjects), fontsize=15, rotation=90)
# axes[0].yaxis.set_ticks()#, fontsize=15, rotation=45)
axes[1].plot(np.array(subjects), np.array(MRint)[:, 1]/10, linewidth=2.0, label="1/10 MR max")
axes[1].plot(np.array(subjects), np.array(MRint)[:, 0], linewidth=1.0, label="MR min")
axes[1].plot(np.array(subjects), np.array(MRint)[:, 3], linewidth=2.0, label="MR norm max")
axes[1].plot(np.array(subjects), np.array(MRint)[:, 2], linewidth=1.0, label="MR norm min")
axes[1].grid()
axes[0].set_title('Intensity range of {} Images'.format(len(subjects)), fontsize=20)
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xticks(fontsize=15, rotation=90)
plt.yticks(fontsize=15, rotation=45)
axes[1].legend(fontsize=15, loc='best')
plt.xlabel("Subject IDs",fontsize=15)
plt.savefig(os.path.join(figPath, 'MRCT-batch2.png'))
plt.show()
