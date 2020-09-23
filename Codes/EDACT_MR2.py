import nibabel as nib
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pickle
import h5py

# foldDir = r'/media/banikr2/DATA/emModels/MatlabTestRun/Temp'
# foldDir = r'/media/banikr2/DATA/emModels/CT-MR_batch2'
foldDir = r'/media/banikr2/DATA/banikr_D_drive/DataBackup/MRCT_result'
foldPath = glob(os.path.join(foldDir, '*'))
# print(foldPath[:) # gets all the folders.
figPath = r'/home/banikr2/MATLAB/MatlabProjects/MRCTRegSeg/Figures'
if __name__ != '__main__':
    # print(foldPath[155:157])
    ind = 0
    CTint = []  #np.array([])
    MRint = []  #np.array([])
    subjects = []
    for fp in foldPath:
        print(ind+1, fp[-4:])   # todo: change to -3
        # if __name__ != '__main__':
        imgPath = glob(os.path.join(fp, '*.nii.gz'))
        # print(imgPath[1])
        # if len(imgPath) != 3:   # 0 - CT, 1 - MR, 2 - MR_norm #todo: comment out the loop
        #     print("File missing!!, {}".format(fp[-3:]))
        #     continue
        ind += 1
        # if ind == 5:
        #     break
        subjects.append(fp[-4:])
        # CTobj = nib.load(imgPath[0])
        # CT = CTobj.get_data()
        # CTint.append([CT.min(), CT.max()])  # np.append(CTint, (CT.min(), CT.max()), axis=1)
        # MRobj = nib.load(imgPath[1])
        # MR = MRobj.get_data()
        MRnormobj = nib.load(imgPath[1]) #todo: change to 2 for batch2
        MRnorm = MRnormobj.get_data()
        # MRint.append([MR.min(), MR.max(), MRnorm.min(), MRnorm.max()])
        MRint.append([MRnorm.min(), MRnorm.max()])

    # print(MRint, np.shape(MRint))
    with h5py.File(os.path.join(figPath, 'MRnormIntRangeBatch1.h5'), 'w') as pfile: #todo: change the name
        pfile['subID'] = np.array(subjects, dtype=h5py.string_dtype(encoding='utf-8'))
        # pfile['CT'] = CTint
        pfile['MR_norm1'] = MRint
# ax = plt.figure().gca()


    f = h5py.File(os.path.join(figPath, 'MRnormIntRangeBatch1.h5'), 'r') #todo: grab the new name
    print(f.keys())
    print(f['subID'])
    # CTint = f['CT']
    MRint = f['MR_norm1'] #todo: norm
    subjects = f['subID']
    print(np.array(subjects))

    fig, ax = plt.subplots(1, figsize=(10, 5), dpi=150) # 2, 1
    # axes = ax.ravel()
    # axes[0].plot(np.array(subjects), np.array(CTint)[:, 0], linewidth=2.0, label="CT norm min", linestyle='--')
    # axes[0].plot(np.array(subjects), np.array(CTint)[:, 1], linewidth=2.0, label="CT norm max", linestyle='--')
    # axes[0].grid()
    # axes[0].legend(fontsize=15, loc='best')
    # axes[0].xaxis.set_ticks(np.array(subjects))
    # axes[0].xaxis.set_ticklabels(np.array(subjects), fontsize=15, rotation=90)
    # axes[0].yaxis.set_ticks()#, fontsize=15, rotation=45)
    # axes[1].plot(np.array(subjects), np.array(MRint)[:, 1]/10, linewidth=2.0, label="1/10 MR max")
    # axes[1].plot(np.array(subjects), np.array(MRint)[:, 0], linewidth=1.0, label="MR min")
    ax.plot(np.array(subjects), np.array(MRint)[:, 1], linewidth=2.0, label="MR norm max") #todo: 1 -> 3
    ax.plot(np.array(subjects), np.array(MRint)[:, 0], linewidth=2.0, label="MR norm min") #todo: 0 -> 2
    ax.grid()
    ax.set_title('Intensity range of {} MR norm Images'.format(len(subjects)), fontsize=20) # todo: don't change 0
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(fontsize=15, rotation=90)
    plt.yticks(fontsize=15, rotation=45)
    ax.legend(fontsize=15, loc='best')
    plt.xlabel("Subject IDs",fontsize=15)
    plt.savefig(os.path.join(figPath, 'MR-batch1_norm.png'))
    plt.show()

if __name__ == '__main__':
    abnormalMaxMR2 = [264, 287, 305, 319, 419, 432, 437, 509]
    filePath = os.path.join(figPath, 'MRCTIntRange.h5')
    f = h5py.File(filePath, 'r')
    # print(f.keys())
    MR = f['MR']
    subID = f['subID']
    ind = 0
    print(MR[:, 1])
    MRmaxsum = 0
    # subID = np.array(subID)
    print(type(subID[0]))
    for i in subID:
        if np.int(i) in abnormalMaxMR2:
            continue
        ind += 1
        MRmaxsum += MR[ind-1, 1]

    print("mean max MR intensity without abnormals: ", MRmaxsum/ind) # 690.8620689655172

    plt.plot(np.array(subID), MR[:, 1], linewidth=2.0, label="1/10 MR max")
    plt.show()