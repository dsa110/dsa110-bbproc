#s script launches beamnsplice on all triggers
# verifies that all voltage files exist, and that trigger has not yet been processed

import os
import os.path
import sys
import numpy as np
import glob

dirname = sys.argv[1];

## lists all specnums that have already been processed
listfiles = glob.glob('/mnt/data/dsa110/findpulse/*.fil');
allspecnum = [];
for fname in listfiles:
    fi = fname.find('_');
    se = fname[fi+1:].find('_');
    if fname[:fi] == dirname:
        allspecnum.append(int(fname[fi+1:fi+1+se]));

## lists triggers from T2 (not all of them have voltage dumps associated with them)
flname = '/mnt/data/dsa110/T2/' + dirname + '/fl_output.dat';
with open(flname, "r") as f:
    trilist = f.readlines();

idxtoremove = [];
for kk in range(1,len(trilist)):
    # reads list of fl_ names for each trigger
    flnames = trilist[kk].split(' ')[1:];
    flnames[-1] = flnames[-1][:-1];
    specnum  = int(flnames[0][flnames[0].find('.out')+5:]);
    
    # check that specnum has not already been processed
    if specnum not in allspecnum:
        # if not all files listed, remove existing voltage files (verification in case of overwritten files)
        if len(flnames) < 16:
            os.system('rm -r /mnt/data/dsa110/T3/corr*/' + dirname + '/*' + str(specnum) + '*');
            print('removed specnum '+ str(specnum));
            idxtoremove.append(kk);
        else:
            # verify that all files exist
            fullnames = [];
            for k in range(len(flnames)):
                if k == 0:
                    fullnames.append('/home/user/data/T3/corr00/' + dirname + '/'+flnames[k]);
                else:
                    fullnames.append('/home/user/data/T3/corr' + str(k+1).zfill(2) + '/' + dirname + '/'+flnames[k]);
            # delete triggers if not all files exist
            fileexist = 0;
            for fn in fullnames:
                if not os.path.isfile(fn):
                    fileexist += 1;
            if fileexist != 0:
                os.system('rm -r /mnt/data/dsa110/T3/corr*/' + dirname + '/*' + str(specnum) + '*');
                print('removed specnum '+ str(specnum));
                idxtoremove.append(kk);
    else:
        idxtoremove.append(kk);

trilist = np.delete(trilist,idxtoremove);

for k in range(1,len(trilist)):
    try:
        os.system('python beamnsplice.py ' + dirname + ' /mnt/data/dsa110/T3/calibs/beamformer_weights_corrXX_25feb21.dat ' + trilist[k].split(' ')[0] + ' &');
        print('processing tigger #'+trilist[k].split(' ')[0]);
    except:
        print('could not process trigger ' + trilist[k].split(' ')[0]+'\n');
