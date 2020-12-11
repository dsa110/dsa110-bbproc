import sys
import os

datadir = sys.argv[1];
val = str((int(sys.argv[2])-477)*16);

f = open(os.path.join(datadir,"dumps.dat"),'r');

tmp = f.readline();
while tmp != '':
    tmp = f.readline();
    if tmp.find(val)!=-1:
        break;

f.close();

if tmp == '':
    print('file not found.')
else:
    filnum = int(tmp[tmp.find('DUMPNUM')+7:tmp.find('DUMPNUM')+7+tmp[tmp.find('DUMPNUM')+7:].find(' ')]);

    print('file is:');
    print(os.path.join(datadir,'fl_' + str(filnum) +'.out'));
