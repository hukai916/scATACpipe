import sys
import shutil

keepfile = sys.argv[1]
allfile  = sys.argv[2]

keep_list = [x.split()[0] for x in open(keepfile)]
allfile   = [x.strip() for x in open(allfile)]

for file in allfile:
    todel = 1
    for keep in keep_list:
        if keep in file:
            todel = 0
    if "singularity" in file:
        todel = 0
    if todel:
        print("Deleting: ", file)
        shutil.rmtree(file)
