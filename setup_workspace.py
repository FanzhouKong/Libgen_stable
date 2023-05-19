import os
import sys
workspace_parent_dir = sys.argv[1]
library_name = sys.argv[2]

workspace_dir = os.path.join(workspace_parent_dir, library_name)
existance_status = os.path.exists(workspace_dir)
if existance_status==False:
    # print(workspace_dir)
    os.mkdir(workspace_dir)
    os.chdir(workspace_dir)
    subdirs = ['standard_list','spectra', 'features', 'curated_library', 'figures']
    for sub in subdirs:
        if os.path.exists(sub)==False:
            os.makedirs(sub)
        else:
            print('the workspace has been created')
            continue
else:
    # os.chdir(workspace_dir)
    print('the workspace is built')


print(workspace_dir)