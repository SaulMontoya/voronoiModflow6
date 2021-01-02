import os, flopy

name = 'Model'
workspace = '../model'
mf_exe_name = '../exe/mf6'

fname = os.path.join(workspace, name + '.hds')
hdobj = flopy.utils.HeadFile(fname, precision='double')
head = hdobj.get_data()
fname = os.path.join(workspace, name + '.cbc')
bdobj = flopy.utils.CellBudgetFile(fname, precision='double', verbose=False)

print(head)
