import os

for i in range(4,20):
    cmd = ('./output ' + str(i) + ' > output' + str(i) + '.csv')
    os.system(cmd)

