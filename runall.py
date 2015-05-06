import os

for i in range(3,20):
    cmd = ('./output ' + str(i) + ' > output' + str(i) + '.csv')
    os.system(cmd)

