import os

for i in range(20,60,10):
    cmd = ('./output ' + str(i) + ' > output' + str(i) + '.csv')
    os.system(cmd)

