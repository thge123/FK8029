from os import listdir,remove

folders = ['X','V','O','R']
for folder in folders:
    for i in listdir(folder):
        remove('{}/{}'.format(folder,i))
