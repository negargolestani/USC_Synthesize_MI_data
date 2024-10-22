import os
import numpy as np
import re
from pathlib import Path
import csv
import matplotlib.pyplot as plt

separator = '\t'
MHAD_dict = [
    #1: 
    'JumpingInPlace',
    # 2: 
    'JumpingJacks',
    # 3: 
    'Bending',
    # 4: 
    'Punching',
    # 5: 
    'WavingTwoHands',
    # 6: 
    'WavingOneHand',
    # 7: 
    'Clapping',
    # 8: 
    'Throwing',
    # 9: 
    'SitdownAndStandup',
    # 10: 
    'Sitdown',
    # 11: 
    'Standup'
    ]


#######################################################################################################
def read_txt(txtFile_path):
    txt = dict({ 'Time':list(), 'Data':list()  })

    with open(txtFile_path, 'r') as txtFile:

        # Read header
        header_list = txtFile.readline().strip().split('\t') 

        # Read data
        while True:
            line = txtFile.readline()
            if not line: break

            fields = line.strip().replace('i','j').split(separator)
            
            txt['Time'].append( float(fields[0]) )
            
            dt = list()
            for field in fields[1:]: dt.append( complex(field) )
            txt['Data'].append( dt )

    return txt
#######################################################################################################
def get_info(databaseName, fileName):
    # This Function returns dictionary with key values of activity types and values of their fileNames 
    #---------------------------------------------------------------
    if databaseName == 'BML': 
        parts = fileName.split('_')
        subject = parts[0]
        activity = parts[1]
    #---------------------------------------------------------------
    elif databaseName == 'MHAD': 
        parts = fileName.split('_')
        subject = parts[1][1:]
        activity = MHAD_dict[ int(parts[2][1:])-1 ]
    #---------------------------------------------------------------
    info = dict( {'Activity':activity, 'Subject':subject })
    return info
#######################################################################################################    
def load_file(databaseName, txtFile_path):
    txt = read_txt(txtFile_path)
    fileName = os.path.basename(txtFile_path)
    info = get_info(databaseName, fileName)
    return {**txt, **info}
#######################################################################################################
def database2numpy(databaseName, txtFolder_path):
    file_list = list()
    for fileName in os.listdir(txtFolder_path):
        if (fileName[-3:]!='txt'): continue
        txtFile_path = txtFolder_path  + '/' + fileName
        file_list.append( load_file(databaseName, txtFile_path) )
    return np.array(file_list)
#######################################################################################################



#######################################################################################################
if __name__ == "__main__":  # Generate Numpy Data
    
    database = 'MHAD'
    
    databaseFolder_path = './results/synth_MImotion/' + database
    numpyFile_path = databaseFolder_path+'/'+database+'.npy'

    database_np = database2numpy(database, databaseFolder_path)
    np.save( numpyFile_path, database_np )     # Save numpy file
#######################################################################################################


