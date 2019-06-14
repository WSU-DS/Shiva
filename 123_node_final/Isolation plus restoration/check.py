# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 11:33:04 2018

@author: spoudel
"""

import restoration as rs
import Fault_Isolation as fi
import numpy as np

fault_location=99

op_sw_fi = fi.Fault_Isolation(fault_location)

Sec, Tie=rs.restoration(op_sw_fi)

print ('Sec \n')
print(Sec)
print ('\n Tie \n')
print(Tie)

MRIDS_LIST=['_F5C161E8-198A-51DD-E61E-FCD5D5711602','_A84EF66A-5EC4-46FF-63BC-225E101A5BCD','_55AA66AB-D8C4-50AF-FC21-41DD46D2FBFA',\
           '_18E47F11-2454-B57A-DBA0-E02E8A046FD7','_A3DF5CA3-7C17-0534-9B32-2B554E12FFA6','_4116F5F0-2986-C31E-29E9-89E1315B09AE',\
           '_6A8AC2EC-1D7E-FDCD-71D1-6ABE39BE2294','_C0D74ECE-E5D3-B01D-2B20-EC466686DE3E']
SWITCH_MRIDS= []
nSec=Sec.__len__()
for k in range(0,nSec):
    if Sec[k]==0:
        SWITCH_MRIDS.append(MRIDS_LIST[k])            
        
nTie=Tie.__len__();
for k in range(0,nTie):
    if Tie[k]==1:
        SWITCH_MRIDS.append(MRIDS_LIST[k+5])           
          
for str in SWITCH_MRIDS:
    print(str)


op_sw=Sec.__len__()-np.count_nonzero(Sec)
cls_sw=Tie.__len__()-np.count_nonzero(Tie)
print('The switch statuses are')
print(op_sw)
print(cls_sw)

