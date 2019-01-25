# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 12:48:33 2018

@author: Shiva
"""

import numpy as np

def Zmatrixb(conf):
    Zabc = np.matrix([[0.0+0.0j, 0.0-0.0j, 0.0-0.0j], 
                       [0.0+0.0j, 0.0-0.0j, 0.0-0.0j], 
                       [0.0+0.0j, 0.0-0.0j, 0.0-0.0j]]);
    if conf==1:          
        Zabc = np.matrix([[0.4576 + 1.0780j,   0.1560 + 0.5017j,   0.1535 + 0.3849j], 
                          [0.1560 + 0.5017j,   0.4666 + 1.0482j,   0.1580 + 0.4236j], 
                          [0.1535 + 0.3849j,  0.1580 + 0.4236j,   0.4615 + 1.0651j]]);
    if conf==2:          
        Zabc = np.matrix([[0.4666 + 1.0482j,   0.1580 + 0.4236j,   0.1560 + 0.5017j],
                          [0.1580 + 0.4236j,   0.4615 + 1.0651j,   0.1535 + 0.3849j],
                          [0.1560 + 0.5017j,   0.1535 + 0.3849j,   0.4576 + 1.0780j]]);
    if conf==3:          
        Zabc = np.matrix([[0.4615 + 1.0651j,   0.1535 + 0.3849j,   0.1580 + 0.4236j],
                          [0.1535 + 0.3849j,   0.4576 + 1.0780j,   0.1560 + 0.5017j],
                          [0.1580 + 0.4236j,   0.1560 + 0.5017j,   0.4666 + 1.0482j]]);
    if conf==4:          
        Zabc = np.matrix([[0.4615 + 1.0651j,   0.1580 + 0.4236j,   0.1535 + 0.3849j],
                          [0.1580 + 0.4236j,   0.4666 + 1.0482j,   0.1560 + 0.5017j],
                          [0.1535 + 0.3849j,   0.1560 + 0.5017j,   0.4576 + 1.0780j]]);
    if conf==5:          
        Zabc = np.matrix([[0.4666 + 1.0482j,   0.1535 + 0.3849j,   0.1560 + 0.5017j],
                          [0.1535 + 0.3849j,   0.4615 + 1.0651j,   0.1580 + 0.4236j],
                          [0.1560 + 0.5017j,   0.1580 + 0.4236j,   0.4576 + 1.0780j]]);
    if conf==6:          
        Zabc = np.matrix([[0.4576 + 1.0780j,   0.1560 + 0.5017j,   0.1535 + 0.3849j],
                          [0.1560 + 0.5017j,   0.4666 + 1.0482j,   0.1580 + 0.4236j],
                          [0.1535 + 0.3849j,   0.1580 + 0.4236j,   0.4615 + 1.0651j]]);
    if conf==7:          
        Zabc = np.matrix([[0.4576 + 1.0780j,   0.0000 + 0.0000j,   0.1535 + 0.3849j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.1535 + 0.3849j,   0.0000 + 0.0000j,   0.4615 + 1.0651j]]);
    if conf==8:          
        Zabc = np.matrix([[0.4576 + 1.0780j,   0.1535 + 0.3849j,   0.0000 + 0.0000j],
                          [0.1535 + 0.3849j,   0.4615 + 1.0651j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j]]);
    if conf==9:          
        Zabc = np.matrix([[1.3292 + 1.3475j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j]]);
    if conf==10:          
        Zabc = np.matrix([[0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   1.3292 + 1.3475j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j]]);
    if conf==11:          
        Zabc = np.matrix([[0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   0.0000 + 0.0000j],
                          [0.0000 + 0.0000j,   0.0000 + 0.0000j,   1.3292 + 1.3475j]]); 
    
    if conf==11:          
        Zabc = np.matrix([[1.5209 + 0.7521j,   0.5198 + 0.2775j,   0.4924 + 0.2157j],
                          [0.5198 + 0.2775j,   1.5329 + 0.7162j,   0.5198 + 0.2775j],
                          [0.4924 + 0.2157j,   0.5198 + 0.2775j,   1.5209 + 0.7521j]]);

    r_bb=np.real(Zabc[1,1]);   
    x_bb=np.imag(Zabc[1,1]);
    r_ba=np.real(Zabc[1,0]);  
    x_ba=np.imag(Zabc[1,0]);
    r_bc=np.real(Zabc[1,2]);  
    x_bc=np.imag(Zabc[1,2]);      
    return (r_bb,x_bb,r_ba,x_ba,r_bc,x_bc)