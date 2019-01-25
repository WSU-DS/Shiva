# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 09:35:36 2018
@author: Shiva
"""

import numpy as np
from pulp import *
import Zmatrixa as zma
import Zmatrixb as zmb
import Zmatrixc as zmc

def restoration(F):
    #data = f.get('data/variable1') 
    #data = np.array(data)
    # import Real_Power_FLow_Linear as rpfl   
    #######################################################################################
    nNodes=125;
    nEdges=129;
    DG=[95, 115, 122];
    Tie_Switches=np.matrix([[54, 94], [117, 123], [DG[0], 125],[DG[1], 125],[DG[2], 125]]);
    edges = np.loadtxt('edges.txt'); 
    demand = np.loadtxt('LoadData.txt'); 
    loops = np.loadtxt('cycles.txt'); 
    LineData=np.loadtxt('Line_Config.txt');
    mult=-1*(demand[:,1]+demand[:,3]+demand[:,5]);
    v_i=range(0,nNodes);
    s_i=range(0,nNodes);
    x_ij=range(0,nEdges);
    V_i=range(0,nNodes);
    
    # Different variables for optimization function
    #######################################################################################
    assign_vars_vi = LpVariable.dicts("xv", v_i, 0, 1,  LpBinary)
    assign_vars_si = LpVariable.dicts("xs", s_i, 0, 1,  LpBinary)
    assign_vars_xij = LpVariable.dicts("xl", x_ij, 0, 1,  LpBinary)
    assign_vars_xij0 = LpVariable.dicts("xl0", x_ij, 0, 1,  LpBinary)
    assign_vars_xij1 = LpVariable.dicts("xl1", x_ij, 0, 1,  LpBinary)
    assign_vars_Pija = LpVariable.dicts("xPa", x_ij, -10000, 10000)
    assign_vars_Pijb = LpVariable.dicts("xPb", x_ij, -10000, 10000)
    assign_vars_Pijc = LpVariable.dicts("xPc", x_ij, -10000, 10000)
    assign_vars_Qija = LpVariable.dicts("xQa", x_ij, -10000, 10000)
    assign_vars_Qijb = LpVariable.dicts("xQb", x_ij, -10000, 10000)
    assign_vars_Qijc = LpVariable.dicts("xQc", x_ij, -10000, 10000)
    assign_vars_Via = LpVariable.dicts("xVa", V_i, 0.95, 1.05)
    assign_vars_Vib = LpVariable.dicts("xVb", V_i, 0.95, 1.05)
    assign_vars_Vic = LpVariable.dicts("xVc", V_i, 0.95, 1.05)
    ##########################################################################################
    
    # Indices for tie switches and virtual switches to insert into objetive functions
    N=[124, 125, 126, 127, 128];
    # Optimization problem objective definitions
    prob = LpProblem("Resilient Restoration",LpMinimize)
    prob += lpSum(assign_vars_si[k] * mult[k] for k in s_i)+\
            lpSum(assign_vars_xij[N[k]] for k in range(0,5))
       
    # Constraints (v_i<=0)
    for k in range(0,nNodes):
        prob += assign_vars_vi[k] <= 1
        
    # Constraints (s_i<=v_i)
    for k in range(0,nNodes):
        prob += assign_vars_si[k] <= assign_vars_vi[k]
        
    # Constraints (x_ij<=v_i*v_j). This is non-linear and is linearized here..
    for k in range(0,nEdges):
        prob += assign_vars_xij[k] <= assign_vars_vi[edges[k,0]-1]
        prob += assign_vars_xij[k] <= assign_vars_vi[edges[k,1]-1]
        
    # Constraints (x_ij0+x_ij1<=x_ij)
    for k in range(0,nEdges):
        prob += assign_vars_xij0[k] + assign_vars_xij1[k] <= assign_vars_xij[k]
        
    ##########################################################################################
    #########################################################################################  
    # Now writing the real power flow equation for Phase A phase B and Phase C
    fr=edges[:,0];
    to=edges[:,1];
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0,pa.__len__());
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0,ch.__len__());  
        # The overall power flow equation for Phase A now can be written as,    
        prob += lpSum(assign_vars_Pija[pa[j]] for j in N) - demand[ed,1] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Pija[ch[j]] for j in M)
     
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0,pa.__len__());
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0,ch.__len__());   
        # The overall power flow equation for Phase B now can be written as,    
        prob += lpSum(assign_vars_Pijb[pa[j]] for j in N)-demand[ed,3] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Pijb[ch[j]] for j in M)
    #    
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0, pa.__len__());
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0, ch.__len__());   
        # The overall power flow equation for Phase C now can be written as,    
        prob += lpSum(assign_vars_Pijc[pa[j]] for j in N)-demand[ed,5] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Pijc[ch[j]] for j in M)
    #    
    ##########################################################################################
    #
    ## Now imposing the big-M method to ensure the real-power flowing in open line is zero
    ## -M * x_ij0 <= Pij <= x_ij1* M 
    M=100000;
    for k in range(0, nEdges):    
        prob += assign_vars_Pija[k] <= M * assign_vars_xij1[k]
        prob += assign_vars_Pijb[k] <= M * assign_vars_xij1[k] 
        prob += assign_vars_Pijc[k] <= M * assign_vars_xij1[k]     
        prob += assign_vars_Pija[k] >= -M * assign_vars_xij0[k]
        prob += assign_vars_Pijb[k] >= -M * assign_vars_xij0[k] 
        prob += assign_vars_Pijc[k] >= -M * assign_vars_xij0[k] 
    
    ##########################################################################################
    ##########################################################################################
    #  
    ## Now writing the reactive power flow equation for Phase A phase B and Phase C
    fr=edges[:,0];
    to=edges[:,1];
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0, pa.__len__());
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0, ch.__len__());   
        # The overall power flow equation for Phase A now can be written as,    
        prob += lpSum(assign_vars_Qija[pa[j]] for j in N)-demand[ed,2] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Qija[ch[j]] for j in M)
    #
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0, pa.__len__());
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0, ch.__len__());   
        # The overall power flow equation for Phase B now can be written as,    
        prob += lpSum(assign_vars_Qijb[pa[j]] for j in N)-demand[ed,4] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Qijb[ch[j]] for j in M)
    #   
    for k in range(0, nEdges):    
        ed=int(edges[k,1]-1);
        node=edges[k,1];
        # Finding the all parent nodes of a particular node
        pa=np.array(np.where(to==edges[k,1]));
        pa=pa.flatten();
        N=range(0, pa.__len__()); 
        # Finding the all children nodes of a particular node
        ch=np.array(np.where(fr==edges[k,1]));
        ch=ch.flatten();
        M=range(0, ch.__len__());  
        # The overall power flow equation for Phase C now can be written as,    
        prob += lpSum(assign_vars_Qijc[pa[j]] for j in N)-demand[ed,6] * assign_vars_si[node-1]== \
        lpSum(assign_vars_Qijc[ch[j]] for j in M)
    #    
    ###########################################################################################
    #
    ## Now imposing the big-M method to ensure the reactive-power flowing in open line is zero
    ## -M * x_ij0 <= Qij <= x_ij1* M
    M=100000;
    for k in range(0, nEdges):    
        prob += assign_vars_Qija[k] <= M * assign_vars_xij[k]
        prob += assign_vars_Qijb[k] <= M * assign_vars_xij[k] 
        prob += assign_vars_Qijc[k] <= M * assign_vars_xij[k] 
        prob += assign_vars_Qija[k] >= -M * assign_vars_xij0[k]
        prob += assign_vars_Qijb[k] >= -M * assign_vars_xij0[k] 
        prob += assign_vars_Qijc[k] >= -M * assign_vars_xij0[k]
    
    ##########################################################################################
    ##########################################################################################
    
    # Now the voltage constraints are written as set of inequality constraints by coupling them with
    # line or switch variable.
    base_Z=4.16**2;
    #prob = LpProblem("Resilient Restoration",LpMinimize)
    #prob += lpSum(assign_vars_si[k] * mult[k] for k in s_i)
    for k in range(0, 124):
        conf=LineData[k,3];
        len=LineData[k,2];
        # Get the Z matrix for a line
        r_aa,x_aa,r_ab,x_ab,r_ac,x_ac = zma.Zmatrixa(conf);
        line=[LineData[k,0], LineData[k,1]];
        prob += assign_vars_Via[int(line[0])-1]-assign_vars_Via[int(line[1])-1] - \
        2*r_aa*len/(5280*base_Z*1000)*assign_vars_Pija[int(line[1])-1]- \
        2*x_aa*len/(5280*base_Z*1000)*assign_vars_Qija[int(line[1])-1]+ \
        (r_ab+np.sqrt(3)*x_ab)*len/(5280*base_Z*1000)*assign_vars_Pijb[int(line[1])-1] +\
        (x_ab-np.sqrt(3)*r_ab)*len/(5280*base_Z*1000)*assign_vars_Qijb[int(line[1])-1] +\
        (r_ac-np.sqrt(3)*x_ac)*len/(5280*base_Z*1000)*assign_vars_Pijc[int(line[1])-1] +\
        (x_ac+np.sqrt(3)*r_ac)*len/(5280*base_Z*1000)*assign_vars_Qijc[int(line[1])-1] ==0
        
    for k in range(0, 124):
        conf=LineData[k,3];
        len=LineData[k,2];
        # Get the Z matrix for a line
        r_bb,x_bb,r_ba,x_ba,r_bc,x_bc = zmb.Zmatrixb(conf);
        line=[LineData[k,0], LineData[k,1]];
        prob += assign_vars_Vib[int(line[0])-1]-assign_vars_Vib[int(line[1])-1] - \
        2*r_bb*len/(5280*base_Z*1000)*assign_vars_Pijb[int(line[1])-1]- \
        2*x_bb*len/(5280*base_Z*1000)*assign_vars_Qijb[int(line[1])-1]+ \
        (r_ba-np.sqrt(3)*x_ba)*len/(5280*base_Z*1000)*assign_vars_Pija[int(line[1])-1] +\
        (x_ba+np.sqrt(3)*r_ba)*len/(5280*base_Z*1000)*assign_vars_Qija[int(line[1])-1] +\
        (r_bc+np.sqrt(3)*x_bc)*len/(5280*base_Z*1000)*assign_vars_Pijc[int(line[1])-1] +\
        (x_bc-np.sqrt(3)*r_bc)*len/(5280*base_Z*1000)*assign_vars_Qijc[int(line[1])-1] ==0
        
    for k in range(0, 124):
        conf=LineData[k,3];
        len=LineData[k,2];
        # Get the Z matrix for a line
        r_cc,x_cc,r_ca,x_ca,r_cb,x_cb = zmc.Zmatrixc(conf);
        line=[LineData[k,0], LineData[k,1]];
        prob += assign_vars_Vic[int(line[0])-1]-assign_vars_Vic[int(line[1])-1] - \
        2*r_cc*len/(5280*base_Z*1000)*assign_vars_Pijc[int(line[1])-1]- \
        2*x_cc*len/(5280*base_Z*1000)*assign_vars_Qijc[int(line[1])-1]+ \
        (r_ca+np.sqrt(3)*x_ca)*len/(5280*base_Z*1000)*assign_vars_Pija[int(line[1])-1] +\
        (x_ca-np.sqrt(3)*r_ca)*len/(5280*base_Z*1000)*assign_vars_Qija[int(line[1])-1] +\
        (r_cb-np.sqrt(3)*x_cb)*len/(5280*base_Z*1000)*assign_vars_Pijb[int(line[1])-1] +\
        (x_cb+np.sqrt(3)*r_cb)*len/(5280*base_Z*1000)*assign_vars_Qijb[int(line[1])-1] ==0
        
    prob += assign_vars_Via[124] == 1
    prob += assign_vars_Vib[124] == 1
    prob += assign_vars_Vic[124] == 1
    ##########################################################################################
    ##########################################################################################
    
    # Insert Cyclic constratins
    nC=loops.__len__();
    for k in range(0,nC):
        Sw=loops[k];
        nSw_C=np.count_nonzero(Sw);
        #print (nSw_C);
        prob += lpSum(assign_vars_xij[Sw[j]-1] for j in range(0,nSw_C)) <= nSw_C-1
            
    # Insert fault in the system
    nF=F.__len__();
    for k in range(0, nF):
        prob += assign_vars_xij[F[k]] == 0
        
    ##########################################################################################
    ##########################################################################################   
    prob.solve()
    
    varsdict = {}
    for v in prob.variables():
        varsdict[v.name] = v.varValue
    
    ##########################################################################################
    ##########################################################################################   

    
    ##########################################################################################
    ##########################################################################################    

    
    Sec=[varsdict['xl_0'], varsdict['xl_120'], varsdict['xl_115'], varsdict['xl_118'],\
       varsdict['xl_123'],varsdict['xl_119']]
    
    Tie=[varsdict['xl_124'], varsdict['xl_125'] ]  
   
    return (Sec, Tie)
        
    
    
    
    
    
    
    
    
