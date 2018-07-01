function cplexmilp
%% 
    try
    tic
    % This is a restoration code for multi-feeder test case incorporating
    % 3-phase power flow
    % In this code the number of switches in 4 Feeder 1069 node are 
    % 40*4 for four feeder
    % 4 switches for connecting feeder to  sub-transmission node
    % 7 tie-switches among the feeder
    % 4 DG modeled as virtual tie-switch
    
    clear all;
    clc;
    nNodes=1317;
    nEdges=1327;
    nSwitches=175;
    load 'Line_OH.txt'; ov_line=Line_OH;
    load 'Line_UG.txt'; ug_line=Line_UG;
    load 'Fuse.txt';    fuse=Fuse;
    load 'Xfrm.txt';    tran=Xfrm;
    load 'Swtch.txt';   swtch=Swtch;
    
    load 'all_edges.txt';
    load 'all_edges_graph.txt';
    edges=all_edges;
    %Defining the x_ij variable for only switch such that PF are written in
    %combination of x_ij for switches only and for other it is written
    %simple
    line=all_edges_graph;  
    load 'load_data.txt'
    %make loadata in KW;
    power=load_data;    
    load Z_mat_UG;
    load Z_matrix_OH;
     
    %% 
    %%   
    %Node assignment variables
    f = [-(power(:,2)+power(:,4)+power(:,6)).*ones(nNodes, 1)];      
    %Branch energized variables ;
    f=[f;zeros((nSwitches-11),1);ones(11,1)];      
    
    % Power flowing on each branch for P_ij;
    f=[f;zeros(nSwitches*3,1)];    
    % Product of varibales for x_ij*P_ij;
    f=[f;zeros(nEdges*3,1)];   
    
    % Power flowing on each branch for Q_ij;
    f=[f;zeros(nSwitches*3,1)];    
    % Product of varibales for x_ij*Q_ij;
    f=[f;zeros(nEdges*3,1)];   
    
    %Voltage equations for all nodes. The voltage equations are written as,
    % x_ij*V_j= [V_i-(r_ij*P_ij+x_ij*Q_ij)]*x_ij. So the number of
    % variables will be nEdges*2
    f=[f;zeros(nNodes*3,1)];
    f=[f;zeros(nSwitches*6,1)];
    %% 
    %%    
    fr=edges(:,1);
    t=edges(:,2);    
    G1 = graph(fr,t);   
    plot(G1,'r','LineWidth',2);       
    
    Source=1317;    
    A4=zeros(nNodes-1,size(f,1));       
    %inequality for each node to restore or not
    for i=1:nNodes
          A4(i,i)=1;
    end
    Aineq=A4;
    bineq=ones(nNodes,1);  
    %%
    %%     
    %Linearzing product of two binary variables for finding the lines
    %energied from nodes energized.     
    edges=G1.Edges;
    edges=table2array(edges);  
    % Adding tie-switches in the test case. Source is the sub-transmission
    % node of the test case. There are 4 DG in the system at nodes
    mul=329; 
    Tie_Switches=[261 263+mul; 248  254+mul*3; 262 258+mul; 244+mul*2 257+mul*3; 236+mul 256+mul*2; 75+mul 252+mul*2;...
        266+mul*2 252+mul*3; Source 264; Source 577; Source 733; Source 1253];
    edges1= flipud(dfsearch(G1,Source,'edgetonew'));
    
    %X_ij<=v_i*v_j is only for switches not all lines  
    swtch=[swtch; Tie_Switches];
    [Aineq, bineq]=edges_from_nodes(Aineq,bineq,f,nNodes,swtch,nSwitches);
        
    %% 
    %%  
    % This is for finding the index of swtch where they are in edges1. The
    % PF equation starts from the begining of edges1.    
    ind=1;
    store_=[];
    edges1=[edges1;Tie_Switches];
    for k=1:nEdges
        line1=edges1(k,:);
        %Find this line1 in swtch for index
        a=line1(1,1); b=line1(1,2);
        see=find(a==swtch(:,1));    
        for m=1:size(see,1)
            if b==swtch(see(m),2)
                vec1(ind)=nNodes+see(m);
                ind=ind+1;
                store_=[store_ k];
            end    
        end        
        see=find(a==swtch(:,2));
        
        for m=1:size(see,1)
            if b==swtch(see(m),1)
                vec1(ind)=nNodes+see(m);
                ind=ind+1;
                store_=[store_ k];
            end    
        end 
    end      
              
    %%   
    %%          
    %Linearzing product of x_ij and P_ij for three phases for formulating the real power flow
    vec2_Pa=nNodes+(nSwitches)+1:nNodes+(nSwitches)*2;
    s=size(f,1);
    store=store_+nNodes+nSwitches*4;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Pa,store,s,nSwitches);
    %%     
    vec2_Pb=nNodes+(nSwitches)*2+1:nNodes+(nSwitches)*3;
    s=size(f,1);
    store=store_+nNodes+nSwitches*4+nEdges;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Pb,store,s,nSwitches);
    %%     
    vec2_Pc=nNodes+(nSwitches)*3+1:nNodes+(nSwitches)*4;
    s=size(f,1);
    store=store_+nNodes+nSwitches*4+nEdges*2;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Pc,store,s,nSwitches);
    %%     
    %Linearzing product of x_ij and Q_ij for three phases for formulating the reactive power flow
    vec2_Qa=nNodes+nEdges*3+nSwitches*4+1:nNodes+nEdges*3+nSwitches*5;
    s=size(f,1);
    store=store_+nNodes+nSwitches*7+nEdges*3;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Qa,store,s,nSwitches);
    %%     
    vec2_Qb=nNodes+nEdges*3+nSwitches*5+1:nNodes+nEdges*3+nSwitches*6;
    s=size(f,1);
    store=store_+nNodes+nSwitches*7+nEdges*4;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Qb,store,s,nSwitches);
    %%     
    vec2_Qc=nNodes+nEdges*3+nSwitches*6+1:nNodes+nEdges*3+nSwitches*7;
    s=size(f,1);
    store=store_+nNodes+nSwitches*7+nEdges*5;
    [Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2_Qc,store,s,nSwitches);
     %% 
     %%   
    %Calculate cummulative active power of "active" nodes based on location of Source;   
     phase=2;
    [P1eq_a, PLeq_a,~]= Real_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    [Q1eq_a, QLeq_a,~]= Reactive_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    
     phase=4;
    [P1eq_b, PLeq_b,~]= Real_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    [Q1eq_b, QLeq_b,~]= Reactive_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    
     phase=6;
    [P1eq_c, PLeq_c,~]= Real_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    [Q1eq_c, QLeq_c,~]= Reactive_Power_Flow_Linear(Source,power,G1,Tie_Switches,phase);
    
    
    %[Volt1, Vp_A, Vq_A,Vp_B, Vq_B,Vp_C, Vq_C,vec2]=Voltage_Constraints_PhaseA(edges1,line,nNodes,nEdges,G1, Source,linepar,z_721,z_722,z_723,z_724); 
    z1=zeros(nNodes-1,nEdges);
    z0=zeros(nNodes-1,nSwitches);
    z2=zeros(nNodes-1,nNodes);
    z3=zeros(nNodes-1,nSwitches*3);
    z4=zeros(nNodes-1,nSwitches*2);
    Real_Power_Phase_A=[PLeq_a z0 z3 P1eq_a z1 z1 z3 z1 z1 z1 z2 z2 z2 z4 z4 z4];   
    Real_Power_Phase_B=[PLeq_b z0 z3 z1 P1eq_b z1 z3 z1 z1 z1 z2 z2 z2 z4 z4 z4];
    Real_Power_Phase_C=[PLeq_c z0 z3 z1 z1 P1eq_c z3 z1 z1 z1 z2 z2 z2 z4 z4 z4];
    
    Reactive_Power_Phase_A=[QLeq_a z0 z3 z1 z1 z1 z3 Q1eq_a z1 z1 z2 z2 z2 z4 z4 z4];   
    Reactive_Power_Phase_B=[QLeq_b z0 z3 z1 z1 z1 z3 z1 Q1eq_b z1 z2 z2 z2 z4 z4 z4];
    Reactive_Power_Phase_C=[QLeq_c z0 z3 z1 z1 z1 z3 z1 z1 Q1eq_c z2 z2 z2 z4 z4 z4];
    %% 
    %%    
       
    Aeq=[Real_Power_Phase_A;Real_Power_Phase_B;Real_Power_Phase_C];
    Aeq=[Aeq;Reactive_Power_Phase_A;Reactive_Power_Phase_B;Reactive_Power_Phase_C];
    beq=[zeros(1316,1);zeros(1316,1);zeros(1316,1);zeros(1316,1);zeros(1316,1);zeros(1316,1)];             
    %% 
    %%  
    z1=zeros(nEdges,nNodes);
    z2=zeros(nEdges,nSwitches);
    z3=zeros(nEdges,nSwitches*3);
    z4=zeros(nEdges,nSwitches*2);    
    [Volt1,Volt2,Vp_A, Vq_A,Vp_B, Vq_B,Vp_C, Vq_C,vec2,store_]=Voltage_Constraints_PhaseA(edges1,line,nEdges,G1,Source,swtch,nNodes,nSwitches);
    Voltage_A=[z1 z2 z3 Vp_A Vp_B Vp_C z3 Vq_A Vq_B Vq_C Volt2 z1 z1 Volt1 z4 z4];
    
    %Linearzing product of x_ij and V_ij for formulating the Voltage
    vec11=repmat(vec1',1,2)';
    vec11=vec11(:)';
    vec2=vec2+nNodes+nSwitches*7+nEdges*6;
    store=store_+nNodes*4+nSwitches*7+nEdges*6;
    [Aineq, bineq]=xij_times_Vij(Aineq,bineq, vec11,vec2,store,nSwitches,s);     %% 
    %%  
    %%
    [Volt1,Volt2,Vp_A, Vq_A,Vp_B, Vq_B,Vp_C, Vq_C,vec2,store_]=Voltage_Constraints_PhaseB(edges1,line,nEdges,G1,Source,swtch,nNodes,nSwitches);
    Voltage_B=[z1 z2 z3 Vp_A Vp_B Vp_C z3 Vq_A Vq_B Vq_C z1 Volt2 z1 z4 Volt1 z4];
    
    %Linearzing product of x_ij and V_ij for formulating the Voltage

    vec2=vec2+nNodes*2+nSwitches*7+nEdges*6;
    store=store_+nNodes*4+nSwitches*7+nEdges*6+nSwitches*2;
    [Aineq, bineq]=xij_times_Vij(Aineq,bineq, vec11,vec2,store,nSwitches,s);
    %% 
    %%     
    [Volt1,Volt2,Vp_A,Vq_A,Vp_B, Vq_B,Vp_C, Vq_C,vec2,store_]=Voltage_Constraints_PhaseC(edges1,line,nEdges,G1,Source,swtch,nNodes,nSwitches);
    Voltage_C=[z1 z2 z3 Vp_A Vp_B Vp_C z3 Vq_A Vq_B Vq_C z1 z1 Volt2 z4 z4 Volt1]; 
    
    %Linearzing product of x_ij and V_ij for formulating the Voltage
    vec2=vec2+nNodes*3+nSwitches*7+nEdges*6;
    store=store_+nNodes*4+nSwitches*7+nEdges*6+nSwitches*4;
    [Aineq, bineq]=xij_times_Vij(Aineq,bineq, vec11,vec2,store,nSwitches,s);
    
    %%
    %%
    
    %Feeder capacity constraints for Phase A, B and C. Real and Reactive
    %Power Capacity Constraints for MG are also included here
    
    [Aineq, bineq]=Feeder_Capacity(f,Aineq,bineq,nEdges);
    
    % Modeling Fault at different locations
    Aeq_switch=zeros(1,size(f,1));
    Aeq_switch(1,1317+154)=1;  % 154 index is for one particular line in Feeder-d   
    Aeq=[Aeq;Aeq_switch];
    beq=[beq;0];
    
    % Cyclic Constraints
    load Aineq_loop;
    load bineq_loop;
    Aineq=[Aineq; Aineq_loop];
    bineq=[bineq; bineq_loop'];
    
    Aeq=[Aeq;Voltage_A;Voltage_B;Voltage_C];
    beq=[beq;zeros(nEdges,1);zeros(nEdges,1);zeros(nEdges,1)];
    Aeq_voltage=zeros(1,size(f,1));
    Aeq_voltage(1,10504+1317)=1;
    Aeq_voltage(2,10504+1317+1317)=1;
    Aeq_voltage(3,10504+1317+1317+1317)=1;
    Aeq=[Aeq;Aeq_voltage];
    beq=[beq;1;1;1]; 
    
    %% 
    %%     
     
    lb1=zeros(nNodes+nSwitches,1);
    ub1=ones(nNodes+nSwitches,1);
    lb2=-100000*ones(size(f,1)-nNodes-nSwitches,1);
    ub2=10000*ones(size(f,1)-nNodes-nSwitches,1);
    lb=[lb1;lb2];
    ub=[ub1;ub2];
    
    %options = cplexoptimset;
    options.Display = 'on';
    size(f);
    for i=1:nNodes+nSwitches;
        ctype(i)='I';
    end
    for i=nNodes+(nSwitches)+1:size(f,1)
        ctype(i)='C';
    end
    
    options = cplexoptimset('MaxIter',900);
    
    disp ('The time taken for model creation in MATLAB is:');
    toc
    tic    
    % use this if it converges within time
    [x,fval] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype);
  
    % Use this function where it do not converges but actually has reached
    % the optimal value.
    
    %[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
   
    fprintf ('Number of nodes active in MG = %f\n', fval);
    disp ('Index =');
    %disp (x');
    active=find(x>0.5);
    store1=[];
       %Real Power
       a=nNodes+nSwitches*4; 
       c=nNodes+nSwitches*7+nEdges*3;
       [(1:nEdges)' edges1 x(a+1:a+nEdges) x(a+nEdges+1:a+nEdges*2) x(a+nEdges*2+1:a+nEdges*3)...
           x(c+1:c+nEdges) x(c+nEdges+1:c+nEdges*2) x(c+nEdges*2+1:c+nEdges*3)]
   
   %% 
   %%     
   %Finding what is the total amount of power supplied and the
   %loadshedding for each case
  F_a_real=[x(a+329) x(a+nEdges+329) x(a+nEdges*2+329)];
  F_b_real=[x(a+329*2) x(a+nEdges+329*2) x(a+nEdges*2+329*2)];
  F_c_real=[x(a+329*3) x(a+nEdges+329*3) x(a+nEdges*2+329*3)];
  F_d_real=[x(a+329*4) x(a+nEdges+329*4) x(a+nEdges*2+329*4)];

  
  F_a_reac=[x(c+329) x(c+nEdges+329) x(c+nEdges*2+329)];
  F_b_reac=[x(c+329*2) x(c+nEdges+329*2) x(c+nEdges*2+329*2)];
  F_c_reac=[x(c+329*3) x(c+nEdges+329*3) x(c+nEdges*2+329*3)];
  F_d_reac=[x(c+329*4) x(c+nEdges+329*4) x(c+nEdges*2+329*4)];
  
  Real_power=4*(1450.724+1499.657+1416.573);
  Reac_power=4*(432.783+444.933+421.489);
  
  DG_a_real=[x(a+1324) x(a+nEdges+1324) x(a+nEdges*2+1324)];   DG_a_reac=[x(c+1324) x(c+nEdges+1324) x(c+nEdges*2+1324)];
  DG_b_real=[x(a+1325) x(a+nEdges+1325) x(a+nEdges*2+1325)];   DG_b_reac=[x(c+1325) x(c+nEdges+1325) x(c+nEdges*2+1325)];
  DG_c_real=[x(a+1326) x(a+nEdges+1326) x(a+nEdges*2+1326)];   DG_c_reac=[x(c+1326) x(c+nEdges+1326) x(c+nEdges*2+1326)];
  DG_d_real=[x(a+1327) x(a+nEdges+1327) x(a+nEdges*2+1327)];   DG_d_reac=[x(c+1327) x(c+nEdges+1327) x(c+nEdges*2+1327)];

  
  Loadshed_Real=Real_power-(sum(F_a_real)+sum(F_b_real)+sum(F_c_real)+sum(F_d_real)+sum(DG_a_real)+sum(DG_b_real)+sum(DG_c_real)+sum(DG_d_real));
  Loadshed_Reac=Reac_power-(sum(F_a_reac)+sum(F_b_reac)+sum(F_c_reac)+sum(F_d_reac)+sum(DG_a_reac)+sum(DG_b_reac)+sum(DG_c_reac)+sum(DG_d_reac));
  Load_Shed_KVA=sqrt(Loadshed_Real^2+Loadshed_Reac^2);
  
  disp('   Total load to be shed during the restoration process is: ')
  disp('   .......................................................')
  disp(   Load_Shed_KVA);
  
  %% 
  %%
    disp('    The sectionalizing switches status for different feeders are: ');
    disp('      S.No     Feeder-a                    Feeder-b                       Feeder-c                    Feeder-d')
    disp('    ......||............................||............................||............................||................................');
    disp('    ......||............................||............................||............................||................................');


    result_sec=[(1:40)' swtch(1:40,:) x(1318:1317+40) swtch(41:80,:)-329 x(1318+40:1317+80) swtch(81:120,:)-329*2 x(1318+80:1317+120) ...
        swtch(121:160,:)-329*3 x(1318+120:1317+160)];
    disp (result_sec)
    
    result_tie=[(1:7)' swtch(165:171,:) x(1318+164:1317+171)];
    disp('....................................................');
    disp('....................................................');
    disp('The tie-switches status for the feeder are: ');
    disp( result_tie)
    disp('.....................................................');
    disp('.....................................................');
    
    result_DG=[(1:4)' swtch(172:175,:) x(1318+171:1317+175)];

    disp('The DG status for the feeder are: ');
    disp( result_DG)
    disp('.....................................................');
    disp('.....................................................');
    
    disp('The time taken for simulation in CPLEX is:')
    toc
    %%
    %%
       
     
    catch m
    disp(m.message);   
    end    

  end











   
 