function cplexmilp
%% 
tic
    try
    clc;
    %lets keep 5 MGs in total
    f = ones(125*5, 1);
    ff=zeros(125*15,1);
    f=[f ;ff];
    nNodes=125;
    %Take the data from .txt files
    load 'branch.txt';
    load 'powerdata.txt';
    load 'linedata.txt'
    line=linedata;
    power=powerdata;
    edges=branch;
    size(edges);
    fr=edges(:,1);
    t=edges(:,2);
    G1 = graph(fr,t);
    plot(G1,'r','LineWidth',1);
    
    % Here we form ineq matrix for all DGs separately and combine them
    % later on. Ineq matrix means no two nodes can belong to same MG, child
    % and parent node constraints... the total number of constraints for
    % branch will be n-1. so the A1-A5 matrix is 124*125.
    
    
    mg=[4 26 46 86 67];
    %% 
    %Start With Mg1
    T1= dfsearch(G1,4,'edgetonew');
    fr=T1(:,1);
    t=T1(:,2);     
    v1=dfsearch(G1,4);
    size(v1);
    A1=formineq(G1,mg(1),nNodes);      
    %Calculate cummulative active power of "active" nodes based on DG at 4.
    %Later these values are used to calculate the voltage at the nodes with
    %critical loads
    P1eq=zeros(125,125);
    P11eq=zeros(125,125);
     insert=1;
     for k=size(v1):-1:1        
         a=v1(k);
         P1(a)=power(a,2);
         P1eq(insert,a)=1;
         P11eq(insert,a)=-1*P1(a);
         add=find(fr==v1(k));
         for m=1:size(add,1)  
              P1eq(insert,t(add(m)))=-1;
         end
         insert=insert+1;          
     end
     
    %Calculate cummulative reactive power of "active" nodes based on DG at 4.
    %Later these values are used to calculate the voltage at node with
    %critical loads 
    Q1eq=zeros(125,125);
    Q11eq=zeros(125,125);
     insert=1;
     for k=size(v1):-1:1        
         a=v1(k);
         Q1(a)=power(a,3);
         Q1eq(insert,a)=1;
         Q11eq(insert,a)=-1*Q1(a);
         add=find(fr==v1(k));
          for m=1:size(add,1)  
              Q1eq(insert,t(add(m)))=-1;
          end
         insert=insert+1;          
     end     
    %Calculate voltage equalities parameters for nodes based on DG at 4
    %V_j+1=V_j-(r*P_i+x*Q_i) 
 
    Volt1=zeros(125,125);
    Vp1=zeros(125,125);
    Vq1=zeros(125,125);
    n=size((v1),1)-1;
    Volt1(1,4)=1;
    for i=2:size(v1)        
        Volt1(i,fr(n))=1;
        Volt1(i,t(n))=-1;
         res=find (line(:,1)==fr(n));        
            for k=1:size(res)
                if line(res(k),2)==t(n)
                     r=line(res(k),3);
                     x=line(res(k),3);
                end
            end         
         res=find (line(:,2)==fr(n));
         for k=1:size(res)
           if line(res(k),1)==t(n)
                r=line(res(k),3);
                 x=line(res(k),3);
           end
         end 
        Vp1(i,t(n))=-r/100;
        Vq1(i,t(n))=-x/100;
        n=n-1;
    end     
    %% 
    %% 
    %Mg#2
    T2= dfsearch(G1,26,'edgetonew');
    fr=T2(:,1);
    t=T2(:,2);     
    v2=dfsearch(G1,26);
    size(v2);
    A2=formineq(G1,mg(2),nNodes);      
    %Calculate cummulative active power of "active" nodes based on DG at 26.
    %Later these values are used to calculate the voltage at the nodes with
    %critical loads
    P2eq=zeros(125,125);
    P22eq=zeros(125,125);
     insert=1;
     for k=size(v2):-1:1        
         a=v2(k);
         P1(a)=power(a,2);
         P2eq(insert,a)=1;
         P22eq(insert,a)=-1*P1(a);
         add=find(fr==v2(k));
         for m=1:size(add,1)  
              P2eq(insert,t(add(m)))=-1;
         end
         insert=insert+1;          
     end
     
    %Calculate cummulative reactive power of "active" nodes based on DG at 26.
    %Later these values are used to calculate the voltage at node with
    %critical loads 
    Q2eq=zeros(125,125);
    Q22eq=zeros(125,125);
     insert=1;
     for k=size(v2):-1:1        
         a=v2(k);
         Q1(a)=power(a,3);
         Q2eq(insert,a)=1;
         Q22eq(insert,a)=-1*Q1(a);
         add=find(fr==v2(k));
          for m=1:size(add,1)  
              Q2eq(insert,t(add(m)))=-1;
          end
         insert=insert+1;          
     end     
    %Calculate voltage equalities parameters for nodes based on DG at 26
    %V_j+1=V_j-(r*P_i+x*Q_i) 
   
    Volt2=zeros(125,125);
    Vp2=zeros(125,125);
    Vq2=zeros(125,125);
    n=size((v2),1)-1;
    Volt2(1,26)=1;
    for i=2:size(v2)        
        Volt2(i,fr(n))=1;
        Volt2(i,t(n))=-1;
        res=find (line(:,1)==fr(n));        
            for k=1:size(res)
                if line(res(k),2)==t(n)
                     r=line(res(k),3);
                     x=line(res(k),3);
                end
            end         
         res=find (line(:,2)==fr(n));
         for k=1:size(res)
           if line(res(k),1)==t(n)
                r=line(res(k),3);
                 x=line(res(k),3);
           end
         end 
        Vp2(i,t(n))=-r/100;
        Vq2(i,t(n))=-x/100;
        n=n-1;
    end     
    %% 
    %% 
    %MG#3
    T3= dfsearch(G1,46,'edgetonew');
    fr=T3(:,1);
    t=T3(:,2);     
    v3=dfsearch(G1,46);
    size(v3);
    A3=formineq(G1,mg(3),nNodes);      
    %Calculate cummulative active power of "active" nodes based on DG at 46.
    %Later these values are used to calculate the voltage at the nodes with
    %critical loads
    P3eq=zeros(125,125);
    P33eq=zeros(125,125);
    insert=1;
     for k=size(v3):-1:1        
         a=v3(k);
         P1(a)=power(a,2);
         P3eq(insert,a)=1;
         P33eq(insert,a)=-1*P1(a);
         add=find(fr==v3(k));
         for m=1:size(add,1)  
              P3eq(insert,t(add(m)))=-1;
         end
         insert=insert+1;          
     end
     
    %Calculate cummulative reactive power of "active" nodes based on DG at 46.
    %Later these values are used to calculate the voltage at node with
    %critical loads 
    Q3eq=zeros(125,125);
    Q33eq=zeros(125,125);
     insert=1;
     for k=size(v3):-1:1        
         a=v3(k);
         Q1(a)=power(a,3);
         Q3eq(insert,a)=1;
         Q33eq(insert,a)=-1*Q1(a);
         add=find(fr==v3(k));
          for m=1:size(add,1)  
              Q3eq(insert,t(add(m)))=-1;
          end
         insert=insert+1;          
     end     
    %Calculate voltage equalities parameters for nodes based on DG at 46
    %V_j+1=V_j-(r*P_i+x*Q_i) 
    
    Volt3=zeros(125,125);
    Vp3=zeros(125,125);
    Vq3=zeros(125,125);
    n=size((v3),1)-1;
    Volt3(1,46)=1;
    for i=2:size(v3)        
        Volt3(i,fr(n))=1;
        Volt3(i,t(n))=-1;
        res=find (line(:,1)==fr(n));        
            for k=1:size(res)
                if line(res(k),2)==t(n)
                     r=line(res(k),3);
                     x=line(res(k),3);
                end
            end         
         res=find (line(:,2)==fr(n));
         for k=1:size(res)
           if line(res(k),1)==t(n)
                r=line(res(k),3);
                 x=line(res(k),3);
           end
         end 
        Vp3(i,t(n))=-r/100;
        Vq3(i,t(n))=-x/100;
        n=n-1;
    end     
    %% 
    %% 
    %MG#4
    T4= dfsearch(G1,86,'edgetonew');
    fr=T4(:,1);
    t=T4(:,2);     
    v4=dfsearch(G1,86);
    size(v4);
    A4=formineq(G1,mg(4),nNodes);      
    %Calculate cummulative active power of "active" nodes based on DG at 86.
    %Later these values are used to calculate the voltage at the nodes with
    %critical loads
    P4eq=zeros(125,125);
    P44eq=zeros(125,125);
     insert=1;
     for k=size(v4):-1:1        
         a=v4(k);
         P1(a)=power(a,2);
         P4eq(insert,a)=1;
         P44eq(insert,a)=-1*P1(a);
         add=find(fr==v4(k));
         for m=1:size(add,1)  
              P4eq(insert,t(add(m)))=-1;
         end
         insert=insert+1;          
     end
     
    %Calculate cummulative reactive power of "active" nodes based on DG at 86.
    %Later these values are used to calculate the voltage at node with
    %critical loads 
    Q4eq=zeros(125,125);
    Q44eq=zeros(125,125);
     insert=1;
     for k=size(v4):-1:1        
         a=v4(k);
         Q1(a)=power(a,3);
         Q4eq(insert,a)=1;
         Q44eq(insert,a)=-1*Q1(a);
         add=find(fr==v4(k));
          for m=1:size(add,1)  
              Q4eq(insert,t(add(m)))=-1;
          end
         insert=insert+1;          
     end     
    %Calculate voltage equalities parameters for nodes based on DG at 86
    %V_j+1=V_j-(r*P_i+x*Q_i) 
  
    Volt4=zeros(125,125);
    Vp4=zeros(125,125);
    Vq4=zeros(125,125);
    n=size((v4),1)-1;
    Volt4(1,86)=1;
    for i=2:size(v4)        
        Volt4(i,fr(n))=1;
        Volt4(i,t(n))=-1;
        res=find (line(:,1)==fr(n));        
            for k=1:size(res)
                if line(res(k),2)==t(n)
                     r=line(res(k),3);
                     x=line(res(k),3);
                end
            end         
         res=find (line(:,2)==fr(n));
         for k=1:size(res)
           if line(res(k),1)==t(n)
                r=line(res(k),3);
                 x=line(res(k),3);
           end
         end 
        Vp4(i,t(n))=-r/100;
        Vq4(i,t(n))=-x/100;
        n=n-1;
    end     
    %% 
    %% 
    %MG#5
    T5= dfsearch(G1,67,'edgetonew');
    fr=T5(:,1);
    t=T5(:,2);     
    v5=dfsearch(G1,67);
    size(v5);
    A5=formineq(G1,mg(5),nNodes);      
    %Calculate cummulative active power of "active" nodes based on DG at 67.
    %Later these values are used to calculate the voltage at the nodes with
    %critical loads
    P5eq=zeros(125,125);
    P55eq=zeros(125,125);
     insert=1;
     for k=size(v5):-1:1        
         a=v5(k);
         P1(a)=power(a,2);
         P5eq(insert,a)=1;
         P55eq(insert,a)=-1*P1(a);
         add=find(fr==v5(k));
         for m=1:size(add,1)  
              P5eq(insert,t(add(m)))=-1;
         end
         insert=insert+1;          
     end
     
    %Calculate cummulative reactive power of "active" nodes based on DG at 67.
    %Later these values are used to calculate the voltage at node with
    %critical loads 
    Q5eq=zeros(125,125);
    Q55eq=zeros(125,125);
     insert=1;
     for k=size(v5):-1:1        
         a=v5(k);
         Q1(a)=power(a,3);
         Q5eq(insert,a)=1;
         Q55eq(insert,a)=-1*Q1(a);
         add=find(fr==v5(k));
          for m=1:size(add,1)  
              Q5eq(insert,t(add(m)))=-1;
          end
         insert=insert+1;          
     end     
    %Calculate voltage equalities parameters for nodes based on DG at 67
    %V_j+1=V_j-(r*P_i+x*Q_i) 

    Volt5=zeros(125,125);
    Vp5=zeros(125,125);
    Vq5=zeros(125,125);
    n=size((v5),1)-1;
    Volt5(1,67)=1;
    for i=2:size(v5)        
        Volt5(i,fr(n))=1;
        Volt5(i,t(n))=-1;
        res=find (line(:,1)==fr(n));        
            for k=1:size(res)
                if line(res(k),2)==t(n)
                     r=line(res(k),3);
                     x=line(res(k),3);
                end
            end         
         res=find (line(:,2)==fr(n));
         for k=1:size(res)
           if line(res(k),1)==t(n)
                r=line(res(k),3);
                 x=line(res(k),3);
           end
         end 
        Vp5(i,t(n))=-r/100;
        Vq5(i,t(n))=-x/100;
        n=n-1;
    end     
    %% 
    %%   
    %Forming the Aineq and Bineq
    fill_zero=zeros(124,125);
     A1=[A1 fill_zero fill_zero fill_zero fill_zero];
     A2=[fill_zero A2 fill_zero fill_zero fill_zero];
     A3=[fill_zero fill_zero A3 fill_zero fill_zero];
     A4=[fill_zero fill_zero fill_zero A4 fill_zero];
     A5=[fill_zero fill_zero fill_zero fill_zero A5];    
     A6=zeros(125,125*5); 
    
    %inequality for each node to belong to either of MG or none of them
    for i=1:125
        if i~=4&&i~=26&&i~=67&&i~=46&&i~=86&&i~=17&&i~=9&&i~=30&&i~=66&&i~=101&&i~=83&&i~=79&&i~=27&&i~=37
          A6(i,i)=1;
          A6(i,i+125)=1;
          A6(i,i+125*2)=1;
          A6(i,i+125*3)=1;
          A6(i,i+125*4)=1;
        end
    end
    Aineq=[A3;A2;A1;A4;A5;A6];
    ad=zeros(745, 125*15);
    Aineq=[Aineq ad];
    size(Aineq);
    bineq1=zeros(124*5,1);
    bineq2=ones(125,1);
    bineq=[bineq1;bineq2];
    size(bineq);  
    
    %% 
    %% 
    %Adding of branch constraints means a branch may be fail during a
    %disaster and can't serve load during the optimization process through
    %it. For example we consider line 18-13 is loss during an event. This
    %means 18 and 13 cannot belong to same MG. 
     br1=zeros(1,125*20); br2=zeros(1,125*20); br3=zeros(1,125*20);
     br4=zeros(1,125*20); br5=zeros(1,125*20);
     br1(1,18)=1; br1(1,13)=1; br2(1,18+125)=1; br2(1,13+125)=1;
     br3(1,18+125*2)=1; br3(1,13+125*2)=1; br4(1,18+125*3)=1;
     br4(1,13+125*3)=1; br5(1,18+125*4)=1; br5(1,13+125*4)=1;
     size(Aineq);
     Aineq=[Aineq;br1]; Aineq=[Aineq;br2]; Aineq=[Aineq;br3];
     Aineq=[Aineq;br4]; Aineq=[Aineq;br5];
     size(Aineq);
     size(bineq);
     bineq=[bineq;1;1;1;1;1];
   
    %% 
    %%   
    %Forming Aeq matrix    
    Aeq=zeros(14,125*5); 
    %equality constraints for CL to be served by only one MG
    %List of Critical Loads
    CL=[9 17 83 79 101 66 27 30 37]';
    for i=1:size(CL)
        Aeq(i,CL(i))=1;
        Aeq(i,CL(i)+125)=1;
        Aeq(i,CL(i)+125*2)=1;
        Aeq(i,CL(i)+125*3)=1;
        Aeq(i,CL(i)+125*4)=1;
    end
    % Node where DG is installed is always one.
    Aeq(10,4)=1; 
    Aeq(11,26+125)=1;
    Aeq(12,46+125*2)=1;
    Aeq(13,86+125*3)=1;
    Aeq(14,60+125*4)=1;
    
    add_for_PV=zeros(14,125*15);
    Aeq=[Aeq add_for_PV];
    z=zeros(125,125);
    size(Aeq)
    mg1p=[P11eq z z z z P1eq z z z z z z z z z z z z z z ]; 
    mg2p=[z P22eq z z z z P2eq z z z z z z z z z z z z z ]; 
    mg3p=[z z P33eq z z z z P3eq z z z z z z z z z z z z ];   
    mg4p=[z z z P44eq z z z z P4eq z z z z z z z z z z z ]; 
    mg5p=[z z z z P55eq z z z z P5eq z z z z z z z z z z ];
    
    mg1q=[Q11eq z z z z z z z z z Q1eq z z z z z z z z z];   
    mg2q=[z Q22eq z z z z z z z z z Q2eq z z z z z z z z]; 
    mg3q=[z z Q33eq z z z z z z z z z Q3eq z z z z z z z];   
    mg4q=[z z z Q44eq z z z z z z z z z Q4eq z z z z z z]; 
    mg5q=[z z z z Q55eq z z z z z z z z z Q5eq z z z z z];
    
    mg1V=[z z z z z Vp1 z z z z Vq1 z z z z Volt1 z z z z];
    mg2V=[z z z z z z Vp2 z z z z Vq2 z z z z Volt2 z z z];
    mg3V=[z z z z z z z Vp3 z z z z Vq3 z z z z Volt3 z z];
    mg4V=[z z z z z z z z Vp4 z z z z Vq4 z z z z Volt4 z];
    mg5V=[z z z z z z z z z Vp5 z z z z Vq5 z z z z Volt5];
        
    beq=[1;1;1;1;1;1;1;1;1;1;1;1;1;1];
    %The following equality matrix contains P and V constraints altogether
    Aeq=[Aeq;mg1p;mg2p;mg3p;mg4p;mg5p;mg1q;mg2q;mg3q;mg4q;mg5q;mg1V;mg2V;mg3V;mg4V;mg5V];
    badd=zeros(125*10,1);
    beq=[beq;badd;1.05;zeros(124,1);1.05;zeros(124,1);1.05;zeros(124,1);1.05;zeros(124,1);1.05;zeros(124,1)];
    size(Aeq);
    size(beq);
    lb=zeros(125*2,1);
    l1=ones(125*18,1)*-2;
    lb=[lb;l1];
    ub=ones(125*5,1);
    ub_add=1000*ones(125*15,1);
    ub=[ub;ub_add]; 
 
    %% 
    %Solver
    options = cplexoptimset;
    options.Display = 'on';
    for i=1:125*5
        ctype(i)='I';
    end
     for i=125*5+1:125*20
        ctype(i)='C';
    end
    [x,fval] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype);
    fprintf ('Number of nodes active in MG = %f\n', fval);
     disp ('Index =');
     disp (x');
    active=find(x==1);
    size(active);
    mg1=[];
    mg2=[];
    mg3=[];
    mg4=[];
    mg5=[];
    for i=1:size(active)
        if active(i)<=125
            mg1=[mg1 active(i)];
        end
        if active(i)>125 && active(i)<=125*2
            active(i)=active(i)-125;
           mg2=[mg2 active(i)];
        end
         if active(i)>125*2&&active(i)<125*3
             active(i)=active(i)-250;
            mg3=[mg3 active(i)];
         end
        if active(i)>125*3 && active(i)<=125*4
            active(i)=active(i)-125*3;
           mg4=[mg4 active(i)];
        end
         if active(i)>125*4&&active(i)<=125*5
             active(i)=active(i)-125*4;
            mg5=[mg5 active(i)];
        end
    end
    mg1
    mg2
    mg3
    mg4
    mg5
    
    catch m
    disp(m.message);   
    end  
    toc
end