
function[Aineq, bineq]=Feeder_Capacity(f,Aineq,bineq,nEdges)
       
    nNodes=1317;
    nSwitches=175;

    % Change these capacities if you want to simulate with/without MGs
    real_P=2000;
    real_PM1=1500;
    real_PM2=700;
    real_PM3=700;
    real_PM4=700;
    
    a=nNodes+nSwitches+nSwitches*3; % Real Power Flow starts from here
   
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,a+329)=1;
    Aineq_feeder(2,a+658)=1;
    Aineq_feeder(3,a+987)=1;
    Aineq_feeder(4,a+1316)=1;
    Aineq_feeder(5,a+1324)=1;
    Aineq_feeder(6,a+1325)=1;
    Aineq_feeder(7,a+1326)=1;
    Aineq_feeder(8,a+1327)=1;
    
    bineq_feeder=[real_P;real_P;real_P;real_P;real_PM1;real_PM2;real_PM3;real_PM4];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
    
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,a+nEdges+329)=1;
    Aineq_feeder(2,a+nEdges+658)=1;
    Aineq_feeder(3,a+nEdges+987)=1;
    Aineq_feeder(4,a+nEdges+1316)=1;
    Aineq_feeder(5,a+nEdges+1324)=1;
    Aineq_feeder(6,a+nEdges+1325)=1;
    Aineq_feeder(7,a+nEdges+1326)=1;
    Aineq_feeder(8,a+nEdges+1327)=1;
    
    bineq_feeder=[real_P;real_P;real_P;real_P;real_PM1;real_PM2;real_PM3;real_PM4];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
    
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,a+nEdges*2+329)=1;
    Aineq_feeder(2,a+nEdges*2+658)=1;
    Aineq_feeder(3,a+nEdges*2+987)=1;
    Aineq_feeder(4,a+nEdges*2+1316)=1;  
    Aineq_feeder(5,a+nEdges*2+1324)=1;
    Aineq_feeder(6,a+nEdges*2+1325)=1;
    Aineq_feeder(7,a+nEdges*2+1326)=1;
    Aineq_feeder(8,a+nEdges*2+1327)=1;
    
    bineq_feeder=[real_P;real_P;real_P;real_P;real_PM1;real_PM2;real_PM3;real_PM4];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
    
    
    b=a+nEdges*3+nSwitches*3; % Reactive power Flow starts from here
    
    
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,b+329)=1;
    Aineq_feeder(2,b+658)=1;
    Aineq_feeder(3,b+987)=1;
    Aineq_feeder(4,b+1316)=1;
    bineq_feeder=[1000;1000;1000;1000];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
    
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,b+nEdges+329)=1;
    Aineq_feeder(2,b+nEdges+658)=1;
    Aineq_feeder(3,b+nEdges+987)=1;
    Aineq_feeder(4,b+nEdges+1316)=1;
    bineq_feeder=[1000;1000;1000;1000];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
    
    Aineq_feeder=zeros(1,size(f,1));
    Aineq_feeder(1,b+nEdges*2+329)=1;
    Aineq_feeder(2,b+nEdges*2+658)=1;
    Aineq_feeder(3,b+nEdges*2+987)=1;
    Aineq_feeder(4,b+nEdges*2+1316)=1;
    bineq_feeder=[1000;1000;1000;1000];
    Aineq=[Aineq;Aineq_feeder];
    bineq=[bineq;bineq_feeder];
end