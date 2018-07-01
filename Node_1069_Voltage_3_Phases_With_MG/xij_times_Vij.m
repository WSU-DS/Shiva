
function[Aineq, bineq]=xij_times_Vij(Aineq,bineq, vec1,vec2,store,nSwitches,s)
% For zp1>=0;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches*2
        zp1(p,store(i))=-1;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; zeros(nSwitches*2,1)];
    %% 
    %%     
    % For zp1<=M*x_ij;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches*2
        zp1(p,store(i))=1;
        zp1(p,vec1(p))=-1000;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; zeros(nSwitches*2,1)];
    %% 
    %%     
    % For zp1<=V_i;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches*2
        zp1(p,store(i))=1;        
        zp1(p,vec2(p))=-1;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; zeros(nSwitches*2,1)];
    %% 
    %% 
    % For zp1>=V_i-(1-x_ij)M;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches*2
        zp1(p,store(i))=-1;        
        zp1(p,vec2(p))=1;
        zp1(p,vec1(p))=1000;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; 1000*ones(nSwitches*2,1)];   
end