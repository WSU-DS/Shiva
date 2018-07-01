
function[Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2,store,s,nSwitches)
    A_low=-2000;
    A_high=2000;
    % For min{0,A_low}<=z
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=-1;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; -A_low*ones(nSwitches,1)];
    %% 
    %% 
    % For z<=A_high
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=1;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; A_high*ones(nSwitches,1)];
    %% 
    %%    
    % For zp1<=A_high*x;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=1;
        zp1(p,vec1(p))=-A_high;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; zeros(nSwitches,1)];

    %% 
    %%    
    % For A_low*x<=z
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=-1;
        zp1(p,vec1(p))=A_low;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; zeros(nSwitches,1)];
    %% 
    %%     
    % For zp1<=A-(1-x)A_low;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=1;        
        zp1(p,vec2(p))=-1;
        zp1(p,vec1(p))=-A_low;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; -A_low*ones(nSwitches,1)];
    %% 
    %%     
    % For zp1>=A-(1-x)A_high;
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=-1;        
        zp1(p,vec2(p))=1;
        zp1(p,vec1(p))=A_high;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; A_high*ones(nSwitches,1)];
    %% 
    %% 
    %For z<=A+(1-x)A_high
    zp1=zeros(1,s);
    p=1;
    for i=1:nSwitches
        zp1(p,store(i))=1;        
        zp1(p,vec2(p))=-1;
        zp1(p,vec1(p))=A_high;
        p=p+1;
    end
    Aineq=[Aineq;zp1];
    bineq=[bineq; A_high*ones(nSwitches,1)];   
end




% function[Aineq, bineq]=Product_bound_AA_(Aineq,bineq, vec1,vec2,store,s,nSwitches)
%     A_high=10000;
%     A_low=-10000;
%     % For min{0,A_low}<=z
%     zp1=zeros(1,s);
%     p=1;
%     %%    
%     % For zp1<=A_high*x;
%     for i=1:nSwitches
%         zp1(p,store(i))=1;
%         zp1(p,vec1(p))=-A_high;
%         p=p+1;
%     end
%     Aineq=[Aineq;zp1];
%     bineq=[bineq; zeros(nSwitches,1)];
% 
%     %% 
%     %%    
%     % For A_low*x<=z
%     zp1=zeros(1,s);
%     p=1;
%     for i=1:nSwitches
%         zp1(p,store(i))=-1;
%         zp1(p,vec1(p))=A_low;
%         p=p+1;
%     end
%     Aineq=[Aineq;zp1];
%     bineq=[bineq; zeros(nSwitches,1)];
%     %% 
%     %%       
% end
% 
% 
% 
% 
