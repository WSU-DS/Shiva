function [Q1eq, Q11eq,v1]= Reactive_Power_Flow_Linear(DER,power,G1,Tie_Switches,phase)
    v1=dfsearch(G1,DER);
    T= dfsearch(G1,DER,'edgetonew');
    fr=T(:,1);
    t=T(:,2);
    Q1eq=zeros(36,36);
    Q11eq=zeros(35,1317);
    insert=1;
    arrange= flipud(v1);
    for k=size(v1):-1:2        
         a=v1(k);
         Q1=power(a,phase+1);
         n=find(a==arrange);
         Q1eq(insert,n)=1;
         Q11eq(insert,a)=-1*Q1;
         add=find(fr==v1(k));
          for m=1:size(add,1)  
              o=t(add(m));
              oo=find(arrange==o);
              Q1eq(insert,oo)=-1;
          end
         insert=insert+1;          
    end    
    
    
    add_new= flipud(dfsearch(G1,DER,'edgetonew'));
    %%
    %%
    %Now insert new variables for inlcuding power flow as per new tie switches
    insert=1317;
    for k=1:size(Tie_Switches,1)
        a=find(Tie_Switches(k,1)==add_new(:,2));
        b=find(Tie_Switches(k,2)==add_new(:,2));
        Q1eq(a,insert)=-1;
        Q1eq(b,insert)=1;
        insert=insert+1;
    end
end
    