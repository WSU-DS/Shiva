function [P1eq, P11eq,v1]= Real_Power_Flow_Linear(DER,power,G1,Tie_Switches,phase)
    v1=dfsearch(G1,DER);
    T= dfsearch(G1,DER,'edgetonew');
    fr=T(:,1);
    t=T(:,2);
    P1eq=zeros(36,36);
    P11eq=zeros(35,1317);
    insert=1;
    arrange= flipud(v1);
    for k=size(v1):-1:2        
        k ;
        a=v1(k);
         P1=power(a,phase);
         n=find(a==arrange);
         P1eq(insert,n)=1;
         P11eq(insert,a)=-1*P1;
         add=find(fr==v1(k));
          for m=1:size(add,1)  
              o=t(add(m));
              oo=find(arrange==o);
              P1eq(insert,oo)=-1;
          end
         insert=insert+1;          
    end 
    
    add_new= flipud(dfsearch(G1,DER,'edgetonew'));
    insert=1317;
    for k=1:size(Tie_Switches,1)
        a=find(Tie_Switches(k,1)==add_new(:,2));
        b=find(Tie_Switches(k,2)==add_new(:,2));
        P1eq(a,insert)=-1;
        P1eq(b,insert)=1;
        insert=insert+1;
    end
end
   