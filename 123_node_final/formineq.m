function A= formineq(G1,mg,nNodes)
    T1=dfsearch(G1,mg,'edgetonew');
    fr=T1(:,1);
    t=T1(:,2);
    v1=dfsearch(G1,mg);
    size(v1);
    A=zeros(size(v1,1)-1,nNodes);
    s=0;
    for i=1:size(fr)
       from=find(fr==v1(i));      
       for k=1:size(from,1)
            s=s+1;
            A(s,v1(i))=-1; 
            A(s,t(from(k)))=1;
       end              
    end
end