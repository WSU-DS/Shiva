    
function [Aineq, bineq]=edges_from_nodes(Aineq,bineq,f,nNodes,edges,nEdges)
    aineqXij=zeros(1,size(f,1));
    ind=nNodes+1;
    ins=1;
    for k=1:nEdges
        aineqXij(ins,ind)=1;
        aineqXij(ins,edges(k,1))=-1;
        ins=ins+1;
        ind=ind+1;
    end      
    ind=nNodes+1;
    for k=1:nEdges
        aineqXij(ins,ind)=1;
        aineqXij(ins,edges(k,2))=-1;
        ins=ins+1;
        ind=ind+1;
    end
%     ind=nNodes+1;
%     for k=1:nEdges
%         aineqXij(ins,ind)=-1;
%         aineqXij(ins,edges(k,1))=1;
%         aineqXij(ins,edges(k,2))=1;
%         ins=ins+1;
%         ind=ind+1;
%     end
     Aineq=[Aineq;aineqXij];
     bineq=[bineq;zeros(nEdges,1);zeros(nEdges,1)];
end