clear all
clc

% This code and attached files is about extracting the parameters for IEEE 8500 node test case.
% Originally looking into Excel files with line.dss and bus_coordinate, we
% extract the information to be used in graph theory related applications.
% The attached files are called upon for building the required graph.

%% 
%% 
%Following are new modified things for regulator and capacitor..

%regulator 2, 3, and 4 of the test cases..
%232	717   --> regxfmr_190-7361-190-7361
%1003	561   --> regxfmr_190-8581-190-8581
%2514	1132  --> regxfmr_190-8593-190-8593

% These are capacitors defined as single phases. However, there is a three
% phase capacitor at node R18242 as well.

% 939-1619  -->L2823592-R20185 Linecode=3PH_H-397_ACSR397_ACSR397_ACSR2/0_ACSR  
% 1797-1972 -->R42247-Q16483   Linecode=3PH_H-397_ACSR397_ACSR397_ACSR2/0_ACSR  
% 451-447   -->R42246-Q16642   Linecode=3PH_H-2/0_ACSR2/0_ACSR2/0_ACSR2_ACSR    


%% 
%% 
load lines
fid2 = fopen('from.txt');
fr=importdata('from.txt');
fid3 = fopen('to.txt');
to=importdata('to.txt');
load 'bus_coor1.txt'
Bus=bus_coor1;
X=Bus(:,2); Y=Bus(:,3);
plot(X,Y,'.r');

l = final_lines;
hold on
for k=1:2520
    plot([X(l(k,1)) X(l(k,2))],[Y(l(k,1)) Y(l(k,2))],'k');
    hold on
end

% Draw the location for regulators 
draw=[232 717 1003 561 2514 1132];
for k =1:6
     text(X(draw(k)),Y(draw(k)),num2str(draw(k)),'fontsize',9);   
end

% Draw location for capacitor connection left behing in final_lines file
draw=[939 1619 1797 1972 451 447];
for k =1:6
     text(X(draw(k)),Y(draw(k)),num2str(draw(k)),'fontsize',9);   
end

% if you want all numbers to be displayed then comment the following.
% for k =1:2514
%      text(X(k),Y(k),num2str(k),'fontsize',9);   
% end

%% 
%% 
% This part of the code is for line parameters extraction of OpenDSS IEEE 8500 node test case. The
% R and X matrix is extracted for getting the topology and line parameters.

% - Start OpenDSS
DSSObj=actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0)
    disp('Unable to start openDSS Engine');
    return
end
DSSText = DSSObj.Text;
DSSCircuit=DSSObj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;

DSSText.Command='Compile (C:\Users\spoudel\Dropbox\Ph.D_WSU\Research\ADMS\8500_restoration\Extract_Nodes\8500_Node_Tom\Master.dss)';
DSSLines=DSSObj.ActiveCircuit.Lines;

fid1   = fopen('Line_names.txt');
name=importdata('Line_names.txt');
% Extract parameters for 2511 lines. There are few more lines added for regulator and capacitor...
for k=1:2511
    line_name=char(name(k));
    if k==47
        line_name='293471';
    end
    DSSLines.Name=char(line_name);
    MyLengthVariable(k) = DSSLines.Length;
    MyRmatrix{k} = DSSLines.Rmatrix;
    MyXmatrix{k} = DSSLines.Xmatrix;
end
% These data are in ohm/km. You can just multiply by length to get line r
% and x...
MyRmatrix{k}
MyXmatrix{k}

%% 
%% 

































% 
% G=graph(final_lines(:,1),final_lines(:,2));
% plot (G)
% for k=1:2514
%     N{k} = neighbors(G,k);
%     a=isempty(N{k});
%     if a==1
%         index(m)=k;
%         m=m+1;
%     end
% end
%% 
%% 
% for m=1:size(final_lines,1)
%     if final_lines(m,1)==0
%         draw(k)=final_lines(m,2);
%         index(p)=m;
%         k=k+1;p=p+1;
%     end    
%     if final_lines(m,2)==0
%         draw(k)=final_lines(m,1);
%         index(p)=m;
%         k=k+1; p=p+1;
%     end    
%     if final_lines(m,1)==0&&final_lines(m,2)==0
%         missing(n)=m;
%         %index(p)=m;
%         n=n+1; %p=p+1;
%     end    
% end
% %% 
% %% 
% for k=1:size(final_lines,1)
%     if final_lines(k,1)==0;
%         final_lines(k,2)=0;
%     end
%      if final_lines(k,2)==0;
%         final_lines(k,1)=0;
%      end
% end
%l = final_lines(any(final_lines,2),:);
%l = final_lines;
%hold on

%for k=1:2520
 %   plot([X(l(k,1)) X(l(k,2))],[Y(l(k,1)) Y(l(k,2))],'k');
  %  hold on
%end
%draw=draw(draw~=0);
% for k=1:2510
%     plot(l(:,1),l(:,2),'k');
%     hold on
% end



