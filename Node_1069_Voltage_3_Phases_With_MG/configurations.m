
function[r_aa,x_aa,r_ab,x_ab,r_ac,x_ac]=configurations(code)
Z_21=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.2355 + 1.0421i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_22=[0.0657 + 0.9096i   0.0288 + 0.3035i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0657 + 0.9096i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0288 + 0.3035i   0.0657 + 0.9096i];

Z_23=[0.2355 + 1.0421i   0.0695 + 0.3373i   0.0695 + 0.3373i
   0.0695 + 0.3373i   0.2355 + 1.0421i   0.0695 + 0.3373i
   0.0695 + 0.3373i   0.0695 + 0.3373i   0.2355 + 1.0421i];

Z_24=[0.2355 + 1.0421i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_25=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.3805 + 1.2414i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_26=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   1.7617 + 1.3083i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_27=[ 0.0657 + 0.9096i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_28=[0.0657 + 0.9096i   0.0288 + 0.3035i   0.0000 + 0.0000i
   0.0288 + 0.3035i   0.0657 + 0.9096i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_29=[0.0657 + 0.9096i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_30=[1.7617 + 1.3083i   0.2617 + 0.5315i   0.2617 + 0.5315i
   0.2617 + 0.5315i   1.7617 + 1.3083i   0.2617 + 0.5315i
   0.2617 + 0.5315i   0.2617 + 0.5315i   1.7617 + 1.3083i];

Z_31=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   1.1891 + 1.1525i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_32=[0.0657 + 0.9096i   0.0288 + 0.3035i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0657 + 0.9096i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0288 + 0.3035i   0.0657 + 0.9096i];

Z_33=[ 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   1.1737 + 1.1933i];

Z_34=[1.1891 + 1.1525i   0.2361 + 0.4416i   0.2361 + 0.4416i
   0.2361 + 0.4416i   1.1891 + 1.1525i   0.2361 + 0.4416i
   0.2361 + 0.4416i   0.2361 + 0.4416i   1.1891 + 1.1525i];

Z_35=[0.1046 + 0.9572i   0.0389 + 0.3155i   0.0389 + 0.3155i
   0.0389 + 0.3155i   0.1046 + 0.9572i   0.0389 + 0.3155i
   0.0389 + 0.3155i   0.0389 + 0.3155i   0.1046 + 0.9572i];

Z_36=[];

Z_37=[0.0657 + 0.9096i   0.0288 + 0.3035i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0657 + 0.9096i   0.0288 + 0.3035i
   0.0288 + 0.3035i   0.0288 + 0.3035i   0.0657 + 0.9096i];

Z_46=[ 0.0586 + 0.8974i   0.0268 + 0.3004i   0.0268 + 0.3004i
   0.0268 + 0.3004i   0.0586 + 0.8974i   0.0268 + 0.3004i
   0.0268 + 0.3004i   0.0268 + 0.3004i   0.0586 + 0.8974i];

Z_47=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0586 + 0.8974i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_38=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   1.3755 + 0.5168i];

Z_39=[1.2714 + 0.5160i   0.1434 - 0.0538i   0.1434 - 0.0538i
   0.1434 - 0.0538i   1.2714 + 0.5160i   0.1434 - 0.0538i
   0.1434 - 0.0538i   0.1434 - 0.0538i   1.2714 + 0.5160i];

Z_40=[ 1.3755 + 0.5168i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_41=[ 1.8450 + 0.7989i   0.3348 + 0.0853i   0.3348 + 0.0853i
   0.3348 + 0.0853i   1.8450 + 0.7989i   0.3348 + 0.0853i
   0.3348 + 0.0853i   0.3348 + 0.0853i   1.8450 + 0.7989i];

Z_42=[0.8147 + 0.5793i   0.1894 - 0.0011i   0.1894 - 0.0011i
   0.1894 - 0.0011i   0.8147 + 0.5793i   0.1894 - 0.0011i
   0.1894 - 0.0011i   0.1894 - 0.0011i   0.8147 + 0.5793i];

Z_43=[0.2825 + 0.6900i   0.1581 + 0.1442i   0.1581 + 0.1442i
   0.1581 + 0.1442i   0.2825 + 0.6900i   0.1581 + 0.1442i
   0.1581 + 0.1442i   0.1581 + 0.1442i   0.2825 + 0.6900i];

Z_44=[ 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   1.3755 + 0.5168i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i];

Z_45=[0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   1.9295 + 0.9188i];
r_aa=0;x_aa=0;r_ab=0;x_ab=0;r_ac=0;x_ac=0;
 switch code
            case 21
                r_aa=real(Z_21(1,1));  x_aa=imag(Z_21(1,1));
                r_ab=real(Z_21(1,2));  x_ab=imag(Z_21(1,2));
                r_ac=real(Z_21(1,3));  x_ac=imag(Z_21(1,3));
                
            case 22
                r_aa=real(Z_22(1,1));  x_aa=imag(Z_22(1,1));
                r_ab=real(Z_22(1,2));  x_ab=imag(Z_22(1,2));
                r_ac=real(Z_22(1,3));  x_ac=imag(Z_22(1,3));
                
            case 23
                r_aa=real(Z_23(1,1));  x_aa=imag(Z_23(1,1));
                r_ab=real(Z_23(1,2));  x_ab=imag(Z_23(1,2));
                r_ac=real(Z_23(1,3));  x_ac=imag(Z_23(1,3));
                
            case 24
                r_aa=real(Z_24(1,1));  x_aa=imag(Z_24(1,1));
                r_ab=real(Z_24(1,2));  x_ab=imag(Z_24(1,2));
                r_ac=real(Z_24(1,3));  x_ac=imag(Z_24(1,3));
                
           case 25
                r_aa=real(Z_25(1,1));  x_aa=imag(Z_25(1,1));
                r_ab=real(Z_25(1,2));  x_ab=imag(Z_25(1,2));
                r_ac=real(Z_25(1,3));  x_ac=imag(Z_25(1,3));
                
            case 26
                r_aa=real(Z_26(1,1));  x_aa=imag(Z_26(1,1));
                r_ab=real(Z_26(1,2));  x_ab=imag(Z_26(1,2));
                r_ac=real(Z_26(1,3));  x_ac=imag(Z_26(1,3));
                
            case 27
                r_aa=real(Z_27(1,1));  x_aa=imag(Z_27(1,1));
                r_ab=real(Z_27(1,2));  x_ab=imag(Z_27(1,2));
                r_ac=real(Z_27(1,3));  x_ac=imag(Z_27(1,3));
                
            case 28
                r_aa=real(Z_28(1,1));  x_aa=imag(Z_28(1,1));
                r_ab=real(Z_28(1,2));  x_ab=imag(Z_28(1,2));
                r_ac=real(Z_28(1,3));  x_ac=imag(Z_28(1,3));
                
            case 29
                r_aa=real(Z_29(1,1));  x_aa=imag(Z_29(1,1));
                r_ab=real(Z_29(1,2));  x_ab=imag(Z_29(1,2));
                r_ac=real(Z_29(1,3));  x_ac=imag(Z_29(1,3));
                
            case 30
                r_aa=real(Z_30(1,1));  x_aa=imag(Z_30(1,1));
                r_ab=real(Z_30(1,2));  x_ab=imag(Z_30(1,2));
                r_ac=real(Z_30(1,3));  x_ac=imag(Z_30(1,3));
                
            case 31
                r_aa=real(Z_31(1,1));  x_aa=imag(Z_31(1,1));
                r_ab=real(Z_31(1,2));  x_ab=imag(Z_31(1,2));
                r_ac=real(Z_31(1,3));  x_ac=imag(Z_31(1,3));
                
            case 32
                r_aa=real(Z_32(1,1));  x_aa=imag(Z_32(1,1));
                r_ab=real(Z_32(1,2));  x_ab=imag(Z_32(1,2));
                r_ac=real(Z_32(1,3));  x_ac=imag(Z_32(1,3));
                
            case 33
                r_aa=real(Z_33(1,1));  x_aa=imag(Z_33(1,1));
                r_ab=real(Z_33(1,2));  x_ab=imag(Z_33(1,2));
                r_ac=real(Z_33(1,3));  x_ac=imag(Z_33(1,3));
                
            case 34
                r_aa=real(Z_34(1,1));  x_aa=imag(Z_34(1,1));
                r_ab=real(Z_34(1,2));  x_ab=imag(Z_34(1,2));
                r_ac=real(Z_34(1,3));  x_ac=imag(Z_34(1,3));
                
            case 35
                r_aa=real(Z_35(1,1));  x_aa=imag(Z_35(1,1));
                r_ab=real(Z_35(1,2));  x_ab=imag(Z_35(1,2));
                r_ac=real(Z_35(1,3));  x_ac=imag(Z_35(1,3));
                
            case 36
                r_aa=real(Z_36(1,1));  x_aa=imag(Z_36(1,1));
                r_ab=real(Z_36(1,2));  x_ab=imag(Z_36(1,2));
                r_ac=real(Z_36(1,3));  x_ac=imag(Z_36(1,3));
                
            case 37
                r_aa=real(Z_37(1,1));  x_aa=imag(Z_37(1,1));
                r_ab=real(Z_37(1,2));  x_ab=imag(Z_37(1,2));
                r_ac=real(Z_37(1,3));  x_ac=imag(Z_37(1,3));
                
            case 38
                r_aa=real(Z_38(1,1));  x_aa=imag(Z_38(1,1));
                r_ab=real(Z_38(1,2));  x_ab=imag(Z_38(1,2));
                r_ac=real(Z_38(1,3));  x_ac=imag(Z_38(1,3));
                
            case 39
                r_aa=real(Z_39(1,1));  x_aa=imag(Z_39(1,1));
                r_ab=real(Z_39(1,2));  x_ab=imag(Z_39(1,2));
                r_ac=real(Z_39(1,3));  x_ac=imag(Z_39(1,3));
                
               
            case 40
                r_aa=real(Z_40(1,1));  x_aa=imag(Z_40(1,1));
                r_ab=real(Z_40(1,2));  x_ab=imag(Z_40(1,2));
                r_ac=real(Z_40(1,3));  x_ac=imag(Z_40(1,3));
                
            case 41
                r_aa=real(Z_41(1,1));  x_aa=imag(Z_41(1,1));
                r_ab=real(Z_41(1,2));  x_ab=imag(Z_41(1,2));
                r_ac=real(Z_41(1,3));  x_ac=imag(Z_41(1,3));
                
            case 42
                r_aa=real(Z_42(1,1));  x_aa=imag(Z_42(1,1));
                r_ab=real(Z_42(1,2));  x_ab=imag(Z_42(1,2));
                r_ac=real(Z_42(1,3));  x_ac=imag(Z_42(1,3));
                
            case 43
                r_aa=real(Z_43(1,1));  x_aa=imag(Z_43(1,1));
                r_ab=real(Z_43(1,2));  x_ab=imag(Z_43(1,2));
                r_ac=real(Z_43(1,3));  x_ac=imag(Z_43(1,3));
                
            case 44
                r_aa=real(Z_44(1,1));  x_aa=imag(Z_44(1,1));
                r_ab=real(Z_44(1,2));  x_ab=imag(Z_44(1,2));
                r_ac=real(Z_44(1,3));  x_ac=imag(Z_44(1,3));
                
            case 45
                r_aa=real(Z_45(1,1));  x_aa=imag(Z_45(1,1));
                r_ab=real(Z_45(1,2));  x_ab=imag(Z_45(1,2));
                r_ac=real(Z_45(1,3));  x_ac=imag(Z_45(1,3));
                
            case 46
                r_aa=real(Z_46(1,1));  x_aa=imag(Z_46(1,1));
                r_ab=real(Z_46(1,2));  x_ab=imag(Z_46(1,2));
                r_ac=real(Z_46(1,3));  x_ac=imag(Z_46(1,3));
                
            case 47
                r_aa=real(Z_47(1,1));  x_aa=imag(Z_47(1,1));
                r_ab=real(Z_47(1,2));  x_ab=imag(Z_47(1,2));
                r_ac=real(Z_47(1,3));  x_ac=imag(Z_47(1,3));
            
 end
end



































% clear all
% clc
% load Z_mat_UG;
% load Z_matrix_OH;
% load 'Line_OH.txt'; ov_line=Line_OH;
% load 'Line_UG.txt'; ug_line=Line_UG;
% 
% config=[];
% ind=1;
% for k=1:88
%     code=ov_line(k,5);
%     a=find(config==code);
%     config=[config; code];
%     if isempty(a)
%         Z_abc_OH{ind}=Z_OH{1,k};
%         ind=ind+1;
%     end
% end
%     
% config=[];
% ind=1;
% for k=1:117
%     code=ug_line(k,5);
%     a=find(config==code);
%     config=[config; code];
%     if isempty(a)
%         Z_abc_UG{ind}=Z_UG{1,k};
%         ind=ind+1;
%     end
% end
    






