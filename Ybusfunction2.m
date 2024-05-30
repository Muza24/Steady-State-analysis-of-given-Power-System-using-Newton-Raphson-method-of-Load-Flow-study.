function [Z,LD,Y_Bus,nbus,Y_b] = Ybusfunction2()

    %  From To  Line impedance  jB/2
LD=[1   2  0.02+1i*0.04     0.0j
    2   3  0.0125+1i*0.025  0.0j
    3   1  0.01+1i*0.03     0.0j];


nbus=max(max(LD(:,1)),max(LD(:,2)));
Z=zeros(nbus);
LD_size=size(LD,1);

for i=1:LD_size
    if Z(LD(i,1),LD(i,2))==0
        Z(LD(i,1),LD(i,2))=LD(i,3);
        Z(LD(i,2),LD(i,1))=LD(i,3);
    else
        %For accomodating the parallel lines between buses
        Z(LD(i,1),LD(i,2))=(Z(LD(i,1),LD(i,2))*LD(i,3))/(Z(LD(i,1),LD(i,2))+LD(i,3));
        Z(LD(i,2),LD(i,1))=(Z(LD(i,2),LD(i,1))*LD(i,3))/(Z(LD(i,2),LD(i,1))+LD(i,3));
    end
end
Y=zeros(nbus);
for i=1:nbus
   for j=1:nbus
      if Z(i,j)~=0
          Y(i,j)=1/Z(i,j);
      end
   end
end
%Updating the shunt admittance values
for i=1:LD_size
   k=LD(i,1);
   j=LD(i,2);
   Y(k,k)=Y(k,k)+LD(i,4);
   Y(j,j)=Y(j,j)+LD(i,4);
end
%Generating the shunt admittance matrix
Y_b=zeros(nbus);
for i=1:LD_size
   Y_b(LD(i,1),LD(i,2))=(LD(i,4));
   Y_b(LD(i,2),LD(i,1))=(LD(i,4));
end

%Forming the Y-bus
Y_Bus=zeros(nbus);
for i=1:nbus
   
   Y_Bus(i,i)=sum(Y(i,:));
   for j=1:nbus
        if i~=j
            Y_Bus(i,j)=-Y(i,j);%for off diagonal terms
        end
   end
end
%Y_Bus
end