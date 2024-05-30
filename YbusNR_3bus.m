clc
clear
format compact
[Z,LD,Y_Bus,nbus,Y_b] = Ybusfunction2();
BaseMVA=100;
Nbus = 3;

%1-Slack bus,2-PQ bus,3-PV bus
%          1     2         3        4    5    6       7     8 
%      bus_no bus_type   vmag     delta  Pi   Qi    Qmax  Qmin 
B=[      1       1       1.05       0    0    0      10   -10;
         2       2       1.00       0   -4   -2.5    10   -10;
         3       3       1.04       0    2    0      1   -1];
V=B(:,3).*exp(B(:,4));
V_mag=abs(V);
del=angle(V);
Y_mag = abs(Y_Bus);
theta=angle(Y_Bus);
error =  1;
tolerance = 1e-10;
iter = 0;
while error>tolerance
    iter = iter + 1;
    
    sl = find(B(:,2)==1);
    pq = find(B(:,2)==2);
    pv = find(B(:,2)==3);
    npq = length(pq);
    npv = length(pv);

    V = B(:,3).* exp(1i*B(:,4));
    V_mag = abs(V);
    del = angle(V);
    
    Pi = zeros(Nbus,Nbus);
    Qi = zeros(Nbus,Nbus);

    for i = 1:Nbus
        for k = 1:Nbus
            Pi(i,k)= Pi(i,k)+Y_mag(i,k)*V_mag(i)*V_mag(k)*cos(theta(i,k)-del(i)+del(k));
            Qi(i,k)= Qi(i,k)-Y_mag(i,k)*V_mag(i)*V_mag(k)*sin(theta(i,k)-del(i)+del(k));
        end
    end
    
    Pi = sum(Pi,2);
    Qi = sum(Qi,2);
    del_P = B(:,5)-Pi;
    del_P(sl,:) = [];
    del_Q = B(:,6) - Qi;
    del_Q = del_Q(pq(1:npq));
    
    delpq = [del_P; del_Q];

   % Jacobian element J1
J1 = zeros(Nbus-1,Nbus-1);     
for j = 1:(Nbus-1)
    z=j+1;
    for m = 1:(Nbus-1)
        i=m+1;
        if j==m
            j1=0;
            for k=1:Nbus
                
                j1 = j1 + Y_mag(i,k)*V_mag(i)*V_mag(k)*sin(theta(i,k)-del(i)+del(k));
                
            end
            
            J1(j,m) = j1 - Y_mag(i,i)*(V_mag(i)^2)*sin(theta(i,i));
        else 
            J1(j,m) = -Y_mag(z,i)*V_mag(z)*V_mag(i)*sin(theta(z,i)-del(z)+del(i));
        end
    end
end

% Jacobian element J2
J2 = zeros(Nbus-1,npq);         %npq is the number of pq buses
for j = 1:(Nbus-1)
    z = j+1;
    
    for m=1:npq
        i=m+1;
        
        if (j==m)
            j2 = 0;
            for k=1:Nbus
%                 j2 = j2 + Y_mag(i,k)*V_mag(i)*V_mag(k)*cos(theta(i,k)-del(i)+del(k));
                j2 = j2 + Y_mag(i,k)*V_mag(k)*cos(theta(i,k)-del(i)+del(k));
            end
            
%             J2(j,m)=j2+Y_mag(i,i)*(V_mag(i)^2)*cos(theta(i,i));
            J2(j,m)=j2+Y_mag(i,i)*cos(theta(i,i));
            
        else
%             J2(j,m)=Y_mag(z,i)*V_mag(z)*cos(theta(z,i)-del(z)+del(i));
            J2(j,m)=Y_mag(z,i)*V_mag(z)*cos(theta(z,i)-del(z)+del(i));
        end
    end
end

% Jacobian element J3            
J3 = zeros(npq,Nbus-1)
for j = 1:npq
    z=j+1;
    for m = 1:(Nbus-1)
        i=m+1;
        if j==m
            j3=0;
            for k=1:Nbus
                
                j3 = j3 + Y_mag(i,k)*V_mag(i)*V_mag(k)*cos(theta(i,k)-del(i)+del(k));
                
            end
            
            J3(j,m) = j3 - Y_mag(i,z)*V_mag(i)*V_mag(z)*cos(theta(i,z)-del(i)+del(z));
        else 
            J3(j,m) = -Y_mag(z,i)*V_mag(z)*V_mag(i)*cos(theta(i,z)-del(i)+del(z));
        end
    end
end

% Jacobian element J4
J4 = zeros(npq,npq);
for j=1:npq
    z=j+1;
    for m=1:npq
        i=m+1;
        if j==m
            j4=0;
            for k=1:Nbus
                j4 = j4 -Y_mag(i,k)*V_mag(k)*sin(theta(i,k)-del(i)+del(k));
            end
            J4 = j4 - Y_mag(i,i)*V_mag(i)*sin(theta(i,i)-del(i)+del(i));
        end
    end
end

%Jacobian matrix
J = [J1 J2 ; J3 J4];
delx = J\delpq;         % the backslash \ serves as inverse. i.e. J inverse * deltapq % delpq is the power mismatch

p=0;
for i=1:(Nbus-1)
        j=i+1;
        B(j,4)=B(j,4)+delx(i);
end
V;

for o=1:npq
    i=pq(o);
    B(i,3)=B(i,3)+delx(j);
    j=j+1;
end

error = abs(max(delpq))
end

abs(V)
angle(V)


















           
                
              