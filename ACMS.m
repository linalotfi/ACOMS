clc
clear all
close all
format long


CHA
%number of adsorbent atoms
nox=length(OX(:,1));
nsi=length(SI(:,1));

%%%number of initial population
minpop=4;
pop=minpop;

%%%critcal point Of CO2
TC=31.2+273;   %K
PC=73.8*.986; %atm
w=0.228;
R=82.056e-6;  %atm.m3/gmol.K


%%%Temprature & Pressure &...
P=.1*.986;%bar
T=303;%K
K=8.314;%J/gmole.K
KB=1.3806503e-23; %J/K
avag=6.023e23;% avagadro number
beta=1/(KB*T);
Rg=8.314;   %J/gmole.K
PMax=1.1*.986;
deltap=.1;


%parameter of ACO
promote=2;
punish=0.0;
evap=0.95;
min=0;
max=100;
%parameters of LJ potential For CO2-SIOX
eps=[82.9848 44.2697 124.9;44.2697 23.6164 66.62;124.9 66.62 188];
sigma=[3.3 3.75 3.09;3.75 4.2 3.55;3.09 3.55 2.9];
eps=eps*KB;   %J/molecule
eps=eps*avag;


%Structure of Adsorbenet
mind=zeros(1,3);
for j=1:3
     %maxd(j)=max(OX(:,j));
     maxd(j)=20;

end

%Accuracy of structure of lattice
accrx=0.1;
accry=0.1;
accrz=0.1;
nx=round((maxd(1,1)-mind(1,1))/accrx);
ny=round((maxd(1,2)-mind(1,2))/accry);
nz=round((maxd(1,3)-mind(1,3))/accrz);

B=ones(nx,ny,nz,6); %pheromone level of paths
B(:,1,:,1)=0;
B(:,:,nz,2)=0;
B(nx,:,:,3)=0;
B(:,ny,:,4)=0;
B(:,:,1,5)=0;
B(1,:,:,6)=0;
maxiter=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Lattice

Vlattice=1e-30;
for i=1:3
    Vlattice=maxd(i)*Vlattice;%m^3
end
dens=1.2e3;%Kg/m^3
Ml=dens*Vlattice*1000;% weight of lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%


% while P<PMax
%     deltaY=-1;
%     while deltaY<0
        
        
        
MIN=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generating initial population%%%%%%%%%




  i=0;
           Iflag=0;
          while i<pop
                  if Iflag==0
                    i=i+1;
                  end
         Iflag=0;
         TEAS=0;
         TEAO=0;
         TEAA=0;
         CTEAA=0;
     %generating  random ANTs
     I(i)=ceil(rand()*nx); %x position 
     J(i)=ceil(rand()*ny);%y position
     K(i)=ceil(rand()*ny);%y position 
     x(i)=mind(1,1)+(maxd(1,1)-mind(1,1))*(I(i)/nx);
     y(i)=mind(1,2)+(maxd(1,2)-mind(1,2))*(J(i)/ny);
     z(i)=mind(1,3)+(maxd(1,3)-mind(1,3))*(K(i)/nz);
     A(i,1)=x(i);
     A(i,2)=y(i);
     A(i,3)=z(i);
    
     %calculating distanse between ONE atome of adsorbate and other
     %adsorbate and adsorbent atoms &&& ENERGY of ONE atom of adsorbate
        for k=1:nox
         raox(i,k)=norm(OX(k,:)-A(i,:));
         Eaox(i,k)=4*eps(1,3)*((sigma(1,3)/raox(i,k))^12-(sigma(1,3)/raox(i,k))^6);
         TEAO=TEAO+Eaox(i,k);
        end
      
        for k=1:nsi
          rasi(i,k)=norm(SI(k,:)-A(i,:));
          Easi(i,k)=4*eps(2,3)*((sigma(2,3)/rasi(i,k))^12-(sigma(2,3)/rasi(i,k))^6);
          TEAS=TEAS+Easi(i,k);
        end
      
        for k=1:i
          raa(i,k)=norm(A(k,:)-A(i,:));
          raa(k,i)=raa(i,k);
            if i==k
               Eaa(i,k)=0;
            end
            if i~=k
            Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6);
            end
            TEAA=TEAA+Eaa(i,k);
        end
%         if i>1
%         for k=1:i-1
%             raa(k,i)=norm(A(i,:)-A(k,:));
%             Eaa(k,i)=4*eps(3,3)*((sigma(3,3)/raa(k,i))^12-(sigma(3,3)/raa(k,i))^6);
%             CTEAA=CTEAA+Eaa(k,i);
%             TE(i-1)=TE(i-1)+CTEAA;
%         end
%         end
        TE(i)=TEAO+TEAS+TEAA;  
         if TE(i)>0 
              Iflag=1;
           end
           if Iflag==1 & i==pop
              i=pop-1;
           end
 end
TE
% break
iter=1;
         %****start ANTCo Algorithm
while iter<maxiter
                           
      for i=1:pop
              path=0;
              
              Sigma=0;
           for ij=1:6
                 Sigma=Sigma+(B(I(i),J(i),K(i),ij));
           end
           
           for j=1:6
                   ipr(j)=B(I(i),J(i),K(i),j)/Sigma;  
           end
           
           cpr(1)=ipr(1);
           for j=2:6
               cpr(j)=ipr(j)+cpr(j-1);
           end

           r=rand();
           if r<cpr(1)
               path=1;
               II(i)=I(i);
               JJ(i)=J(i)-1;
               KK(i)=K(i);
           end
           if r>cpr(1) & r<cpr(2) 
            path=2;
           II(i)=I(i);
           JJ(i)=J(i);
           KK(i)=K(i)+1;
           end
           
           if r>cpr(2) & r<cpr(3)  
             path=3;
             II(i)=I(i)+1;
             JJ(i)=J(i);
             KK(i)=K(i);
           end
     
           if r>cpr(3)   & r<cpr(4)
              path=4;
              II(i)=I(i);
              JJ(i)=J(i)+1;
              KK(i)=K(i);
          end
          
          if r>cpr(4)   & r<cpr(5)
              path=5;
              II(i)=I(i);
              JJ(i)=J(i);
              KK(i)=K(i)-1;
          end
  
          if r>cpr(5)    & I(i)>1
              path=6;
              II(i)=I(i)-1;
              JJ(i)=J(i);
              KK(i)=K(i);
          end
                                          
          x(i)=mind(1,1)+(maxd(1,1)-mind(1,1))*(II(i)/nx);
          y(i)=mind(1,2)+(maxd(1,2)-mind(1,2))*(JJ(i)/ny);
          z(i)=mind(1,3)+(maxd(1,3)-mind(1,3))*(KK(i)/nz);
          A(i,1)=x(i);
          A(i,2)=y(i);
          A(i,3)=z(i);
    
         TEAS=0;
         TEAO=0;
         TEAA=0;
         CTEAA=0;
    
    for k=1:nox
         raox(i,k)=norm(OX(k,:)-A(i,:));
         Eaox(i,k)=4*eps(1,3)*((sigma(1,3)/raox(i,k))^12-(sigma(1,3)/raox(i,k))^6);
         TEAO=TEAO+Eaox(i,k);
        end
      
        for k=1:nsi
          rasi(i,k)=norm(SI(k,:)-A(i,:));
          Easi(i,k)=4*eps(2,3)*((sigma(2,3)/rasi(i,k))^12-(sigma(2,3)/rasi(i,k))^6);
          TEAS=TEAS+Easi(i,k);
        end
      
        for k=1:i
          raa(i,k)=norm(A(k,:)-A(i,:));
          raa(k,i)=raa(i,k);
            if i==k
               Eaa(i,k)=0;
            end
            if i~=k
            Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6);
            end
            TEAA=TEAA+Eaa(i,k);
        end
        if i>1
%         for k=1:i-1
%             raa(k,i)=norm(A(i,:)-A(k,:));
%             Eaa(k,i)=4*eps(3,3)*((sigma(3,3)/raa(k,i))^12-(sigma(3,3)/raa(k,i))^6);
%             CTEAA=CTEAA+Eaa(k,i);
%             NTE(i-1)=NTE(i-1)+CTEAA;
%         end
        end
             NTE(i)=TEAO+TEAS+TEAA;
                                end
    
    
    
    
    for i=1:pop
        
    delta=NTE(i)-TE(i);
   
           if delta<0
         
                     if path==1  
                         if B(I(i),J(i),K(i),1)<max
                           B(I(i),J(i),K(i),1)=promote* B(I(i),J(i),K(i),1);
                         end
                           B(II(i),JJ(i),KK(i),4)=punish*B(II(i),JJ(i),KK(i),4);
                     end
                     
                     if path==2 
                         if  B(I(i),J(i),K(i),2)<max
                        B(I(i),J(i),K(i),2)=promote* B(I(i),J(i),K(i),2);
                         end
                        B(II(i),JJ(i),KK(i),4)=punish*B(II(i),JJ(i),KK(i),5);
                 
                     end
                     
                     if path==3 
                         if  B(I(i),J(i),K(i),3) < max
                     B(I(i),J(i),K(i),3)=promote* B(I(i),J(i),K(i),3);
                         end
                     B(II(i),JJ(i),KK(i),6)=punish*B(II(i),JJ(i),KK(i),6);
                     end
                     
                     if path==4 
                         if  B(I(i),J(i),K(i),4)<max
                         B(I(i),J(i),K(i),4)=promote* B(I(i),J(i),K(i),4);
                         end
                         B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                     end
                     if path==5 
                         if B(I(i),J(i),K(i),5)<max
                         B(I(i),J(i),K(i),5)=promote* B(I(i),J(i),K(i),5);
                         end
                         B(II(i),JJ(i),KK(i),2)=punish*B(II(i),JJ(i),KK(i),2);
                     end
                     if path==6 
                         if  B(I(i),J(i),K(i),6)<max
                         B(I(i),J(i),K(i),6)=promote* B(I(i),J(i),K(i),6);
                         end
                         B(II(i),JJ(i),KK(i),3)=punish*B(II(i),JJ(i),KK(i),3);
                     end
                     
%                      I(i)=II(i);
%                      J(i)=JJ(i);
%                      K(i)=KK(i);
%                      TE(i)=NTE(i);
          end
     if delta>0
                     if path==1  
                           B(I(i),J(i),K(i),1)=punish* B(I(i),J(i),K(i),1);
                                    if  B(II(i),JJ(i),KK(i),4)<max
                                          B(II(i),JJ(i),KK(i),4)=promote*B(II(i),JJ(i),KK(i),4);
                                     end
                       end
                     
                     
                     if path==2 
                        
                        B(I(i),J(i),K(i),2)=punish* B(I(i),J(i),K(i),2);
                                      if  B(II(i),JJ(i),KK(i),4)<max
                                               B(II(i),JJ(i),KK(i),4)=promote*B(II(i),JJ(i),KK(i),5);
                                        end
                     end
                     
                     if path==3 
                     B(I(i),J(i),K(i),3)=punish* B(I(i),J(i),K(i),3);
                               if B(II(i),JJ(i),K(i),6)<max
                                  B(II(i),JJ(i),K(i),6)=promote*B(II(i),JJ(i),KK(i),6);
                                   end
                     end
                     
                     if path==4 
                         B(I(i),J(i),K(i),4)=punish* B(I(i),J(i),K(i),4);
                             if   B(II(i),JJ(i),KK(i),1)<max
                                   B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                             end
                     end
                     
                     if path==5
                        B(I(i),J(i),K(i),5)=punish* B(I(i),J(i),K(i),5);
                                 if    B(II(i),JJ(i),K(i),2)<max
                                 B(II(i),JJ(i),KK(i),2)=promote*B(II(i),JJ(i),KK(i),2);
                                  end
                     end
                     
                     if path==6
                         B(I(i),J(i),K(i),6)=punish* B(I(i),J(i),K(i),6);
                         if  B(II(i),JJ(i),K(i),3)<max
                             B(II(i),JJ(i),K(i),3)=promote*B(II(i),JJ(i),K(i),3);
                         end
                     end
     end
                         I(i)=II(i);
                         J(i)=JJ(i);
                         K(i)=KK(i);
                         TE(i)=NTE(i);

    end
  xx(iter)=iter;
  E(iter)=sum(TE);
     if E(iter)<MIN
             best=A;
             MIN=E(iter);
    end
             
  iter=iter+1;
  B=B*evap;
  
                       end
                       plot(xx,E)