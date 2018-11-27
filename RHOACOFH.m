clc
clear all
close all
format long
tic
SIOX1
%number of adsorbent atoms
nox=length(OX(:,1));
nsi=length(SI(:,1));

%%%number of initial population
minpop=2;
pop=minpop;

%%%critcal point Of H2
TC=33;   %K
PC=12.8; %atm
w=-0.218;
R=82.056e-6;  %atm.m3/gmol.K



%%%Temprature & Pressure &...
P=.005*.986;%bar
T=80;%K
avag=6.023e23;% avagadro number
KA=8.314/avag;%J/gmole.K
KB=1.3806503e-23; %J/K

beta=1/(KB*T);
Rg=8.314;   %J/gmole.K
PMax=1;
deltap=.05;


%parameter of ACO
promote=1.1;
punish=.9;
evap=0.9999;
min=0.1;
max=10;
%parameters of LJ potential For H-SIOX
eps=[82.9848 44.2697 56.6702;44.2697 23.6164 30.2317;56.6702 30.2317 38.7];
sigma=[3.3 3.7229 3.03;3.7229 4.2 3.4182;3.03 3.4182 2.782];%K

eps=eps*KB;   %J/molecule
eps=eps*avag;

%%%%broglie wavelength for 
h=6.62*10^-34;
a=log((sqrt(h^2/(2*pi*2*KB*T))^3)/KB*T);
h=6.62*10^-34;%plank constant
MH2=2;
MSI=28;
MOX=16;
mr=.5;
%Structure of Adsorbenet
mind=zeros(1,3);
for j=1:3
     %maxd(j)=max(OX(:,j));
    maxd(j)=30;

end

%Accuracy of structure of lattice
accrx=.5;
accry=.5;
accrz=.5;
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
maxiter=500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Lattice

Vlattice=1e-30;
for i=1:3
    Vlattice=maxd(i)*Vlattice;%m^3
end
dens=.6e3;%Kg/m^3
Ml=dens*Vlattice*1000;% weight of lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%

ii=1;
 while P<PMax
deltaY=-1;
   while deltaY<0
        
        
%          jflag=0;
%         while jflag==0
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
     J(i)=ceil(rand()*ny); %y position
     K(i)=ceil(rand()*ny); %y position 
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
                   Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6+(h^2/24*mr*Rg*T)*(132/raa(i,k)^2)*(sigma(3,3)/raa(i,k))^12-(30/raa(i,k)^2)*(sigma(3,3)/raa(i,k))^6);
                end
          TEAA=TEAA+Eaa(i,k);
        end

        TE(i)=TEAO+TEAS+TEAA;  
         if TE(i)>0 
                Iflag=1;
         end
              if Iflag==1 & i==pop
                 i=pop-1;
              end
 end
TE;
A;
 
iter=1;
         %****start ANTCo Algorithm
         %Jflag=0;
         STE1=0;
 %while Jflag==0
      while iter<maxiter
                               
          for i=1:pop
       
              
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
        
          Bflag=0;
           
          while Bflag==0
               r=rand();
          
              if r<cpr(1)  
                 path=1;
                 II(i)=I(i);
                 JJ(i)=J(i)-1;
                 KK(i)=K(i);
                 Bflag=1;
              end
          
           if r>cpr(1) & r<cpr(2) 
               path=2;
               II(i)=I(i);
               JJ(i)=J(i);
               KK(i)=K(i)+1;
               Bflag=1;
           end
           
            
           if r>cpr(2) & r<cpr(3) 
               path=3;
               II(i)=I(i)+1;
               JJ(i)=J(i);
               KK(i)=K(i);
               Bflag=1;
           end
            
             
           if r>cpr(3) & r<cpr(4) 
              path=4;
              II(i)=I(i);
              JJ(i)=J(i)+1;
              KK(i)=K(i);
              Bflag=1;
          end
           
          if r>cpr(4) & r<cpr(5) 
               path=5;
               II(i)=I(i);
               JJ(i)=J(i);
               KK(i)=K(i)-1;
               Bflag=1;
          end
             
   
          if r>cpr(5)  
              path=6;
              II(i)=I(i)-1;
              JJ(i)=J(i);
              KK(i)=K(i);
              Bflag=1;
          end
         
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
                    Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6+(h^2/24*mr*Rg*T)*(132/raa(i,k)^2)*(sigma(3,3)/raa(i,k))^12-(30/raa(i,k)^2)*(sigma(3,3)/raa(i,k))^6);
                 end
         TEAA=TEAA+Eaa(i,k);
        end
  
             NTE(i)=TEAO+TEAS+TEAA;
           
        
    delta=NTE(i)-TE(i);
    
  %%%%%%%%%%%%%%update B/respect to delat%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if delta<0
         
                     if path==1  
                                
                                       B(I(i),J(i),K(i),1)=promote*B(I(i),J(i),K(i),1);
                                
                                
                            
                                        B(II(i),JJ(i),KK(i),4)=punish*B(II(i),JJ(i),KK(i),4);
                             
                     end
                     
                     if path==2 
                               
                                  B(I(i),J(i),K(i),2)=promote*B(I(i),J(i),K(i),2);
                              
                             
                                
                                     B(II(i),JJ(i),KK(i),5)=punish*B(II(i),JJ(i),KK(i),5);
                             
                     end
                     
                     if path==3 
                                
                                     B(I(i),J(i),K(i),3)=promote*B(I(i),J(i),K(i),3);
                             
                                 
                          
                                     B(II(i),JJ(i),KK(i),6)=punish*B(II(i),JJ(i),KK(i),6);
                                
                     end
                     
                     if path==4 
                            
                                    B(I(i),J(i),K(i),4)=promote*B(I(i),J(i),K(i),4);
                             
                              
                              
                                    B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                           
                     end
                     if path==5 
                    
                              B(I(i),J(i),K(i),5)=promote* B(I(i),J(i),K(i),5);
                      
                         
                          
                                B(II(i),JJ(i),KK(i),2)=punish*B(II(i),JJ(i),KK(i),2);
                          
                     end
                     if path==6 
                              
                                   B(I(i),J(i),K(i),6)=promote*B(I(i),J(i),K(i),6);
                    
                              
                            
                                      B(II(i),JJ(i),KK(i),3)=punish*B(II(i),JJ(i),KK(i),3);
                             
                      end

           end
    
     if delta>0
                     if path==1  
                                  
                                         B(I(i),J(i),K(i),1)=punish*B(I(i),J(i),K(i),1);
                                 
                                  
                       
                                         B(II(i),JJ(i),KK(i),4)=promote*B(II(i),JJ(i),KK(i),4);
                                
                         end
                     
                     
                     if path==2 
                        
                                   B(I(i),J(i),K(i),2)=punish*B(I(i),J(i),K(i),2);
                              
                              
                              
                                  B(II(i),JJ(i),KK(i),5)=promote*B(II(i),JJ(i),KK(i),5);
                              
                     end
                     
                     if path==3 
                             
                                    B(I(i),J(i),K(i),3)=punish*B(I(i),J(i),K(i),3);
                              
                              
                              
                                  B(II(i),JJ(i),KK(i),6)=promote*B(II(i),JJ(i),KK(i),6);
                              
                     end
                     
                     if path==4 
                              
                                  B(I(i),J(i),K(i),4)=punish*B(I(i),J(i),K(i),4);
                              
                              
                              
                                  B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                              
                     end
                  
                     
                     if path==5
                               
                                    B(I(i),J(i),K(i),5)=punish*B(I(i),J(i),K(i),5);
                               
                               
                               
                                    B(II(i),JJ(i),KK(i),2)=promote*B(II(i),JJ(i),KK(i),2);
                               
                     end
                     
                     if path==6
                              
                                    B(I(i),J(i),K(i),6)=punish*B(I(i),J(i),K(i),6);
                              
                       
                                  B(II(i),JJ(i),KK(i),3)=promote*B(II(i),JJ(i),KK(i),3);
                              
                     end
     end
%       B(II(i),JJ(i),KK(i),:)
                         I(i)=II(i);
                         J(i)=JJ(i);
                         K(i)=KK(i);
                         TE(i)=NTE(i);
                   
          end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
  xx(iter)=iter;
   E(iter)=sum(TE);
     if E(iter)<MIN
             best=A;
             MIN=E(iter);
    end
        
  iter=iter+1;
   B=B*evap;
    end
    MIN
    
% STE2=E(iter-1);
% deltaS=abs(STE2-STE1);
 TotalE(pop)=MIN;
% e=abs(STE2*.01);
%    if deltaS<e
%       Jflag=1;
%    else
%       STE1=STE2;
%       iter=1;
%    end

       
    
  % end   
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%peng-robinson equation%%%%%%%%%%
  PP=P/(1.01325);%bar to atm
  apr=(0.457235*(R*TC)^2)/PC;
  bpr=0.077796*R*TC/PC;
  kpr=0.37464+1.54226*w-0.26992*w^2;
  alfapr=(1+kpr*(1-sqrt(T/TC)))^2;

  Apr=apr*alfapr*PP/((R*T)^2);
  Bpr=bpr*(PP/(R*T));

  a1=-(1-Bpr);
  a2=Apr-2*Bpr-3*Bpr^2;
  a3a=-Apr*Bpr+Bpr^2+Bpr^3;

  zzz=solve('zz^3+a1*zz^2+a2*zz+a3a=0');
  zzzz=subs(zzz)
  zg=0;
  for dd=1:3
      zg1(dd)=real(zzzz(dd));
      if zg1(dd)>zg
  zg=zg1(dd);
      end
  end
zg;
  lnfc=(zg-1)-log(zg-Bpr)-((Apr)/(2*sqrt(2)*Bpr))*log((zg+(1+sqrt(2))*Bpr)/(zg+(1-sqrt(2))*Bpr))
  miux=KB*T*lnfc;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  miu=miux+KB*T*log(pop);
  Y(pop)=-pop*miu+KB*T*log(factorial(pop))+TotalE(pop)
     if pop>minpop
        deltaY=Y(pop)-Y(pop-1);
     end
  pop=pop+1
   end
  
  mol(ii)=pop-1;
  MH2=(pop-1)*2/avag;
  Mlattice=MH2+dens*Vlattice*1000;
  teta(ii)=(MH2/Mlattice)*100
  Pressure(ii)=P
  P=P*2;
  ii=ii+1;
  minpop=pop-1;
  pop=minpop;
 end
 
 
 mmol=(mol/avag)*1000;
TETA=mmol/Ml;
plot(Pressure,TETA)
xlabel('pressur(atm)');
ylabel('teta(mmol/g)');
TETA'
toc