clc
clear all
close all
format long
tic

NCHA
%number of adsorbent atoms
nox=length(OX(:,1));
nsi=length(SI(:,1));

%%%number of initial population
minpop=2;
pop=minpop;
avag=6.023e23;% avagadro number
%%%critcal point Of CO2
TC=-82.59+273;   %K
PC=45.99*.986; %atm
w=0.288;
R=82.056e-6;  %atm.m3/gmol.K



%%%Temprature & Pressure &...
P=.1*.986;%bar
T=300;%K
K=8.314;%J/gmole.K
KB=1.3806503e-23; %J/K
beta=1/(KB*T);
Rg=8.314;   %J/gmole.K
PMax=1.1*.986;
deltap=.3;


%parameter of ACO
promote=2;
punish=0.5;
evap=0.999;
min=0.5;
pmax=1000; %max phermone level
%max=10;
%parameters of LJ potential For CO2-SIOX
eps=[82.9848 44.2697 112.679;44.2697 23.6164 60.1027;112.676 60.1027 153];

sigma=[3.3 3.75 3.51;3.75 4.2 3.96;3.51 3.96 3.72];%K

eps=eps*KB;
eps=eps*avag;

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

B=(ones(nx,ny,nz,6))*1000; %pheromone level of paths
B(:,1,:,1)=0;
B(:,:,nz,2)=0;
B(nx,:,:,3)=0;
B(:,ny,:,4)=0;
B(:,:,1,5)=0;
B(1,:,:,6)=0;
maxiter=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Lattice

Vlattice=1e-30;
for i=1:3
    Vlattice=maxd(i)*Vlattice;%m^3
end
dens=1.382e3;%Kg/m^3
Ml=dens*Vlattice*1000;% weight of lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generating initial population%%%%%%%%%


ii=1;
 while P<PMax
deltaY=-1;
   while deltaY<0
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
                   Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6);
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
MIN=sum(TE)
best=A;

iter=1;
         %****start ANTCo Algorithm
        
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
                    Eaa(i,k)=4*eps(3,3)*((sigma(3,3)/raa(i,k))^12-(sigma(3,3)/raa(i,k))^6);
                 end
         TEAA=TEAA+Eaa(i,k);
        end
  
             NTE(i)=TEAO+TEAS+TEAA;
           
        
    delta=NTE(i)-TE(i);
    
  %%%%%%%%%%%%%%update B/respect to delat%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if delta<0
         
                     if path==1  
                                if B(I(i),J(i),K(i),1)<pmax
                                       B(I(i),J(i),K(i),1)=promote*B(I(i),J(i),K(i),1);
                                end
                                
                                if B(II(i),JJ(i),KK(i),4)>min
                                        B(II(i),JJ(i),KK(i),4)=punish*B(II(i),JJ(i),KK(i),4);
                                end
                     end
                     
                     if path==2 
                               if B(I(i),J(i),K(i),2)<pmax
                                  B(I(i),J(i),K(i),2)=promote*B(I(i),J(i),K(i),2);
                               end
                             
                                if B(II(i),JJ(i),KK(i),5)>min
                                     B(II(i),JJ(i),KK(i),5)=punish*B(II(i),JJ(i),KK(i),5);
                                end
                     end
                     
                     if path==3 
                               if B(I(i),J(i),K(i),3)<pmax
                                     B(I(i),J(i),K(i),3)=promote*B(I(i),J(i),K(i),3);
                               end
                                 
                                 if B(II(i),JJ(i),KK(i),6)>min
                                     B(II(i),JJ(i),KK(i),6)=punish*B(II(i),JJ(i),KK(i),6);
                                 end
                     end
                     
                     if path==4 
                             if B(I(i),J(i),K(i),4)<pmax
                                    B(I(i),J(i),K(i),4)=promote*B(I(i),J(i),K(i),4);
                             end
                              
                                if B(II(i),JJ(i),KK(i),1)>min
                                    B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                               end
                     end
                     if path==5 
                       if B(I(i),J(i),K(i),5)<pmax 
                              B(I(i),J(i),K(i),5)=promote* B(I(i),J(i),K(i),5);
                       end
                         
                            if B(II(i),JJ(i),KK(i),2)>min
                                B(II(i),JJ(i),KK(i),2)=punish*B(II(i),JJ(i),KK(i),2);
                            end
                     end
                     if path==6 
                           if B(I(i),J(i),K(i),6)<pmax    
                                   B(I(i),J(i),K(i),6)=promote*B(I(i),J(i),K(i),6);
                           end
                              
                              if B(II(i),JJ(i),KK(i),3)>min
                                      B(II(i),JJ(i),KK(i),3)=punish*B(II(i),JJ(i),KK(i),3);
                              end
                      end

           end
    
     if delta>0
                     if path==1  
                                  if B(I(i),J(i),K(i),1)>min
                                         B(I(i),J(i),K(i),1)=punish*B(I(i),J(i),K(i),1);
                                  end
                                  
                                 if B(II(i),JJ(i),KK(i),4)<pmax
                                         B(II(i),JJ(i),KK(i),4)=promote*B(II(i),JJ(i),KK(i),4);
                                 end
                         end
                     
                     
                     if path==2 
                              if  B(I(i),J(i),K(i),2)>min
                                   B(I(i),J(i),K(i),2)=punish*B(I(i),J(i),K(i),2);
                              end
                              
                          if B(II(i),JJ(i),KK(i),5)<pmax
                                  B(II(i),JJ(i),KK(i),5)=promote*B(II(i),JJ(i),KK(i),5);
                          end
                     end
                     
                     if path==3 
                              if  B(I(i),J(i),K(i),3)>min
                                    B(I(i),J(i),K(i),3)=punish*B(I(i),J(i),K(i),3);
                              end
                              
                             if B(II(i),JJ(i),KK(i),6)<pmax
                                  B(II(i),JJ(i),KK(i),6)=promote*B(II(i),JJ(i),KK(i),6);
                             end
                     end
                     
                     if path==4 
                              if B(I(i),J(i),K(i),4)>min
                                  B(I(i),J(i),K(i),4)=punish*B(I(i),J(i),K(i),4);
                              end
                              
                              if B(II(i),JJ(i),KK(i),1)<pmax
                                  B(II(i),JJ(i),KK(i),1)=punish*B(II(i),JJ(i),KK(i),1);
                              end
                     end
                  
                     
                     if path==5
                               if B(I(i),J(i),K(i),5)>min
                                    B(I(i),J(i),K(i),5)=punish*B(I(i),J(i),K(i),5);
                               end
                               
                             if B(II(i),JJ(i),KK(i),2)<pmax
                                    B(II(i),JJ(i),KK(i),2)=promote*B(II(i),JJ(i),KK(i),2);
                             end
                     end
                     
                     if path==6
                              if  B(I(i),J(i),K(i),6)>min
                                    B(I(i),J(i),K(i),6)=punish*B(I(i),J(i),K(i),6);
                              end
                              
                             if B(II(i),JJ(i),KK(i),3)<pmax
                                  B(II(i),JJ(i),KK(i),3)=promote*B(II(i),JJ(i),KK(i),3);
                             end
                     end
     end
%      B(I(i),J(i),K(i),:)
%      i
%      I 
%      J
%      K
%      pause

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
 TotalE(pop)=MIN;


     
        
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
  zzzz=subs(zzz);
  zg=0;
  for dd=1:3
      zg1(dd)=real(zzzz(dd));
      if zg1(dd)>zg
  zg=zg1(dd);
      end
  end

  lnfc=(zg-1)-log(zg-Bpr)-((Apr)/(2*sqrt(2)*Bpr))*log((zg+(1+sqrt(2))*Bpr)/(zg+(1-sqrt(2))*Bpr));
  miux=KB*T*lnfc;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  miu=miux+KB*T*log(pop);
  Y(pop)=(-pop*miu+KB*T*log(factorial(pop)))*avag+TotalE(pop);
     if pop>minpop
        deltaY=Y(pop)-Y(pop-1);
     end
  pop=pop+1
 B=(ones(nx,ny,nz,6))*1000; %pheromone level of paths
 B(:,1,:,1)=0;
 B(:,:,nz,2)=0;
 B(nx,:,:,3)=0;
 B(:,ny,:,4)=0;
 B(:,:,1,5)=0;
 B(1,:,:,6)=0;

   end
  
  mol(ii)=pop-1;
  MH2=(pop-1)*2/avag;
  Mlattice=MH2+dens*Vlattice*1000;
  teta(ii)=(MH2/Mlattice)*100
  Pressure(ii)=P
  P=P+deltap;
  ii=ii+1;
  minpop=pop-1;
  pop=minpop;
 end
 
 

TETA=(mol/18304)*1000;
plot(Pressure,TETA)
xlabel('pressur(atm)');
ylabel('teta(mmol/g)');
TETA'
toc
  % end   
