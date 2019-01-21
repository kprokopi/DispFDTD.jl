#1D FDTD simulation
#A wave impinges upon a gold layer (20 nm thick)
#The gold is modelled by the Drude-CP model
#and the FDTD method is an ADE algorithm

#Code written in Julia programming language
#Konstantinos Prokopidis, August, 24,2015
#Revised January, 16, 2018

#The code is tested in Julia 1.0.1

using DelimitedFiles

const NZ=8000; #computational domain
const Nt=20000; #time steps

#matrices
global EX=zeros(Float64,3,NZ);
global HY=zeros(Float64,1,NZ);
global TCP_ADE=zeros(Float64,1,Nt);
global PCP_ADE=zeros(Float64,1,Nt);

const m0=12.566e-7; #permeability of free space
const e0=8.854187817620e-12; #permittivity of free space
const c0=1/sqrt(e0*m0); #velocity of light

function fdtd22_DCP_ADE(Nt) #Nt is th number of time steps

#-----------------------------------------
Q=0.3; #Courant number
dz=1e-9; #space step
dt=Q*dz/c0; #time step

#-----------------------------------------

#slab of 20nm
za=2000;
zb=2019;

# Parameters of Au from paper of Vial 2011
 epsilon_inf=1.1431;
 omega_D=1.3202e16;
 gamma=1.0805e14;

 A=[0.26698 3.0834];
 phi_m=[-1.2371 -1.0968];
 Omega=[3.8711e15 4.1684e15];
 Gamma=[4.4642e14 2.3555e15];


#matrices
#EX and HY are defined as global variables
DX=zeros(Float64,1,NZ);

PX_CP1=zeros(Float64,3,NZ);
PX_CP2=zeros(Float64,3,NZ);
PX_CP3=zeros(Float64,3,NZ);

# source constants
zsource=3; #the source location

#source
f_high=3e15;
tw=2/(pi*f_high);
t0=4*tw;

func=zeros(Float64,Nt);
func[1:Nt]=exp.(-((1:Nt).*dt-t0*ones(Nt)).^2.0 ./tw^2.0);

#FDTD coefficients
#Coefficients for the Drude-CP model
#FDTD coefficients for the critical point model
a_0=e0*2*A.*Omega.*(Omega.*cos.(phi_m)-Gamma.*sin.(phi_m));
a_1=-e0*2*A.*Omega.*sin.(phi_m);
b_0=Omega.^2+Gamma.^2;
b_1=2*Gamma;
b_2=[1 1];


#We extend the arrays of the coefficients in order to include the Drude
#parameters
a_0=[a_0 e0*omega_D^2];
a_1=[a_1 0];
b_0=[b_0 0];
b_1=[b_1 gamma];
b_2=[b_2 1];

#FDTD coefficients

C=b_2./(dt^2) + b_1 ./(2*dt) + b_0 ./4;

C1=(2*b_2/(dt^2) - b_0/2)./C;
C2=(b_1./(2*dt)-b_2./(dt^2)-b_0./4)./C;
C3=(a_0./4+a_1./(2*dt))./C;
C4=a_0./(2*C);
C5=(a_0./4-a_1/(2*dt))./C;

c_1=e0*epsilon_inf+(C3[1]+C3[2]+C3[3]);
c_2=C4[1]+C4[2]+C4[3];
c_3=C5[1]+C5[2]+C5[3];

cH=dt/(m0*dz);
cE=dt/(e0*dz);

c_Mur=(c0*dt-dz)/(c0*dt+dz);

for n=1:Nt
   println(n);
   #Mur's ABC at k=NZ (part a)
   EX[1,NZ]=EX[1,NZ-1]-c_Mur*EX[1,NZ];

   #Mur's ABC at k=1 (part a)
   EX[1,1]=EX[1,2]-c_Mur*EX[1,1];

   #free space
   EX[1,2:za-1]=EX[1,2:za-1]-cE*(HY[1,2:za-1]-HY[1,1:za-2]);

  # slab of gold [za,zb]

  DX[1,za:zb]=DX[1,za:zb]-(dt/dz)*(HY[1,za:zb]-HY[1,za-1:zb-1]);

  #back-storing
  EX[3,za:zb]=EX[2,za:zb];
  EX[2,za:zb]=EX[1,za:zb];

  #update EX

  EX[1,za:zb]=(1/c_1)*(DX[1,za:zb]
      -C1[1]*PX_CP1[1,za:zb]
      -C1[2]*PX_CP2[1,za:zb]
      -C1[3]*PX_CP3[1,za:zb]
      -C2[1]*PX_CP1[2,za:zb]
      -C2[2]*PX_CP2[2,za:zb]
      -C2[3]*PX_CP3[2,za:zb]
      -c_2*EX[2,za:zb]-c_3*EX[3,za:zb]);


  #back-storing
  PX_CP1[3,za:zb]=PX_CP1[2,za:zb];
  PX_CP1[2,za:zb]=PX_CP1[1,za:zb];

  PX_CP2[3,za:zb]=PX_CP2[2,za:zb];
  PX_CP2[2,za:zb]=PX_CP2[1,za:zb];

  PX_CP3[3,za:zb]=PX_CP3[2,za:zb];
  PX_CP3[2,za:zb]=PX_CP3[1,za:zb];

  #update PX_CP
  PX_CP1[1,za:zb]=C1[1]*PX_CP1[2,za:zb]+C2[1]*PX_CP1[3,za:zb]+C3[1]*EX[1,za:zb]+C4[1]*EX[2,za:zb]+C5[1]*EX[3,za:zb];
  PX_CP2[1,za:zb]=C1[2]*PX_CP2[2,za:zb]+C2[2]*PX_CP2[3,za:zb]+C3[2]*EX[1,za:zb]+C4[2]*EX[2,za:zb]+C5[2]*EX[3,za:zb];
  PX_CP3[1,za:zb]=C1[3]*PX_CP3[2,za:zb]+C2[3]*PX_CP3[3,za:zb]+C3[3]*EX[1,za:zb]+C4[3]*EX[2,za:zb]+C5[3]*EX[3,za:zb];

   #free space

   EX[1,zb+1:NZ-1]=EX[1,zb+1:NZ-1]-cE*(HY[1,zb+1:NZ-1]-HY[1,zb:NZ-2]);

   #Mur's ABC at k=NZ (part b)
   EX[1,NZ]=EX[1,NZ]+c_Mur*EX[1,NZ-1];

   #Mur's ABC at k=1 (part b)
   EX[1,1]=EX[1,1]+c_Mur*EX[1,2];

   #source (soft)
   EX[1,zsource]=EX[1,zsource]+func[n];

   #Observation point

   PCP_ADE[n]=EX[1,za-1];
   TCP_ADE[n]=EX[1,zb+1];

   #H -component
   HY[1,1:NZ-1]=HY[1,1:NZ-1]-cH*(EX[1,2:NZ]-EX[1,1:NZ-1]);


end


end #of function
@time fdtd22_DCP_ADE(Nt);

#write into file
writedlm("total.dat",TCP_ADE);
