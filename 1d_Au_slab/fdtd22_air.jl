#1D FDTD simulation
#Wave propagation in air toward z-direction
#Mur's absorbing boundary conditions are applied

#Code written in Julia programming language
#Konstantinos Prokopidis, August, 24,2015
#Revised January, 16, 2018

#The code is tested in Julia 1.0.1

using DelimitedFiles

const NZ=8000; #computational domain
const Nt=20000; #time steps

#arrays
global EX=zeros(Float64,NZ+1);
global HY=zeros(Float64,NZ+1);
global INC_ADE=zeros(Float64,1,Nt);

#constants
const m0=12.566e-7; #permability of free space
const e0=8.854e-12; #permittivity of free space
const c0=1/sqrt(e0*m0); #velocity of light

function fdtd22_air(Nt) #Nt the number of time steps

Q=0.3; #Courant number
dz=1e-9; #space step
dt=Q*dz/c0; #time step

#slab position of the total field
za=2000;
zb=2019;

# source constants
zsource=3; #the source location

#source
f_high=3e15;
tw=2.0/(pi*f_high);
t0=4.0*tw;

func=zeros(Float64,Nt+1);
func[1:Nt]=exp.(-((1:Nt).*dt-t0*ones(Nt)).^2.0 ./tw^2.0);

#FDTD coefficients
cH=dt/(m0*dz);
cE=dt/(e0*dz);
c_Mur=(c0*dt-dz)/(c0*dt+dz);

for n=1:Nt

   println(n);

   #Mur's ABC at k=NZ (part a)
   EX[NZ]=EX[NZ-1]-c_Mur*EX[NZ];

   #Mur's ABC at k=1 (part a)
   EX[1]=EX[2]-c_Mur*EX[1];

   EX[2:NZ-1]=EX[2:NZ-1]-cE*(HY[2:NZ-1]-HY[1:NZ-2]);

   #Mur's ABC at k=NZ (part b)
   EX[NZ]=EX[NZ]+c_Mur*EX[NZ-1];

   #Mur's ABC at k=1 (part b)
   EX[1]=EX[1]+c_Mur*EX[2];

   #source
   EX[zsource]=EX[zsource]+func[n]; #soft source

   HY[1:NZ-1]=HY[1:NZ-1]-cH*(EX[2:NZ]-EX[1:NZ-1]);

   #Observation point
   INC_ADE[n]=EX[zb+1];


 end #of for loop

end #of function

@time fdtd22_air(Nt);

#write into file
writedlm("inc.dat",INC_ADE);
