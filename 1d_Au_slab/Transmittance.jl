# Calculation of the transmittance
# using the proposed FDTD method vs analytical

#Code written in Julia programming language
#Konstantinos Prokopidis, August, 24,2015
#Revised in January, 21, 2019

#The code is tested in Julia 1.0.1

#To use fft function first install FFTW package:
#using Pkg
#Pkg.add("FFTW")

using DelimitedFiles
using FFTW
using PyPlot

#solution of the total field
tot_field=readdlm("total.dat");
#solution of the incident field
inc_field=readdlm("inc.dat");

Nt_zeros=80000; #zero padding
tot_field=[tot_field zeros(1,Nt_zeros)];
inc_field=[inc_field zeros(1,Nt_zeros)];

#We take the Fourier transform
totf=fft(tot_field);
incf=fft(inc_field);

#calculation of the transmittance
transm=abs.(totf)./abs.(incf);

#the following parameters are taken from the FDTD codes
Q=0.3; #Courant number
m0=12.566E-7; #permability of free space
e0=8.854187817620E-12; #permittivity of free space
c0=1/sqrt(e0*m0); #velocity of light
dz=1e-9; #space step
dt=Q*dz/c0; #time step
Nt=length(tot_field); #Number of time steps

#Plot the reflection coefficient calculated using FDTD method
df=1/(Nt*dt);
f=0:df:df*(Nt/4);
lambda_FDTD=c0./f;

#----------------------Calculation of the analytical Transmittance
#Parameters of the Drude-CP model for Gold (paper Vial 2011)
epsilon_inf=1.1431;
omega_D=1.3202e16;
gamma=1.0805e14;

A1=0.26698;
phi1=-1.2371;
Omega1=3.8711e15;
Gamma1=4.4642e14;

A2=3.0834;
phi2=-1.0968;
Omega2=4.1684e15;
Gamma2=2.3555e15;

lambda=200e-9:10e-9:1000e-9;
omega=2*pi*(c0./lambda);
Id=ones(length(omega),1);
epsilon=epsilon_inf*Id+(omega_D^2) ./(-omega.^2+im*gamma*omega) +
    A1*Omega1*(exp(im*phi1) ./((Omega1-im*Gamma1)*Id+omega) + exp(-im*phi1) ./((Omega1+im*Gamma1)*Id-omega)) +
    A2*Omega2*(exp(im*phi2) ./((Omega2-im*Gamma2)*Id+omega) + exp(-im*phi2) ./((Omega2+im*Gamma2)*Id-omega));

ri= sqrt.(epsilon);

#complex propagation constant
g=im*(omega/c0).*ri;

b0=omega./c0; #wavelength in free space

#impendances
n0=sqrt(m0/e0); #of free space
n=n0./ri; #of the medium

Z1=n0./n;
Z2=n./n0;

#the thickness of the slab
d=20e-9;
#Transmission coefficient
T=4*exp.(im*d*b0) ./ ((ones(length(n))-Z1).*(ones(length(n))-Z2).*exp.(-g*d)
+(ones(length(n))+Z1).*(ones(length(n))+Z2).*exp.(g*d));


# Plotting
#FDTD solution
plot(lambda_FDTD./1e-9,(transm[1:length(f)]).^2,color="g");
#Analytical solution
plot(lambda./1e-9,abs.(T.^2),"o");
xlabel("Wavelength (nm)")
ylabel("Transmittance")
axis([200, 1000, 0.05, 0.6])
legend(["FDTD method","Analytical"])
