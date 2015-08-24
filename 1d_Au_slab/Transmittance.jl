# Calculation of the transmittance
# using the proposed FDTD method vs analytical

#Code written in Julia programming language
#Konstantinos Prokopidis, August, 24,2015

#clear all variables
workspace();

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
transm=abs(totf)./abs(incf);

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
lamda=c0./f;

#Calculation of the analytical Transmittance
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
epsilon=epsilon_inf+omega_D^2./(-omega.^2+im*gamma*omega)+
    A1*Omega1*(exp(im*phi1)./(Omega1+omega-im*Gamma1)+exp(-im*phi1)./(Omega1-omega+im*Gamma1))+
    A2*Omega2*(exp(im*phi2)./(Omega2+omega-im*Gamma2)+exp(-im*phi2)./(Omega2-omega+im*Gamma2));

ri= sqrt(epsilon);

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
T=4.*exp(im.*b0.*d)./((1-Z1).*(1-Z2).*exp(-g.*d)+(1+Z1).*(1+Z2).*exp(g.*d));


# Plotting
plot(lamda./1e-9,transm[1:length(f)].^2);
draw
plot(lambda./1e-9,abs(T.^2),"o");
xlabel("Wavelength (nm)");
ylabel("Transmission");
axis([200, 1000, 0.05, 0.6]);
legend(["FDTD method","Analytical"]);









