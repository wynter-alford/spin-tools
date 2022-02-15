% Order of storage in resultarray: pulse, tau, delta, coupling

function FIfH2H4 = fidelityMagnusConvert(fidFile,magFile, pulse,tau)

%% Set pulse/tau values of interest
%pulse = 0.4e-6;
%tau = 4e-6;

%% Grid setup
meshing = 30;  % based on meshing from fullAnalysis and fullMagnus

maxTau = 14;
maxPulse = 1.4;
maxDelta = 2000;
maxCoupling = 4000;

taus=zeros(meshing,1);
pulses=zeros(meshing,1);
deltas = zeros(meshing,1);
couplings=zeros(meshing,1);

for i=1:meshing
    taus(i)=(maxTau/meshing)*i*10^(-6);
    pulses(i)=(maxPulse/meshing)*i*10^(-6);
    deltas(i)=(-1) * maxDelta + 2*(maxDelta/meshing)*i;
    couplings(i)=(maxCoupling/meshing)*i;
end

%% Find chosen pulse/tau values in the array

i = 1;
while pulses(i)<pulse
    i = i + 1;
end

pulseIndex = i;

i = 1;
while taus(i)<tau
    i = i + 1;
end

tauIndex = i;

%% Load results, create a 2D array for chosen pulse & tau

ff = load(fidFile);
fidelity4D = ff.resultArray;

fidelity2D = zeros(meshing,meshing);

hh = load(magFile);
fn = fieldnames(hh);
h2D = getfield(hh,fn{1});

for i=1:meshing
    for j=1:meshing
        fidelity2D(i,j) = fidelity4D(pulseIndex,tauIndex,i,j);
    end
end


%% Convert 2D arrays to matched 1D arrays
mm = meshing^2;
fidelities = zeros(mm,1);
lFidelities = zeros(mm,1);
hterms = zeros(mm,1);

for i = 1:meshing
    for j = 1:meshing  
        position = 30*(i-1) + j;
        fidelities(position) = fidelity2D(i,j);
        hterms(position) = h2D(i,j);
        
        lFidelities(position) = -log10(1 - fidelities(position));
    end
end

FIfH2H4 = struct('fidelities',fidelities,'lFidelities',lFidelities,'hterms',10.^hterms);

end
