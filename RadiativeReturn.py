# RadiativeReturn.py
#
# Simple 3-body kinematics in the center-of-momentum system
# See Madison and Wilson, 2109.03281 for this angles setup.
#
# e+ e- -> f fbar gamma with various values for the fermion and difermion mass
#
import random
import math
from ROOT import TFile, TH1D, TF1, TMath
from histUtils import histInfo
from kinUtils import OpeningAngle, SupplementaryAngle

# Histogram and histogram file initialization
hcosthcm = TH1D("hcosthcm","hcosthcm; costhstar",100,-1.0,1.0)
hphicm = TH1D("hphicm","hphicm; phicm",100,0.0,2.0*math.pi)
hdev = TH1D("hdev","hdev; Deviation (ppm)",1000,0.0,200.0)
hdev2 = TH1D("hdev2","hdev2; Deviation (ppm)",1000,0.0,5000.0)
hdev3 = TH1D("hdev3","hdev3; Deviation (ppm)",1000,0.0,1.0e-3)
hpt = TH1D("hpt","hpt; pT (GeV)",250,0.0,50.0)
f  = TFile("RadiativeReturn.root", "recreate")

Nexp = 10

ECM = 250.0                    # Center of mass energy (GeV)
Mff = 91.1876                  # the invariant mass of the di-fermion (GeV)
M = 0.0
ML1 = 0.51099895000e-3         # electron mass (GeV)
ML2 = 0.1056583755             # muon mass (GeV)
ML3 = 1.77686                  # tau mass (GeV)

ilepton = 2                    # Choose lepton generation

if ilepton == 1:
    M = ML1
elif ilepton == 2:
    M = ML2
else:
    M = ML3

ECMSQ = ECM**2
XGAM = 1.0 - (Mff**2/ECMSQ)
EGAM = XGAM*ECM/2.0

print('fermion mass set to ',M)
print('xgamma set to ',XGAM)
print('mff^2/s       ',1.0-XGAM)
print('Egamma set to ',EGAM)

plab = EGAM
Elab = math.sqrt(plab**2 + Mff**2)
b=plab/Elab
g=Elab/Mff
gb=g*b
print('Boost beta:   ',b)
print('Boost gamma:  ',g)
print('Boost gam*bet ',gb)

random.seed(200)
for i in range(Nexp):
# Choose random scattering angle in the di-muon center-of-mass
# according to (1+costhcm**2) distribution by hit/miss MC
    tryAgain = True
#    costhcm = 0.0
    while tryAgain:
        costhcm = random.uniform(-1.0,1.0)
        testr = 2.0*random.uniform(0.0,1.0)
        yvalue = 1.0 + costhcm**2
        if testr < yvalue:
            tryAgain = False                 # Keep this one
    
    print('costhcm = ',costhcm)
    thcm = math.acos(costhcm)    
    phicm = random.uniform(0.0,2.0*math.pi)
    hcosthcm.Fill(costhcm)
    hphicm.Fill(phicm)

# 4-vectors oF each fermion in the di-fermion rest frame        
    E1cm = Mff/2.0
    E2cm = E1cm
    pcm = math.sqrt(E1cm**2 - M**2)
    px1cm = pcm*math.sin(thcm)*math.cos(phicm)
    py1cm = pcm*math.sin(thcm)*math.sin(phicm)
    pz1cm = pcm*math.cos(thcm)
    px2cm = -px1cm
    py2cm = -py1cm
    pz2cm = -pz1cm
    
# Now boost these into the lab (assuming a recoiling photon of energy EGAM along the -z axis). 
# This leads to a difermion system with mass of mff of course, and a momentum of EGAM.
# Boost the muons along the +z axis.    
#
    E1lab  = g*E1cm + gb*pz1cm
    px1lab = px1cm
    py1lab = py1cm
    pz1lab = gb*E1cm + g*pz1cm
    p1lab = math.sqrt(px1lab**2 + py1lab**2 + pz1lab**2)
    costh1 = pz1lab/p1lab
    pt1 = math.sqrt(px1lab**2+py1lab**2)
    
    E2lab  = g*E2cm + gb*pz2cm
    px2lab = px2cm
    py2lab = py2cm
    pz2lab = gb*E2cm + g*pz2cm
    p2lab = math.sqrt(px2lab**2 + py2lab**2 + pz2lab**2)
    costh2 = pz2lab/p2lab
    pt2 = math.sqrt(px2lab**2+py2lab**2)    
    
    E3 = EGAM
    p3 = EGAM
    px3 = 0.0
    py3 = 0.0
    pz3 = -EGAM
    costh3 = pz3/p3    
    
    print('Event   ',i)
    print('Muon1:  ',E1lab, px1lab, py1lab, pz1lab, p1lab, costh1, pt1)
    print('Muon2:  ',E2lab, px2lab, py2lab, pz2lab, p2lab, costh2, pt2)
    print('Photon: ',E3, px3, py3, pz3, p3, costh3)

# Calculate opening angles
    ang12 = OpeningAngle( px1lab, py1lab, pz1lab, px2lab, py2lab, pz2lab)
    ang23 = OpeningAngle( px2lab, py2lab, pz2lab, px3, py3, pz3)
    ang31 = OpeningAngle( px3, py3, pz3, px1lab, py1lab, pz1lab) 

# Check
    angsum = ang12 + ang23 + ang31
    print('angsum: ',angsum,' = ',angsum*180.0/math.pi,' degrees ')
    
# Calculate supplementary angles (for the triangle construction)
    sang12 = SupplementaryAngle( px1lab, py1lab, pz1lab, px2lab, py2lab, pz2lab)
    sang23 = SupplementaryAngle( px2lab, py2lab, pz2lab, px3, py3, pz3)
    sang31 = SupplementaryAngle( px3, py3, pz3, px1lab, py1lab, pz1lab) 

# Check
    sangsum = sang12 + sang23 + sang31
    print('sangsum: ',sangsum,' = ',sangsum*180.0/math.pi,' degrees ')  
    
# Now evaluate the interior angle formula
    sinesum = math.sin(sang12) + math.sin(sang23) + math.sin(sang31)
#    numerator = 2.0*math.sin(sang23)*math.sin(sang31)*(1.0 - math.cos(ang12)) #NB psi12 is NOT the opening angle
    numerator = 2.0*math.sin(sang23)*math.sin(sang31)*(1.0 + math.sin(sang12))
    
    ratio = numerator/(sinesum**2)            # evaluate equation 13 from 2209.03281
    ecmestimate = Mff/math.sqrt(ratio)        # equation 18 from 2209.03281
    deviation = 1.0e6*(ecmestimate-ECM)/ECM   # deviation in ppm
    hdev.Fill(deviation)
    hdev2.Fill(deviation) 
    hdev3.Fill(deviation)
    hpt.Fill(pt1)       
    
    print('sinesum:          ',sinesum)
    print('numerator:        ',numerator)
    print('ratio:            ',ratio)
    print('ecmest:           ',ecmestimate)
    print('deviation (ppm) = ',deviation)
         
histList = [ hcosthcm, hphicm, hdev, hdev2, hdev3, hpt]
for h in histList:
    histInfo(h)

f.Write()
