import math

def OpeningAngle(px1, py1, pz1, px2, py2, pz2):
#
# Opening angle in radians between two 3-vectors
#    
    p1  = math.sqrt(px1**2 + py1**2 + pz1**2)
    p2  = math.sqrt(px2**2 + py2**2 + pz2**2)    
    costh = (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2)
    angle = math.acos(costh)
    
    print('Calculated angle in radians is ',angle,'( ',angle*180.0/math.pi,' deg )')
    
    return angle
    
def SupplementaryAngle(px1, py1, pz1, px2, py2, pz2):
#
# Supplementary angle of the opening angle in radians between two 3-vectors
#    
    p1  = math.sqrt(px1**2 + py1**2 + pz1**2)
    p2  = math.sqrt(px2**2 + py2**2 + pz2**2)    
    costh = (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2)
    angle = math.acos(costh)
    
    sangle = math.pi - angle
    print('Supplementary angle in radians is ',sangle,'( ',sangle*180.0/math.pi,' deg )')
    
    return sangle
