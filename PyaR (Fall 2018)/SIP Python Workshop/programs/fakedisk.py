"""
NAME: fakedisk.py

PURPOSE: create a fake disk consisting of perfectly cold disks
at different layers. 

"""

import numpy as np
import pdb

class fakedisk:
    def __init__(self, incl = 77, pa0 = 37.7, v = 220., N = 1000, shape = 'constant', scaleFactor = 1,
                 dispersion = 0):
        """generate a single layer perfectly cold disk
        with inclination and pa as specified.

        incl = inclination (degrees)
        pa0 = position angle of major axis (degrees)
        v = rotation speed
        N = number of points
        shape = shape of rotation curve. Only "constant" currently supported
        scale_factor = radius of disk in kpc
        dispersion = intrinsic radial velocity dispersion
        """
        self.dispersion = dispersion
        self.incl = incl * np.pi / 180.#Inclination of disk, degrees
        self.pa0 = pa0  * np.pi / 180. #position angle of major axis of fdisk, degrees
        self.vrot = v #rotation speed of disk (assume all 1 # for now)
        self.N = N #number of stars in this layer
        self.vrotshape = shape
        self.scaleFactor = scaleFactor
        self.z = np.zeros(0)
        self.x = np.zeros(0)
        self.y = np.zeros(0)
        self.vx = np.zeros(0)
        self.vy = np.zeros(0)
        self.vz = np.zeros(0)
        self.dpa = np.zeros(0)
        self.rdisk = np.zeros(0)
        self.original_dpa = np.zeros(0)
        
        self.add_layer()
        
        return

    def rotate(self, x1, x2, rotAngle):
        #private method to rotate axes. Angle is in radians.

        x1New = x1 * np.cos(rotAngle) - x2 * np.sin(rotAngle)
        x2New = x1 * np.sin(rotAngle) + x2 * np.cos(rotAngle)

        return x1New, x2New

    def rotationCurve(self, r, shape = 'linear'):
        if shape == 'constant':
            v = self.vrot * np.ones_like(r)
        else:
            print('Not yet configured for other rotation curve shapes')
            
        return v
    def add_layer(self, z0 = 0):
        """Public method to add a layer of stars at height z above the midplane. 
        """

        #positions of new layer. 
        x, y = np.random.uniform(low = -self.scaleFactor, high =  self.scaleFactor, size = (2, self.N))
       
        #pull stars within a circle of radius 1. eliminate stars in crowded center
        r = x*x + y*y
        ok = (r <= self.scaleFactor**2) 
        x = x[ok]
        y = y[ok]
        r = r[ok]
        z = z0 * np.ones(len(x))

        #add coordinates for "original" positions.
        self.rdisk = np.concatenate((self.rdisk, np.sqrt(r)))
        pa = np.arctan2(y, -x) - np.pi / 2.
        dpa = pa - np.pi/2.  #this is fine. 
        self.original_dpa = np.concatenate((self.original_dpa, dpa))
        
        #velocities of new layer. sign is somewhat ambiguous
        if self.dispersion == 0: 
            vrot = self.vrot * np.ones_like(r) #phi component
            vz = np.zeros_like(r) #z component
            vr = np.zeros_like(r) #r component
        else: #have sigmaR = sigmaPHI = 2 * sigmaZ
            vrot = np.random.normal(loc = self.vrot, scale = self.dispersion, size = len(r))
            vr = np.random.normal(loc = 0, scale = self.dispersion, size = len(r))
            vz = np.random.normal(loc = 0, scale = self.dispersion/2., size = len(r))
        vx = vrot * np.cos(pa) - vr * np.sin(pa)
        vy = vrot * np.sin(pa) + vr * np.cos(pa)
        

        #rotate
        y, z = self.rotate(y, z, np.pi/2. - self.incl)
        x, z = self.rotate(x, z, np.pi/2. -self.pa0)

        vy, vz = self.rotate(vy, vz, np.pi/2. - self.incl)
        vx, vz = self.rotate(vx, vz,  np.pi/2. -self.pa0)

        #relative position angle in radians. Origin is at (0, 1)
        pa = np.arctan2(z, -x) - np.pi / 2. #this is identical to the original pa. 
        #dpa = pa - self.pa0 #weird that i have to subtract this again. 
        dpa[dpa < 0] += 2 * np.pi
        dpa[dpa > 2 * np.pi] -= 2 * np.pi

        
        #concatenate. for now, just save current coords. 
        self.x = np.concatenate((self.x, x))
        self.y = np.concatenate((self.y, y))
        self.z = np.concatenate((self.z, z))
        self.vx = np.concatenate((self.vx, vx))
        self.vy = np.concatenate((self.vy, vy))
        self.vz = np.concatenate((self.vz, vz))
        self.dpa = np.concatenate((self.dpa, dpa))

        return


    def get_sv(self, index, radius):
        """Return LOS component of dispersion for stars near star index
        """

        d = (self.x - self.x[index])**2 + (self.z - self.z[index])**2
        indices = (d < radius**2)
        v = self.vy[indices]
        sv = np.std(v)
        #get uncertainty
        mu = np.average(v)
        I11 = np.sum(1/sv**2)
        I22 = I11 + np.sum(3 * (v - mu)**2/sv**4+ 2/sv**2)
        I12 = np.sum(2 * (v - mu)/sv**3)

        sig_sv = np.sqrt(I11/(I11 * I22 - I12**2))
        return sv, sig_sv



    def get_sv_coords(self, coords, radius):
        """Return LOS component of dispersion for stars near coords = (xi, eta)
        """
        d = (self.x - coords[0])**2 + (self.z - coords[1])**2
        indices = (d < radius**2)

        v = self.vy[indices]
        return np.std(v), len(v)


    
    
