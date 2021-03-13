from neuron import h
from neuron.units import ms,mV

class Cell:
    def __init__(self,c_id,x,y,z,theta):
        self.c_id = c_id
        self.x=self.y=self.z=0
        self._set_morphology()
        self.all = self.soma.wholetree()
        self._set_biophysics()
        h.define_shape()       
        self._set_position(x,y,z)
        self._rotate_z(theta)
        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self.c_id)
    
    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
    def _rotate_z(self, theta):
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))