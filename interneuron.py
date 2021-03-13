from neuron import h
from neuron.units import ms,mV
from cell import Cell

class Interneuron(Cell):
    name = 'Interneuron'
    def _set_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.soma.L= 10/3.14
        self.soma.diam = 10 #SA= 100 um2
        
    def _set_biophysics(self):
        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1