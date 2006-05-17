import numarray as num
import pySources
from FunctionBaseDecorator import FunctionBaseDecorator

params = "flux=0.1, tstart=0., tstop=1e4, templateFile=testTemplate.dat, useLogParabola=1, z=1, eblModel=20"

my_src = pySources.SpectralTransient(params)

#grb_src = pySources.GRBmanager("60, 600, 30")

def generateData(src, sim_time):
    time = 0
    arrTimes = []
    energies = []
    glons = []
    glats = []
    while time < sim_time:
        dt = src.interval(time)
        time += dt
        arrTimes.append(time)
        energies.append(src.energy(time))
        try:
            l, b = src.dir(energies[-1])
        except AttributeError:
            l, b = 0, 0
        glons.append(l)
        glats.append(b)
    return (num.array(arrTimes), num.array(energies),
            num.array(glons), num.array(glats))

import hippoplotter as plot

nt = plot.newNTuple(generateData(my_src, 1e4),
                    ('time', 'energy', 'glon', 'glat'))

plot.Histogram(nt, 'time', xrange=(0, 1e4))
spectrum = plot.Histogram(nt, 'energy', xlog=1, ylog=1)
plot.XYHist(nt, 'glon', 'glat')

@FunctionBaseDecorator
def f(e, e1=300, a=1, b=2, k=1):
    return k*(e/e1)**(-(a + b*num.log(e/e1)))

spectrum.addFunction(f)

#nt2 = plot.newNTuple(generateData(grb_src, 3600),
#                     ('time', 'energy', 'glon', 'glat'))
#
#plot.Histogram(nt2, 'time')
#plot.Histogram(nt2, 'energy', xlog=1, ylog=1)
#plot.XYHist(nt2, 'glon', 'glat')
