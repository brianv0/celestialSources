import numarray as num
import pySources

params = "0.1, 0., 1e4, $(GENERICSOURCESROOT)/data/testTemplate.dat"

my_src = pySources.SpectralTransient(params)

grb_src = pySources.GRBmanager("60, 600, 30")

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

plot.Histogram(nt, 'time')
plot.Histogram(nt, 'energy', xlog=1, ylog=1)
plot.XYHist(nt, 'glon', 'glat')

nt2 = plot.newNTuple(generateData(grb_src, 3600),
                     ('time', 'energy', 'glon', 'glat'))

plot.Histogram(nt2, 'time')
plot.Histogram(nt2, 'energy', xlog=1, ylog=1)
plot.XYHist(nt2, 'glon', 'glat')
