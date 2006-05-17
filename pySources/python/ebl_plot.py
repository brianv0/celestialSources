import numarray as num
from FunctionWrapper import FunctionWrapper
import pySources
import hippoplotter as plot
from log_array import log_array

energies = log_array(100, 1e4, 1e7-1e3)
zvalues = log_array(100, 0.01, 6)

blazars = {}
blazars['Mrk421'] = 0.03
blazars['3C273'] = 0.158
blazars['3C279'] = 0.536
blazars['q0906'] = 5.47

def eblPlot(model=pySources.Primack05):
    atten = pySources.EblAtten(model)
    tau = []
    for z in zvalues:
        foo = FunctionWrapper(lambda energy : atten(energy, z))
        tau.append(foo(energies))
    tau = num.array(tau)

    disp = plot.image(tau, num.log10(energies)-3, num.log10(zvalues),
                      xname='log10(E/GeV)', yname='log10(z)', zname='tau',
                      aspect=1)
    disp.setPointRep(plot.prf.create('Contour'))
    disp.setContourLevels((0.5, 1, 2, 3, 5))
    disp.setRange('z', 0, 6)
    plot.canvas.selectDisplay(disp)
    for id in blazars:
        plot.hline(num.log10(blazars[id]))
    plot.vline(num.log10(20./2e3))
    plot.vline(num.log10(2e5/2e3))
    return disp

if __name__ == '__main__':
    eblPlot().setTitle('Primack et al. 2005')
    eblPlot(pySources.Stecker05).setTitle('Stecker et al. 2005')
