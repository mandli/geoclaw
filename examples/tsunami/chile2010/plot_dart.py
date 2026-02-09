"""
Compare DART buoy results to NOAA results.
"""

import matplotlib.pyplot as plt
from pyclaw.plotters.data import ClawPlotData

# from matplotlib import image

dartpng = plt.image.imread('dart32412_comp-2.png')
plt.figure(400)
plt.clf()
plt.imshow(dartpng)
plt.hold(True)

# origin of plot in pixels on image:
t0 = 195.
y0 = 221.
tscale = 120.0
yscale = -446.66666666666669
plotdata = ClawPlotData()
plotdata.outdir = "_output"
gaugedata = plotdata.getgauge(32412)
t = gaugedata.t
eta = gaugedata.q[:,3]
t = t0 + tscale * (t / 3600.)
y = y0 + yscale*eta

plt.plot(t,y,'b')
plt.plot([160,225],[550,550],'b')
plt.text(240,560,'GeoClaw (added to NOAA original)',fontsize=8)
plt.text(100,600,'NOAA Original from http://nctr.pmel.noaa.gov/chile20100227/dart32412_comp-2.pdf',fontsize=8)


plt.xlim([50,1300])
plt.ylim([640,0])
plt.axis('off')

plt.savefig('dart.png')
