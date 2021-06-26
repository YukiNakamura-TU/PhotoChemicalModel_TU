import re #Regular expression
import os #operating system interface
import matplotlib.pyplot as plt #plot library
from matplotlib import ticker, cm, colors
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np #mathematical library

#   w/o metal
hi_sigma_P_0 = [[[0.0 for i in range(11)] for lat in range(181)] for lt in range(361) ]
hi_sigma_H_0 = [[[0.0 for i in range(11)] for lat in range(181)] for lt in range(361) ]
#   w/  metal
hi_sigma_P_1 = [[[0.0 for i in range(11)] for lat in range(181)] for lt in range(361) ]
hi_sigma_H_1 = [[[0.0 for i in range(11)] for lat in range(181)] for lt in range(361) ]

# read conductivity data
for ib in range(11):
    path = './no_metal_Hill/output/conductivity/B_Jupiter_i='+str(ib)+'/hi_sigma_HP0.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            data = np.loadtxt(path, comments='!')
            for i in range(len(data)):
                y = int(data[i][0]) # latitude -90 ~ +90
                x = int(data[i][1]) # longitude 1 ~ 361
                hi_sigma_H_0[x-1][y+90][ib] = data[i][2]
                hi_sigma_P_0[x-1][y+90][ib] = data[i][3]
                #print(data[i][2],hi_sigma_H_0[x-1][y+90][ib])
        f.close()

    path = './metal_Hill/output/conductivity/B_Jupiter_i='+str(ib)+'/hi_sigma_HP0.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            data = np.loadtxt(path, comments='!')
            for i in range(len(data)):
                y = int(data[i][0]) # latitude -90 ~ +90
                x = int(data[i][1]) # longitude 1 ~ 361
                hi_sigma_H_1[x-1][y+90][ib] = data[i][2]
                hi_sigma_P_1[x-1][y+90][ib] = data[i][3]
        f.close()


sigma_P_ratio = [[0.0 for lt in range(25)] for ib in range(11)]
sigma_H_ratio = [[0.0 for lt in range(25)] for ib in range(11)]

# sigma ratio  -----------------------------------------------------------

lat = 30
ltarr = []
for lt in range(25): # 0-24
    ltarr.append(lt)
ibarr = []
for ib in range(11): # 0-10
    bfactor = 10.0**(float(ib-5)/5.0)
    ibarr.append(bfactor)

X, Y = np.meshgrid(ltarr, ibarr)

for lt in range(25): # 0-24
    for ib in range(11): # 0-10
        x = lt*15
        y = lat+90
        sigma_P_ratio[ib][lt] = hi_sigma_P_1[x][y][ib]/hi_sigma_P_0[x][y][ib]
        sigma_H_ratio[ib][lt] = hi_sigma_H_1[x][y][ib]/hi_sigma_H_0[x][y][ib]

path = './figures/sigma_ratio.dat'
with open(path, mode = 'w') as f:
    for i in range(25):
        for j in range(11):
            f.write(str(ltarr[i])+' '+str(ibarr[j])+' '+str(sigma_H_ratio[j][i])+' '+str(sigma_P_ratio[j][i])+'\n')
        f.write('\n')
f.close()


#---------------------------------------------------
#                   PLOT
#---------------------------------------------------

# Pedersen Condustance ratio ---------------
fig = plt.figure(dpi=100, figsize=(16.0, 16.0))
plt.subplots_adjust(wspace=0.6, hspace=0.5)
fs_title=20
fs_ticks=20
fs_label=20

ax1 = plt.subplot(321)
im1 = ax1.pcolormesh(X, Y, sigma_P_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax1)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax1.set_xlabel('Local time [hour]', fontsize=fs_label)
ax1.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax1.xaxis.set_ticks([0, 6, 12, 18, 24])
ax1.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax1.yaxis.set_ticks([0.1, 1.0, 10.0])
ax1.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax1.set_ylim([0.08,12])

# Colorbar
im1.set_clim(1.0,100.0)
cbar1 = fig.colorbar(im1, cax=ax_cb, ticks=[1,10,100])
cbar1.ax.set_yticklabels(['1', '10', '100'],fontsize=fs_ticks)
cbar1.set_label('Ratio', size=fs_label)

# title
ax1.set_title("$\Sigma_P$ ratio (Latitude = 30$^\circ$)", loc='center',size=fs_title)

#plt.show()
#plt.savefig("sigma_P_ratio_LAT"+str(lat)+".png")


# Hall Condustance ratio ---------------
ax2 = plt.subplot(322)
im2 = ax2.pcolormesh(X, Y, sigma_H_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax2)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax2.set_xlabel('Local time [hour]', fontsize=fs_label)
ax2.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax2.xaxis.set_ticks([0, 6, 12, 18, 24])
ax2.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax2.yaxis.set_ticks([0.1, 1.0, 10.0])
ax2.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax2.set_ylim([0.08,12])

# Colorbar
im2.set_clim(1.0,1000.0)
cbar2 = fig.colorbar(im2, cax=ax_cb, ticks=[1,10,100,1000])
cbar2.ax.set_yticklabels(['1', '10', '100', '1000'],fontsize=fs_ticks)
cbar2.set_label('Ratio', size=fs_label)

# title
ax2.set_title("$\Sigma_H$ ratio (Latitude = 30$^\circ$)", loc='center',size=fs_title)

#plt.show()
#plt.savefig("sigma_H_ratio_LAT"+str(lat)+".png")


# sigma ratio  -----------------------------------------------------------

lat = 60

for lt in range(25): # 0-24
    for ib in range(11): # 0-10
        x = lt*15
        y = lat+90
        sigma_P_ratio[ib][lt] = hi_sigma_P_1[x][y][ib]/hi_sigma_P_0[x][y][ib]
        sigma_H_ratio[ib][lt] = hi_sigma_H_1[x][y][ib]/hi_sigma_H_0[x][y][ib]

path = './sigma_ratio.dat'
with open(path, mode = 'w') as f:
    for i in range(25):
        for j in range(11):
            f.write(str(ltarr[i])+' '+str(ibarr[j])+' '+str(sigma_H_ratio[j][i])+' '+str(sigma_P_ratio[j][i])+'\n')
        f.write('\n')
f.close()


#---------------------------------------------------
#                   PLOT
#---------------------------------------------------

# Pedersen Condustance ratio ---------------
ax3 = plt.subplot(323)
im3 = ax3.pcolormesh(X, Y, sigma_P_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax3)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax3.set_xlabel('Local time [hour]', fontsize=fs_label)
ax3.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax3.xaxis.set_ticks([0, 6, 12, 18, 24])
ax3.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax3.yaxis.set_ticks([0.1, 1.0, 10.0])
ax3.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax3.set_ylim([0.08,12])

# Colorbar
im3.set_clim(1.0,100.0)
cbar3 = fig.colorbar(im3, cax=ax_cb, ticks=[1,10,100])
cbar3.ax.set_yticklabels(['1', '10', '100'],fontsize=fs_ticks)
cbar3.set_label('Ratio', size=fs_label)

# title
ax3.set_title("$\Sigma_P$ ratio (Latitude = 60$^\circ$)", loc='center',size=fs_title)

#plt.show()
#plt.savefig("sigma_P_ratio_LAT"+str(lat)+".png")


# Hall Condustance ratio ---------------
ax4 = plt.subplot(324)
im4 = ax4.pcolormesh(X, Y, sigma_H_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax4)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax4.set_xlabel('Local time [hour]', fontsize=fs_label)
ax4.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax4.xaxis.set_ticks([0, 6, 12, 18, 24])
ax4.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax4.yaxis.set_ticks([0.1, 1.0, 10.0])
ax4.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax4.set_ylim([0.08,12])

# Colorbar
im4.set_clim(1.0,1000.0)
cbar4 = fig.colorbar(im4, cax=ax_cb, ticks=[1,10,100,1000])
cbar4.ax.set_yticklabels(['1', '10', '100', '1000'],fontsize=fs_ticks)
cbar4.set_label('Ratio', size=fs_label)

# title
ax4.set_title("$\Sigma_H$ ratio (Latitude = 60$^\circ$)", loc='center',size=fs_title)

#plt.show()
#plt.savefig("sigma_ratio.png")



# sigma ratio  -----------------------------------------------------------

lat = 73

for lt in range(25): # 0-24
    for ib in range(11): # 0-10
        x = lt*15
        y = lat+90
        sigma_P_ratio[ib][lt] = hi_sigma_P_1[x][y][ib]/hi_sigma_P_0[x][y][ib]
        sigma_H_ratio[ib][lt] = hi_sigma_H_1[x][y][ib]/hi_sigma_H_0[x][y][ib]

path = './sigma_ratio.dat'
with open(path, mode = 'w') as f:
    for i in range(25):
        for j in range(11):
            f.write(str(ltarr[i])+' '+str(ibarr[j])+' '+str(sigma_H_ratio[j][i])+' '+str(sigma_P_ratio[j][i])+'\n')
        f.write('\n')
f.close()


#---------------------------------------------------
#                   PLOT
#---------------------------------------------------

# Pedersen Condustance ratio ---------------
ax5 = plt.subplot(325)
im5 = ax5.pcolormesh(X, Y, sigma_P_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax5)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax5.set_xlabel('Local time [hour]', fontsize=fs_label)
ax5.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax5.xaxis.set_ticks([0, 6, 12, 18, 24])
ax5.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax5.yaxis.set_ticks([0.1, 1.0, 10.0])
ax5.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax5.set_ylim([0.08,12])

# Colorbar
im5.set_clim(1.0,100.0)
cbar5 = fig.colorbar(im5, cax=ax_cb, ticks=[1,10,100])
cbar5.ax.set_yticklabels(['1', '10', '100'],fontsize=fs_ticks)
cbar5.set_label('Ratio', size=fs_label)

# title
ax5.set_title("$\Sigma_P$ ratio (Latitude = 73$^\circ$)", loc='center',size=fs_title)

#plt.show()
#plt.savefig("sigma_P_ratio_LAT"+str(lat)+".png")


# Hall Condustance ratio ---------------
ax6 = plt.subplot(326)
im6 = ax6.pcolormesh(X, Y, sigma_H_ratio, cmap='hot',norm=LogNorm(),shading="auto")
plt.yscale('log')
#plt.colorbar()
divider = make_axes_locatable(ax6)
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
fig.add_axes(ax_cb)
#plt.colorbar(im, cax=ax_cb, ticks=[1,10,100])

# Labels
ax6.set_xlabel('Local time [hour]', fontsize=fs_label)
ax6.set_ylabel('k=B/B$_J$ at surface', fontsize=fs_label)

# X axis
ax6.xaxis.set_ticks([0, 6, 12, 18, 24])
ax6.set_xticklabels(['0', '6', '12', '18', '24'],fontsize=fs_ticks)

# Y axis
ax6.yaxis.set_ticks([0.1, 1.0, 10.0])
ax6.set_yticklabels(['0.1', '1', '10'],fontsize=fs_ticks)
ax6.set_ylim([0.08,12])

# Colorbar
im6.set_clim(1.0,1000.0)
cbar6 = fig.colorbar(im6, cax=ax_cb, ticks=[1,10,100,1000])
cbar6.ax.set_yticklabels(['1', '10', '100', '1000'],fontsize=fs_ticks)
cbar6.set_label('Ratio', size=fs_label)

# title
ax6.set_title("$\Sigma_H$ ratio (Latitude = 73$^\circ$)", loc='center',size=fs_title)

#plt.show()
plt.savefig("figures/sigma_ratio.png")