#!~/anaconda3/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


def update_lines(dataLines, lines):
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2,])
        line.set_3d_properties(data[2,:])
    return lines


def main(filename):
    paths_df = pd.read_csv(filename, header=None)
    ncolumns = len(paths_df.columns)
    nparticles = int((ncolumns-1)/3)

    particlepos_dict = {}
    for i in range(nparticles):
        particlepos_dict['%i'%i] = paths_df.iloc[:,(3*i+1):(3*i+4)]
        data = paths_df.iloc[:,1:4]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(nparticles):
        x = particlepos_dict['%i'%i].iloc[:,0]
        y = particlepos_dict['%i'%i].iloc[:,1]
        z = particlepos_dict['%i'%i].iloc[:,2]
        ax.plot(x, y, z)

    ax.set_xlim3d([-1e-6, 1e-6])
    ax.set_xlabel('x')

    ax.set_ylim3d([-1e-6, 1e-6])
    ax.set_ylabel('y')

    ax.set_zlim3d([-1e-6, 1e-6])
    ax.set_zlabel('z')

    plt.show()

    """ax = p3.Axes3D(fig)
    lines = []
    for index, dat in data.iterrows():
        lines.append(ax.plot(dat.iloc[0], dat.iloc[1], dat.iloc[2])[0])
    line_ani = animation.FuncAnimation(fig, update_lines, fargs=(data, lines),
                                        interval=50, blit=False)
    """




if __name__ == '__main__':
    #main(sys.api_versionargv([1]))
    main('../output0000000.csv')