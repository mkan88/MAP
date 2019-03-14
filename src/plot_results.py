#!~/anaconda3/bin/python

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import ticker


def update_lines(i, nparticles, particlepos_dict, lines):
    #len(lines) == on of particles
    for n, line in enumerate(lines):
        x = particlepos_dict['%i'%n].iloc[i:(i+3),0]
        y = particlepos_dict['%i'%n].iloc[i:(i+3),1]
        z = particlepos_dict['%i'%n].iloc[i:(i+3),2]
        line.set_data([x,y])
        line.set_3d_properties(z)
    return lines


def main(filename, animate=False, beforeafter=False):
    paths_df = pd.read_csv(filename, header=None)
    ncolumns = len(paths_df.columns)
    nparticles = int((ncolumns-1)/3)

    particlepos_dict = {}
    for n in range(nparticles):
        particlepos_dict['%i'%n] = paths_df.iloc[:,(3*n+1):(3*n+4)]

    if animate:
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        lines = []
        for n in range(nparticles):
            x = particlepos_dict['%i'%n].iloc[:3,0]
            y = particlepos_dict['%i'%n].iloc[:3,1]
            z = particlepos_dict['%i'%n].iloc[:3,2]
            lines.append(ax.plot(x, y, z)[0])

        ax.set_xlim3d([-1.5e-6, 1.5e-6])
        ax.set_xlabel('x')

        ax.set_ylim3d([-1.5e-6, 1.5e-6])
        ax.set_ylabel('y')

        ax.set_zlim3d([-1.5e-6, 1.5e-6])
        ax.set_zlabel('z')

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

        line_ani = animation.FuncAnimation(fig, update_lines, fargs=(nparticles, particlepos_dict, lines),
                                        frames=range(0, len(particlepos_dict['%i'%0])-4),
                                        interval=1, blit=False)

        plt.show()
        # plt.hold(False) #Prevents loop
    
    elif beforeafter:
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        lines = []
        for n in range(nparticles):
            x = particlepos_dict['%i'%n].iloc[:3,0]
            y = particlepos_dict['%i'%n].iloc[:3,1]
            z = particlepos_dict['%i'%n].iloc[:3,2]
            lines.append(ax.plot(x, y, z)[0])

        ax.set_xlim3d([-1.5e-6, 1.5e-6])
        ax.set_xlabel('x')

        ax.set_ylim3d([-1.5e-6, 1.5e-6])
        ax.set_ylabel('y')

        ax.set_zlim3d([-1.5e-6, 1.5e-6])
        ax.set_zlabel('z')

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

        line_ani = animation.FuncAnimation(fig, update_lines, fargs=(nparticles, particlepos_dict, lines),
                                        frames=[1, len(particlepos_dict['%i'%0])-4],
                                        interval=500, blit=False)

        plt.show()

    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for n in range(nparticles):
            x = particlepos_dict['%i'%n].iloc[:,0]
            y = particlepos_dict['%i'%n].iloc[:,1]
            z = particlepos_dict['%i'%n].iloc[:,2]
            ax.plot(x, y, z)

        ax.set_xlim3d([-1.5e-6, 1.5e-6])
        ax.set_xlabel('x')

        ax.set_ylim3d([-1.5e-6, 1.5e-6])
        ax.set_ylabel('y')

        ax.set_zlim3d([-1.5e-6, 1.5e-6])
        ax.set_zlabel('z')

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

        plt.show()

        




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
            """Plots particle paths.
            """)
    parser.add_argument('filename', type=str, default='../output0000000.csv')
    parser.add_argument('--animate', '-a', action='store_true')
    parser.add_argument('--beforeafter', '-b', action='store_true')
    args = parser.parse_args()

    main(args.filename, animate=args.animate, beforeafter=args.beforeafter)