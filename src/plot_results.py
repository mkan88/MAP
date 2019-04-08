#!~/anaconda3/bin/python

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import ticker


def update_lines(i, nparticles, particlepos_dict, time_array, lines, time_text=None):
    if time_text is not None:
        time_text.set_text(str(time_array[i]))
    #len(lines) == on of particles
    for n, line in enumerate(lines):
        x = particlepos_dict['%i'%n].iloc[i:(i+2),0]
        y = particlepos_dict['%i'%n].iloc[i:(i+2),1]
        z = particlepos_dict['%i'%n].iloc[i:(i+2),2]
        line.set_data([x,y])
        line.set_3d_properties(z)

    return lines


def plot_particle_path(filename, animate=False, frameskip=1, beforeafter=False):
    paths_df = pd.read_csv(filename, header=None)
    nrows = len(paths_df.index)
    ncolumns = len(paths_df.columns)
    nparticles = int((ncolumns-1)/3)

    time_array = paths_df.iloc[:,0].values

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

        time_text = ax.text(1.6e-6, 1.6e-6, 1.6e-6, str(particlepos_dict['0'].iloc[0,0]))

        #Append adds last frame to frames
        frame_array = np.append(np.arange(0, nrows-1, frameskip), nrows-1)

        line_ani = animation.FuncAnimation(fig, update_lines,
                                        fargs=(nparticles, particlepos_dict, time_array, lines, time_text),
                                        frames=frame_array,
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

        line_ani = animation.FuncAnimation(fig, update_lines,
                                        fargs=(nparticles, particlepos_dict, time_array, lines),
                                        frames=[1, nrows-2],
                                        interval=500, blit=False)

        plt.show()

    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        xmin = 1
        xmax = -1
        ymin = 1
        ymax = -1
        zmin = 1
        zmax = -1
        for n in range(nparticles):
            x = particlepos_dict['%i'%n].iloc[:,0]
            y = particlepos_dict['%i'%n].iloc[:,1]
            z = particlepos_dict['%i'%n].iloc[:,2]
            ax.plot(x, y, z)

            if x.min() < xmin:
                xmin = x.min()
            if x.max() > xmax:
                xmax = x.max()
            if y.min() < ymin:
                ymin = y.min()
            if y.max() > ymax:
                ymax = y.max()
            if z.min() < zmin:
                zmin = z.min()
            if z.max() > zmax:
                zmax = z.max()
            
        ax.set_xlim3d([xmin, xmax])
        ax.set_xlabel('x')

        ax.set_ylim3d([ymin, ymax])
        ax.set_ylabel('y')

        ax.set_zlim3d([zmin, zmax])
        ax.set_zlabel('z')

        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

        plt.show()


def plot_com_path(filename):
    df = pd.read_csv(filename, header=None)
    nrows = len(df.index)
    ncolumns = len(df.columns)
    nparticles = int((ncolumns-1)/3)

    particlepos_array = np.zeros((nrows, nparticles, 3), dtype=np.float)
    for i in range(nrows):
        for j in range(nparticles):
            for k in range(3):
                particlepos_array[i,j,k] = df.iloc[i,(3*j+k+1)]

    compos_array = np.zeros((nrows, 3), dtype=np.float)
    for i in range (nrows):
        compos_array[i][0] = np.average(particlepos_array[i,:,0])
        compos_array[i][1] = np.average(particlepos_array[i,:,1])
        compos_array[i][2] = np.average(particlepos_array[i,:,2])
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(compos_array[:,0], compos_array[:,1], compos_array[:,2])

    ax.set_xlim3d([compos_array[:,0].min(), compos_array[:,0].max()])
    ax.set_xlabel('x')

    ax.set_ylim3d([compos_array[:,1].min(), compos_array[:,1].max()])
    ax.set_ylabel('y')

    ax.set_zlim3d([compos_array[:,2].min(), compos_array[:,2].max()])
    ax.set_zlabel('z')

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

    plt.show()


def plot_time(pos_filename, angle_filename, time):
    paths_df = pd.read_csv(pos_filename, header=None)
    angle_df = pd.read_csv(angle_filename, header=None)
    nrows = len(paths_df.index)
    ncolumns = len(paths_df.columns)
    nparticles = int((ncolumns-1)/3)

    time_array = paths_df.iloc[:,0].values
    timeindex = (np.abs(time_array-time)).argmin()

    particlepos_dict = {}
    particleangle_dict = {}
    for n in range(nparticles):
        particlepos_dict['%i'%n] = paths_df.iloc[:,(3*n+1):(3*n+4)]
        particleangle_dict['%i'%n] = angle_df.iloc[:,(3*n+1):(3*n+4)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xmin = 1
    xmax = -1
    ymin = 1
    ymax = -1
    zmin = 1
    zmax = -1
    vectorlength = 1e-7
    for n in range(nparticles):
        x = particlepos_dict['%i'%n].iloc[timeindex,0]
        y = particlepos_dict['%i'%n].iloc[timeindex,1]
        z = particlepos_dict['%i'%n].iloc[timeindex,2]
        alpha = particleangle_dict['%i'%n].iloc[timeindex,0]
        beta = particleangle_dict['%i'%n].iloc[timeindex,1]
        # gamma = particleangle_dict['%i'%n].iloc[timeindex,2]
        u = vectorlength*np.cos(alpha)*np.sin(beta)
        v = vectorlength*np.sin(alpha)*np.sin(beta)
        w = vectorlength*np.cos(beta)
        ax.quiver(x, y, z, u, v, w)

        if x < xmin:
            xmin = x
        if x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        if y > ymax:
            ymax = y
        if z < zmin:
            zmin = z
        if z > zmax:
            zmax = z

    ax.set_xlim3d([xmin, xmax])
    ax.set_xlabel('x')

    ax.set_ylim3d([ymin, ymax])
    ax.set_ylabel('y')

    ax.set_zlim3d([zmin, zmax])
    ax.set_zlabel('z')

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))

    plt.show()


def plot_velocity_graph(filename):
    df = pd.read_csv(filename, header=None)
    nrows = len(df.index)
    ncolumns = len(df.columns)
    nparticles = int((ncolumns-1)/3)

    time_array = df.iloc[:,0].values
    deltat = time_array[1]-time_array[0]

    particlepos_dict = {}
    for n in range(nparticles):
        particlepos_dict['%i'%n] = df.iloc[:,(3*n+1):(3*n+4)]

    particlevelocity_array = np.zeros(nrows-1, dtype=np.float)
    particlevelocity_temp = np.zeros(nparticles, dtype=np.float)
    for t in range(nrows-1):
        for n in range(nparticles):
            particlevelocity_temp[n] = np.sqrt((particlepos_dict['%i'%n].iloc[t+1,0]-particlepos_dict['%i'%n].iloc[t,0])**2
                                    +(particlepos_dict['%i'%n].iloc[t+1,1]-particlepos_dict['%i'%n].iloc[t,1])**2
                                    +(particlepos_dict['%i'%n].iloc[t+1,2]-particlepos_dict['%i'%n].iloc[t,2])**2)/deltat
        particlevelocity_array[t] = np.average(particlevelocity_temp)

    # Skip last time steap 
    plt.plot(time_array[:-1], particlevelocity_array)
    # plt.errorbar(time_array, alphaavg_list, yerr=alphastd_list, linestyle='None')
    plt.title(r'$\bar{v}$ versus Time')
    plt.xlabel('Time /s')
    plt.ylabel(r'$\bar{v}$')
    plt.show()
    
    
def plot_2d_graph(filename):
    df = pd.read_csv(filename, header=None)
    nrows = len(df.index)
    ncolumns = len(df.columns)
    nparticles = int((ncolumns-1)/3)

    particleangle_array = np.zeros((nrows, nparticles, 3), dtype=np.float)
    for i in range(nrows):
        for j in range(nparticles):
            for k in range(3):
                particleangle_array[i,j,k] = df.iloc[i,(3*j+k+1)]

    alphaavg_list = np.zeros((nrows,), dtype=np.float)
    alphastd_list = np.zeros((nrows,), dtype=np.float)
    betaavg_list = np.zeros((nrows,), dtype=np.float)
    betastd_list = np.zeros((nrows,), dtype=np.float)
    gammaavg_list = np.zeros((nrows,), dtype=np.float)
    gammastd_list = np.zeros((nrows,), dtype=np.float)

    for i in range(nrows):
        alphaavg_list[i] = np.average(particleangle_array[i,:,0])
        alphastd_list[i] = np.std(particleangle_array[i,:,0], ddof=1)
        betaavg_list[i] = np.average(particleangle_array[i,:,1])
        betastd_list[i] = np.std(particleangle_array[i,:,1], ddof=1)
        gammaavg_list[i] = np.average(particleangle_array[i,:,2])
        gammastd_list[i] = np.std(particleangle_array[i,:,2], ddof=1)

    if 'angle' in filename:
        q1 = r'$\alpha$'
        q2 = r'$\beta$'
        q3 = r'$\gamma$'
    elif 'force' in filename:
        q1 = r'$F_x$'
        q2 = r'$F_y$'
        q3 = r'$F_z$'

    plt.plot(df.iloc[:,0], alphaavg_list)
    plt.errorbar(df.iloc[:,0], alphaavg_list, yerr=alphastd_list, linestyle='None')
    plt.title(r'%s versus Time'%q1)
    plt.xlabel('Time /s')
    plt.ylabel(r'%s'%q1)
    plt.show()

    plt.plot(df.iloc[:,0], betaavg_list)
    plt.errorbar(df.iloc[:,0], betaavg_list, yerr=betastd_list, linestyle='None')
    plt.title(r'%s versus Time'%q2)
    plt.xlabel('Time /s')
    plt.ylabel(r'%s'%q2)
    plt.show()

    plt.plot(df.iloc[:,0], gammaavg_list)
    plt.errorbar(df.iloc[:,0], gammaavg_list, yerr=gammastd_list, linestyle='None')
    plt.title(r'%s versus Time'%q3)
    plt.xlabel('Time /s')
    plt.ylabel(r'%s'%q3)
    plt.show()


def plot_rms(filename):
    df = pd.read_csv(filename, header=None)
    nrows = len(df.index)
    ncolumns = len(df.columns)
    
    time_array = df.iloc[:,0].values
    rms_array = df.iloc[:,1].values
    diff_array = df.iloc[:,2].values #diffusionMatrix[0]*currentTime??

    # plt.plot(time_array, rms_array)
    # plt.xlabel('Time /s')
    # plt.show()

    plt.plot(time_array, diff_array)
    plt.xlabel('Time /s')
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
            """Plots results of MAP.
            """)
    parser.add_argument('filename', type=str)
    parser.add_argument('angle_filename', type=str, nargs='?', default=None)
    parser.add_argument('--animate', '-a', action='store_true',
                        help="Animate particle path")
    parser.add_argument('--frameskip', '-s', type=int,
                        help="Animate every n steps")
    parser.add_argument('--beforeafter', '-b', action='store_true',
                        help="Animate the first and last snapshot of particle positions")
    parser.add_argument('--centreofmasspath', '-c', action='store_true',
                        help="Plots path of centre of mass")
    parser.add_argument('--time', '-t', type=float,
                        help="Plot specific time (rounded down if not found)")
    parser.add_argument('--velocity', '-v', action='store_true',
                        help="Plot average velocity")
    args = parser.parse_args()

    if 'angle_output' in args.filename:
        plot_2d_graph(args.filename)
    elif 'forces_output' in args.filename:
        plot_2d_graph(args.filename)
    elif 'rms_output' in args.filename:
        plot_rms(args.filename)
    elif 'details' in args.filename:
        pass
    else:
        if args.centreofmasspath:
            plot_com_path(args.filename)
        elif args.time:
            plot_time(args.filename, args.angle_filename, time=args.time)
        elif args.velocity:
            plot_velocity_graph(args.filename)
        else:
            plot_particle_path(args.filename, animate=args.animate, frameskip=args.frameskip, beforeafter=args.beforeafter)

