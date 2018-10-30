#!/usr/bin/env python

# coding: utf-8
import sys
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from Utils import CSVreader, Pathname

## Classes

class Plot():
    outfile = None
    imgformat = "png"
    title = None
    xlabel = None
    ylabel = None
    xsize = 10
    ysize = 8

    def parseCommonArgs(self, args):
        other = []
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-t":
                self.title = a
                prev = ""
            elif prev == "-xl":
                self.xlabel = a
                prev = ""
            elif prev == "-yl":
                self.ylabel = a
                prev = ""
            elif prev == "-xs":
                self.xsize = float(a)
                prev = ""
            elif prev == "-ys":
                self.ysize = float(a)
                prev = ""
            elif prev == "-f":
                self.imgformat = a
                prev = ""
            elif a in ["-o", "-t", "-xl", "-yl", "-f", "-xs", "-ys"]:
                prev = a
            else:
                other.append(a)
        return other
                
class DensityPlot(Plot):
    infile = None
    log = False
    hasHeader = None
    skipRows = None
    pointSize = 50
    cx = 4
    cy = 3

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage()
        args = self.parseCommonArgs(args)
        prev = ""
        for a in args:
            if prev == "-cx":
                self.cx = int(a)
                prev = ""
            elif prev == "-cy":
                self.cy = int(a)
                prev = ""
            elif prev == "-s":
                self.skipRows = int(a)
                prev = ""
            elif prev == "-p":
                self.pointSize = int(a)
                prev = ""
            elif a in ["-cx", "-cy", "-s", "-p"]:
                prev = a
            elif a == "-l":
                self.log = np.log(2)
            elif a == "-l10":
                self.log = np.log(10)
            elif a == "-t":
                self.hasHeader = 0
            elif self.infile is None:
                self.infile = a
        return (self.infile and self.outfile)

    def run(self):
        df = pd.read_table(self.infile, header=self.hasHeader, skiprows=self.skipRows)
        sys.stderr.write("Data file read.\n")

        if self.log:
            df[self.cy] = df[self.cy].apply(lambda x: np.log(x+1) / self.log)
            df[self.cx] = df[self.cx].apply(lambda x: np.log(x+1) / self.log)
            min_x = min(df[self.cx])-0.1
            min_y = min(df[self.cy])-0.1

        xy = np.vstack([df[self.cx], df[self.cy]])
        z = gaussian_kde(xy)(xy)

        idx = z.argsort()
        x, y, z = df[self.cx][idx], df[self.cy][idx], z[idx]

        fig, ax = plt.subplots(figsize=(self.xsize, self.ysize))
        ax.scatter(x, y, c=z, s=self.pointSize, edgecolor='')
        if self.log:
            plt.xlim(xmin=min_x)
            plt.ylim(ymin=min_y)
        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c="0.3")
        if self.xlabel:
            plt.xlabel(self.xlabel)
        if self.ylabel:
            plt.ylabel(self.ylabel)
        plt.savefig(self.outfile, format=self.imgformat)

    def usage(self):
        sys.stdout.write("""polyplotter.py dscatt - Draw density scatterplots of paired data.

Usage: polyplotter.py dscatt [options] datafile imgfile

Read data from two columns of file `datafile' and draw a density heatmap
of their scatterplot to `imgfile'. 

Options related to input data:

  -cx C | Use column C for X axis coordinates (default: {}).
  -cy C | Use column C for Y axis coordinates (default: {}).
  -s  S | Skip S rows from top of input file (default: {}).
  -l    | If supplied, log-transform data (base 2).
  -l10  | If supplied, log-transform data (base 10).

Graphical options:

  -xs S | Set X dimension of image to S inches (default: {}).
  -ys S | Set Y dimension of image to S inches (default: {}).
  -xl L | Set X axis label to L.
  -yl L | Set Y axis label to L.
  -p P  | Set dot size to P (default: {}).
  -f F  | Set output image format to F (default: {}).

""".format(DensityPlot.cx, DensityPlot.cy, DensityPlot.skipRows, DensityPlot.xsize, DensityPlot.ysize, DensityPlot.pointSize, DensityPlot.imgformat))
        sys.exit(1)


class MethylHist(Plot):
    infile = None

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage()
        args = self.parseCommonArgs(args)
        prev = ""
        for a in args:
            if self.infile == None:
                self.infile = a
        return (self.infile and self.outfile)

    def makeHistogramFromColumn(self, column, edges, normalize=False):
        nbins = len(edges) - 1
        bins = np.zeros(nbins)
        for line in CSVreader(self.infile):
            x = float(line[column])
            if x < edges[0]:
                continue
            for i in range(nbins):
                if x < edges[i+1]:
                    bins[i] += 1
                    break
        if normalize:
            bins = bins / np.sum(bins)
        return bins

    def run(self):
        path = Pathname(self.infile)
        edges = np.linspace(0.0, 1.0, num=11, endpoint=True)
        bins = self.makeHistogramFromColumn(3, edges, normalize=True)
        xc = [ int((x - 0.05) * 100) for x in edges[1:] ]
        fig, ax = plt.subplots(1, 1, figsize=(self.xsize, self.ysize))
        ax.bar(xc, bins, width=8)
        if self.title:
            ax.set_title(self.title)
        else:
            ax.set_title("{} - histogram of methylation values".format(path.name))
        ax.set_xlabel("% Methylation")
        ax.set_ylabel("% Sites")
        fig.savefig(self.outfile)

    def usage(self):
        sys.stdout.write("""polyplotter.py mhist - Draw histogram of methylation data.

Usage: polyplotter.py mhist [options] datafile

Read data from file `datafile' and draw a histogram.

Options related to input data:

- TODO

Graphical options:

- TODO

""".format(MethylHist.infile))
        sys.exit(1)

class DistPlot(Plot):
    infile = None
    log = False
    plot_type = None
    skipRows = None
    pointSize = 50
    index_col = 0
    cx = 4
    cy = 3

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage()
        args = self.parseCommonArgs(args)
        prev = ""
        for a in args:
            if prev == "-i":
                self.index_col = a
                prev = ""
            elif prev == "-p":
                self.plot_type = a
                prev = ""
            elif a in ["-i","-p"]:
                prev = a
            if a == "-l":
                self.log = np.log(2)
            elif a == "-l10":
                self.log = np.log(10)
        if self.infile == None:
            self.infile = a
        if self.plot_type == None:
            self.plot_type = "box"
        return (self.infile and self.outfile)


    def run(self):
        df = pd.read_table(self.infile, index_col=self.index_col)
        sys.stderr.write("Data file read\n")
        df = df.melt()
        if self.log:
            df.value = df["value"].apply(lambda x: np.log(x+1) / self.log)
        if self.plot_type == "box":
            ax = sns.boxplot(x="variable", y="value", data=df)
        elif self.plot_type == "violin":
            ax = sns.violinplot(x="variable", y="value", data=df)
        elif self.plot_type == "boxen":
            ax = sns.boxenplot(x="variable", y="value", data=df)
        else:
            sys.stderr.write("Invalid plot_type\n")
            self.usage()
        if self.title:
            plt.title(self.title)
        if self.xlabel:
            plt.xlabel(self.xlabel)
        if self.ylabel:
            plt.ylabel(self.ylabel)
        plt.savefig(self.outfile, format=self.imgformat)


    def usage(self):
        sys.stdout.write("""polyplotter.py distr - Draw distribution plot ("box", "violin" or "boxen") for each column in data matrix.

Usage: polyplotter.py distr [options] datafile

Read data from n x m matrix `datafile' and draw a boxplot for each column
onto single plot.

Options related to input data:

    -l    | If supplied, log-scale values (base 2).
    -l10  | If supplied, log-scale values (base 10).

Graphical options:

    -xs S | Set X dimension of image to S inches (default: {}).
    -ys S | Set Y dimension of image to S inches (default: {}).
    -xl L | Set X axis label to L.
    -yl L | Set Y axis label to L.
    -i N  | Set index column number to N. (default: {})
    -p S  | Set plot type to S. (default: {})
    -f F  | Set output image format to F (default: {}).

""".format(DistPlot.xsize, DistPlot.ysize, DistPlot.index_col, DistPlot.plot_type, DistPlot.imgformat))
        sys.exit(1)


class SwarmPlot(Plot):
    def parseArgs():
        pass
    def run():
        pass
    def usage():
        pass

   
class JointPlot(Plot):
    pass

class HeatmapPlot(Plot):
    infile = None
    log = False
    plot_type = None
    skipRows = None
    pointSize = 50
    index_col = 0
    cx = 4
    cy = 3
    
    
    def parseArgs(self, args):
        pass
    

    def run(self):
        df = pd.read_table(self.infile, index_col=self.index_col)
        sys.stderr.write("Data file read\n")
        df = df.melt()
        if self.log:
            df.value = df["value"].apply(lambda x: np.log(x+1) / self.log)
        if self.cluster:
            ax = sns.clustermap(x="variable", y="value", data=df)
        else:
            ax = sns.heatmap(x="variable", y="value", data=df)
        if self.title:
            plt.title(self.title)
        if self.xlabel:
            plt.xlabel(self.xlabel)
        if self.ylabel:
            plt.ylabel(self.ylabel)
        plt.savefig(self.outfile, format=self.imgformat)


    def usage(self):
        sys.stdout.write("""polyplotter.py heat - Draw heatmap

Usage: polyplotter.py heat [options] datafile

Read data from n x m matrix `datafile' and draw a heat map.

Options related to input data:

    -l    | If supplied, log-scale values (base 2).
    -l10  | If supplied, log-scale values (base 10).

Graphical options:

    -xs S | Set X dimension of image to S inches (default: {}).
    -ys S | Set Y dimension of image to S inches (default: {}).
    -xl L | Set X axis label to L.
    -yl L | Set Y axis label to L.
    -i N  | Set index column number to N. (default: {})
    -p S  | Set plot type to S. (default: {})
    -f F  | Set output image format to F (default: {}).

""".format(HeatmapPlot.xsize, HeatmapPlot.ysize, HeatmapPlot.index_col, HeatmapPlot.plot_type, HeatmapPlot.imgformat))
        sys.exit(1)


class ScatterPlot(Plot):
    pass

class KDEPlot(Plot):
    pass

class LinearModelPlot(Plot):
    pass

class LinePlot(Plot):
    pass

class BarPlot(Plot):
    pass

def mainUsage():
    sys.stdout.write("""polyplotter.py - command-line tool to generate a variety of useful plots
    
    Usage: polyplotter.py [command] [command-specific_arguments]

Available subcommands:
    
    dscatt
    mhist
    dist
    heat  -- coming soon

    Ex. polyplotter.py mhist -h

""")

## Main    
if __name__ == "__main__":
    try:
        cmd = sys.argv[1]
        args = sys.argv[2:]
    except IndexError:
        sys.stderr.write("Too few parameters\n")
        mainUsage()
        sys.exit(1)
    if cmd == "dscatt":
        P = DensityPlot()
    elif cmd == "mhist":
        P = MethylHist()
    elif cmd == "dist":
        P = DistPlot()
    else:
        mainUsage()
        sys.exit(1)
    try:    
        if P.parseArgs(args):
            P.run()
        else:
            P.usage()
    except NameError as error:
        sys.stderr.write(error)
        mainUsage()
        sys.exit(1)
