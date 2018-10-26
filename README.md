# polyplot
Command-line tool to generate useful plots

## Basic usage

polyplotter.py - command-line tool to generate a variety of useful plots
    
    Usage: polyplotter.py [command] [command-specific_arguments]

Available subcommands:
    
    dist    - Draw distribution plot ("box", "violin" or "boxen") for each column in data matrix
    dscatt  - Draw density scatterplots of paired data
    mhist   - Draw histogram of methylation data

## Example

polyplotter.py dscatt -h
    
    polyplotter.py dscatt - Draw density scatterplots of paired data.
    
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
    
