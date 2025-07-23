#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-07-22 at 10:18:17 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to plot the transmission eigenvalue distribution.
import sys, os, datetime
import compile_tikz

def plot_distrib(args):
    """
    Plot the distribution given by the first argument: filename=args[1].
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 2):
        print(compile_tikz.TAG_ERROR + "No input file, doing nothing...")
        print(compile_tikz.TAG_USAGE + args[0] + " DISTRIB_FILE")
        return 1
    
    distrib_file = args[1]
    file_path = os.path.splitext(distrib_file)[0]  ## The file path is the filename without its extension (used to write new files). 
    
    tikz_code = """%% Generated on {timestamp} by {my_program} {my_copyright}
\\begin{{tikzpicture}}%
\\pgfmathsetmacro\\tavg{{0.239}}%% Average transmission, Tavg ~ 1/[1 + (2/pi) * (L/lscat)]
\\begin{{axis}}[%
    title={{\\detokenize{{{distrib_file}}}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    xmin={xmin}, xmax={xmax},
    ymin=0.02, ymax=200,
    ymode=log,
    scatter/classes={{%%
        success={{black, mark=*, mark size=0.7}},
        failure={{red, mark=triangle*, mark size=1.7}}%% No comma.
    }},
    scatter src=explicit symbolic,
    axis on top=true,
    unbounded coords=jump,  %% Discard NaN's and negative entries.
    clip marker paths=true, %% Clips the marks out of the axis frame.
    clip mode=individual,   %% Ensure the marks do not overlay the other curves.
    table/col sep=comma,
]%
\\addplot[black!30, smooth, domain=0.001:0.999, samples=64] ({{sin(90*\\x)^2}}, {{\\tavg/(2 * sin(90*\\x)^2 * cos(90*\\x))}}); %% Plot the bimodal distribution.
\\addplot[scatter, black] table[x=tval,y=rho,meta=converged]{{\\jobname.csv}};
\\end{{axis}}%
\\end{{tikzpicture}}%""".format(
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = compile_tikz.MY_COPYRIGHT,
        distrib_file = distrib_file,
        xlabel = "Transmission eigenvalue $T$",
        ylabel = "Distribution $\\rho(T)$",
        xmin   = 0,
        xmax   = 1
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = file_path + '.tikz'
    print(compile_tikz.TAG_INFO + "Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_code)
    compile_tikz.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    return 0

if (__name__ == '__main__'):
    exit(plot_distrib(sys.argv))
