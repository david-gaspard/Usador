#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2025-07-17 at 12:48:24 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Python script to test the plotting of a scalar field defined over a square lattice with holes. The input data must have the form [x, y, f(x, z)].
import sys, os, datetime, csv
import numpy as np
import matplotlib.pyplot as mplt
import matplotlib.colors as mcol
import compile_tikz

## Create the 'sunset' colormap originally created by David Gaspard in June 2024:
MY_COPYRIGHT = "(c) 2025 David GASPARD (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr>"
SUNSET_COLORS = [[1.00000, 1.00000, 1.00000],
                 [1.00000, 0.99337, 0.65453],
                 [1.00000, 0.96893, 0.46622],
                 [1.00000, 0.93848, 0.36430],
                 [1.00000, 0.90595, 0.29694],
                 [1.00000, 0.87228, 0.24748],
                 [1.00000, 0.83776, 0.20892],
                 [0.99999, 0.80247, 0.17768],
                 [0.99999, 0.76641, 0.15171],
                 [0.99998, 0.72953, 0.12971],
                 [0.99996, 0.69173, 0.11084],
                 [0.99994, 0.65285, 0.09448],
                 [0.99991, 0.61270, 0.08024],
                 [0.99985, 0.57099, 0.06781],
                 [0.99977, 0.52734, 0.05700],
                 [0.99964, 0.48119, 0.04774],
                 [0.99942, 0.43165, 0.04006],
                 [0.99904, 0.37733, 0.03425],
                 [0.99829, 0.31562, 0.03121],
                 [0.99653, 0.24108, 0.03388],
                 [0.98995, 0.14079, 0.05718],
                 [0.95565, 0.04309, 0.16091],
                 [0.89867, 0.01496, 0.28248],
                 [0.83663, 0.00638, 0.38039],
                 [0.77159, 0.00295, 0.45924],
                 [0.70427, 0.00141, 0.52191],
                 [0.63527, 0.00068, 0.56978],
                 [0.56522, 0.00032, 0.60344],
                 [0.49490, 0.00015, 0.62312],
                 [0.42519, 0.00006, 0.62886],
                 [0.35720, 0.00003, 0.62072],
                 [0.29217, 0.00001, 0.59895],
                 [0.23145, 0.00000, 0.56414],
                 [0.17636, 0.00000, 0.51734],
                 [0.12806, 0.00000, 0.46007],
                 [0.08734, 0.00000, 0.39419],
                 [0.05461, 0.00000, 0.32169],
                 [0.02990, 0.00000, 0.24450],
                 [0.01290, 0.00000, 0.16431],
                 [0.00313, 0.00000, 0.08248],
                 [0.00000, 0.00000, 0.00000]]
SUNSET_NSAMPLE = len(SUNSET_COLORS)
SUNSET_NODES = np.linspace(0., 1., SUNSET_NSAMPLE)
SUNSET_CMAP = mcol.LinearSegmentedColormap.from_list("sunset_cmap", list(zip(SUNSET_NODES, SUNSET_COLORS)))
SUNSET_CMAP.set_bad('white', 0.) ## Set the color when nan is encountered. Args: (color, opacity).

def colormap_to_tikz_code(cmap, nsample):
    """
    Returns a TikZ code version of the given colormap "cmap" using a given number of samples "nsample".
    Example output: 'colormap={temperature}{rgb255=(0,0,128) rgb255=(0,0,255) rgb255=(255,255,255) rgb255=(255,0,0) rgb255=(128,0,0)}'
    """
    ##print("[INFO] cmap.name =", cmap.name, ", cmap.N =", cmap.N, ", cmap(0.5) =", cmap(0.5))
    
    string = "colormap={" + cmap.name + "}{"
    
    for i in range(nsample):
        val = i/(nsample - 1)
        rgb255 = [int(np.round(255*c)) for c in cmap(val)[0:3]]
        string += "rgb255=(" + str(rgb255[0]) + ","+ str(rgb255[1]) +","+ str(rgb255[2]) +") "
    
    return string + "}"

def boundary_to_tikz_code(mesh_file):
    """
    Returns a TikZ code version of the given mesh file.
    The mesh file has the format [x, y, north, south, east, west]. The components 'x' and 'y' are assumed to be integers, 
    and the direction components can be either point indices or boundary conditions ('mirror', 'open', 'input', or 'output').
    """
    try:
        fp = open(mesh_file, 'r')
    except IOError as e:
        print("[WARN] Mesh file '" + mesh_file + "' not found, skipping...")
        return "%% Mesh file '" + mesh_file + "' not found..."
    
    ls = csv.DictReader((line for line in fp if not line.startswith('%')), skipinitialspace=True)
    
    string = "\\begin{scope}[mirror/.style={black, thick}, open/.style={opacity=0}, input/.style={red!80, thick}, output/.style={blue!80, thick}]%\n"
    
    for p in ls:
        if (not p['north'].isdigit()):
            x0 = int(p['x']) - 0.5;
            y0 = int(p['y']) + 0.5;
            x1 = int(p['x']) + 0.5;
            y1 = int(p['y']) + 0.5;
            string += "\\draw[{bnd}] (axis cs:{x0}, {y0}) -- (axis cs:{x1}, {y1});\n".format(bnd=p['north'], x0=x0, y0=y0, x1=x1, y1=y1)
        if (not p['south'].isdigit()):
            x0 = int(p['x']) + 0.5;
            y0 = int(p['y']) - 0.5;
            x1 = int(p['x']) - 0.5;
            y1 = int(p['y']) - 0.5;
            string += "\\draw[{bnd}] (axis cs:{x0}, {y0}) -- (axis cs:{x1}, {y1});\n".format(bnd=p['south'], x0=x0, y0=y0, x1=x1, y1=y1)
        if (not p['east'].isdigit()):
            x0 = int(p['x']) + 0.5;
            y0 = int(p['y']) + 0.5;
            x1 = int(p['x']) + 0.5;
            y1 = int(p['y']) - 0.5;
            string += "\\draw[{bnd}] (axis cs:{x0}, {y0}) -- (axis cs:{x1}, {y1});\n".format(bnd=p['east'], x0=x0, y0=y0, x1=x1, y1=y1)
        if (not p['west'].isdigit()):
            x0 = int(p['x']) - 0.5;
            y0 = int(p['y']) - 0.5;
            x1 = int(p['x']) - 0.5;
            y1 = int(p['y']) + 0.5;
            string += "\\draw[{bnd}] (axis cs:{x0}, {y0}) -- (axis cs:{x1}, {y1});\n".format(bnd=p['west'], x0=x0, y0=y0, x1=x1, y1=y1)
    
    fp.close()
    
    return string + "\\end{scope}"

def plot_map(args):
    """
    Plots the given component args[1]='column_name' of the given CSV file args[2] + "_field.csv".
    The data in the field file are assumed to have the form [x, y, f1, f2, ..., fn], and the components 'x' and 'y' are assumed to be integers.
    It also incorporates the boundary conditions using the mesh file args[2] + "_mesh.csv".
    The data in the mesh file are assumed to have the form [x, y, north, south, east, west]. See also: boundary_to_tikz_code().
    """
    ## Check if the number of arguments is correct:
    if (len(args) != 3):
        print("[USAGE] " + args[0] + " COLUMN_NAME FIELD_PATH")
        return 1
    
    column_name = args[1]  ## Interpret arg #1 as the name of the column in the field file.
    field_path  = args[2]  ## Interpret arg #1 as the name of the field file.
    field_file = field_path + ".csv"
    mesh_file = field_path + "_mesh.csv"
    
    try:
        fp = open(field_file, 'r')
    except IOError as e:
        print("Field file '" + field_file + " not found, aborting now...")
        return 1
    
    ls = list(csv.reader(line for line in fp if not line.startswith('%')))
    fp.close()
    
    columns = [elem.strip() for elem in ls[0]]
    idx = columns.index(column_name)
    ##print("[INFO] columns =", columns, ", idx =", idx)
    
    point = np.asarray([line[0:2] for line in ls[1:]], dtype=int)
    field = np.asarray([line[idx] for line in ls[1:]], dtype=float)
    npt = point.shape[0]
    ##print("[INFO] npoint =", npt)
    
    xmin = point[:, 0].min()
    xmax = point[:, 0].max()
    ymin = point[:, 1].min()
    ymax = point[:, 1].max()
    
    ##print("[INFO] xmin =", xmin, ", xmax =", xmax, ", ymin =", ymin, ", ymax =", ymax)
    ##print("[INFO] field =", field[0:5])
    
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1
    matrix = np.full((ny, nx), np.nan)
    
    ## Loop on the points to write in 'matrix':
    for ip in range(npt):
        i = ymax - point[ip, 1]
        j = point[ip, 0] - xmin
        matrix[i, j] = field[ip]
    
    matrix = np.ma.array(matrix, mask=np.isnan(matrix))  ## Use a mask to escape nan values.
    
    cmap = SUNSET_CMAP  ## Use custom 'sunset' colormap.
    ##cmap = mplt.cm.jet ## Use 'jet' colormap.
    ##mplt.imshow(matrix, cmap=cmap)  ## Show the plot in live (optional).
    ##mplt.colorbar()
    ##mplt.show()
    
    vmin = matrix.min() ## Extract the depth range of the field [vmin, vmax].
    vmax = matrix.max()
    norm = mplt.Normalize(vmin=vmin, vmax=vmax)
    image = cmap(norm(matrix)) ## Create the bitmap image.
    ##print("[INFO] vmin = ", vmin, ", vmax = ", vmax)
    
    ##pre, ext = os.path.splitext(field_file)
    bitmap_file = field_path + '_map.png'
    mplt.imsave(bitmap_file, image)  ## Save the raw pixel-constrained bitmap to a PNG file.
    
    tikz_string = '''%% Generated on {timestamp} by {my_program} {my_copyright}
\\begin{{tikzpicture}}%
\\begin{{axis}}[%
    title={{\\detokenize{{{title}}}}},
    xlabel={{{xlabel}}},
    ylabel={{{ylabel}}},
    colorbar, %% Enable colobar.
    point meta min={vmin}, %% Set colorbar range.
    point meta max={vmax},
    {colorbar_string},
    enlargelimits=true, %% Allow for larger view.
    axis equal image, %% Unit aspect ratio.
    axis on top=true, %% Prevent the image from hiding the axes.
]%
\\addplot graphics[xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax}]{{{bitmap_file}}};
{boundary_code}
\\end{{axis}}%
\end{{tikzpicture}}%'''.format(#
        timestamp = datetime.datetime.now().astimezone().strftime("%F at %T %z"),
        my_program = args[0],
        my_copyright = MY_COPYRIGHT,
        title  = column_name,
        xlabel = "$x$",
        ylabel = "$y$",
        colorbar_string = colormap_to_tikz_code(cmap, 41),
        vmin   = vmin,
        vmax   = vmax,
        xmin   = xmin-0.5,
        xmax   = xmax+0.5,
        ymin   = ymin-0.5,
        ymax   = ymax+0.5,
        bitmap_file = "\\jobname_map.png",
        boundary_code = boundary_to_tikz_code(mesh_file)
    )
    
    ## Export the TikZ code to a file and compile it:
    tikz_file = field_path + '.tikz'
    print("[INFO] Writing TikZ file: '" + tikz_file + "'...")
    open(tikz_file, 'w').write(tikz_string)
    compile_tikz.compile_tikz(tikz_file) ## Compile the TikZ file.
    
    return 0

if (__name__ == '__main__'):
    exit(plot_map(sys.argv))
