#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat 1st August 2020
Python 2.7 Anaconda Distribution recommended
Written using Spyder and/or Sublime Text

@author:    Bruce J MacLachlan
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: bruce.maclachlan@gmail.com
            
            Artwork of 6R0E for ALGW71
            
"""

description = "A script which uses PyMOL to generate an image of the complex structure of F11 TCR bound to HLA-DR1-PKY. The PyMOL image is raytraced and placed on\
a canvas generated in matplotlib. Biopython is then used to parse info from the structure deposition to annotate the image upon a plaque."

#########################################################################################
###########################          IMPORTS          ###################################
#########################################################################################

import os
import sys
import time
import subprocess
import pymol
import argparse
import itertools

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

######## SOME DEFAULT SETTINGS AT THE TOP FOR EASY ACCESS ########

colourSet =     {   'top'      :   '#acc4ce',
                    'bottom'   :   '#306b64',

                    'HLAA'     :   '#bfbfbf',
                    'HLAB'     :   '#f2f2f2',
                    'peptide'  :   '#305A78',
                    'TCRA'     :   '#1e454d',
                    'TCRB'     :   '#306b64'
                 }

chains_dict =   {   "HLAA"     :    "AAA",
                    "HLAB"     :    "BBB",
                    "peptide"  :    "CCC",
                    "TCRA"     :    "DDD",
                    "TCRB"     :    "EEE"
                }      

#########################################################################################
###########################      ARGUMENT PARSER      ###################################
#########################################################################################

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ray', dest = 'do_ray',   action='store_true', required = False, help ='Do you want to render and save images? Add flag to generate images (time consuming).')
    parser.add_argument('--no_ray', dest = 'do_ray',   action='store_false', required = False, help ='Do you want to render and save images? Add flag to generate no images (quicker).')
    parser.add_argument('--view', dest = 'do_view',   action='store_true', required = False, help ='Do you want to view the pymol session on end? Add flag to leave open at end.')

    parser.add_argument('--palette', dest = 'palette', action='store', required = False, help='Supply a colour palette of the palettable library. example input is palettable.wesanderson.Darjeeling3_5. See https://jiffyclub.github.io/palettable/')

    parser.set_defaults(do_ray=False)
    parser.set_defaults(do_view=False)
    parser.set_defaults(palette=None)
    
    args = parser.parse_args()
    return args

#########################################################################################
###############################      FUNCTIONS      #####################################
#########################################################################################

######## These are my standard PyMOL functions  ############
def initialisePymol():
    '''
    Asks python to start a new pymol session and apply a set of parameters related pymol renders the molecules.
    i.e. I don't like shadows, so they are turned off.
    This helps to keep all figures consistent.
    '''
    print "\nInitialising pymol...\n"
    pymol.finish_launching(['pymol', '-c'])
    pymol.cmd.reinitialize()
    # set PyMOL parameters
    # these are my favourite parameters for making figures
    pymol.cmd.set("ray_shadows","0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("valence", 0)
    pymol.cmd.set("ray_opaque_background", "0")
    pymol.cmd.set("ray_trace_mode", 1)
    pymol.cmd.set("transparency_mode", 2)
    pymol.cmd.set("dash_round_ends", 0)
    return None

def wait4ray(query):  
    '''
    Helps with an issue in rendering PyMOL images where the script proceeds without PyMOL finishing the render job.
    It checks every second whether the image has been created and then allows to proceed.
    '''
    counter = 0
    spinner = itertools.cycle(['-', '/', '|', '\\'])
    while not os.path.exists(query):
        toWrite=spinner.next() + " Time elapsed: "+str(counter)+" seconds"
        sys.stdout.write(toWrite)  # write the next character
        sys.stdout.flush()                # flush stdout buffer (actual character display)
        sys.stdout.write(len(toWrite)*"\b")           # erase the last written chars
        time.sleep(1)
        counter += 1
    return None

def rayTime(saveas, do_ray):
    '''
    Helps with an issue in rendering PyMOL images where the script proceeds without PyMOL finishing the render job.
    It removes any existing file of the same file name, pauses, begins the ray trace and uses wait4ray to check whether the job has been finished.
    '''
    if do_ray == 0:
        return None
    else:
        print "Outputting image.. This may take a few seconds.."
        if os.path.exists(saveas):
            print "Removing "+saveas+" as it already exists!"
            os.remove(saveas)
        time.sleep(5)
        pymol.cmd.png(saveas,ray=do_ray,width="20cm",height="30cm", dpi=300)
        wait4ray(saveas) 
        print "Done! "+str(saveas)+ " was outputted" 
        return None

#### Some colour functions #####
    
def rgb255_to_fraction(rgb255list):  
    '''
    pymol deals with colours as rgb values in fractions from 0-1.0
    Use to convert rgb values in 255 format to fraction format.
    i.e.
    
    [255,255,255] to [1.0, 1.0, 1.0] e.g. white
    '''
    output = []
    for v in rgb255list:
        output.append(v/255.0)
    return output


def hex_to_fraction(hex_str):
    '''
    pymol deals with colours as rgb values in fractions from 0-1.0
    Use to convert hex values  to fraction format.
    i.e.
    
    "#FFFFFF" to [1.0, 1.0, 1.0] e.g. white
    '''
    if hex_str.startswith('#'):
        hex_str = hex_str[1:]
    rgb255list = tuple([int(hex_str[i:i + 2], 16) for i in xrange(0, len(hex_str), 2)])
    output = []
    for v in rgb255list:
        output.append(v/255.0)
    print output
    return output


def set_new_colour(name, rgblist):
    '''
    Adds a new colour to the pymol colourset list.
    Bit obsolete because it ended up being a one-liner.
    '''
    pymol.cmd.set_color(name, rgblist, mode=0, quiet=1)
    return None


def show_palette(palette):
    '''
    Makes a small plot in matplotlib and outputs png file of the colour palette requested from palettable. It pops out in the root directory.
    Select your colours from bottom [0] to top [-1]
    '''

    print "Selecting from palettable palette:"
    print palette.name
    print "Number of colours in palette:"
    print palette.number
    fig = plt.figure(frameon=True)

    ax = fig.add_axes([0,0, palette.number, 1])
    ax.axis('on')

    i = 0
    while i < palette.number:
        ax.axhspan(i, i+1.0, facecolor=palette.mpl_colors[i], zorder=1)
        ax.annotate(str(i), xy=(0.5, float(i)+0.5), zorder=10)
        i+=1

    print "Saving a colour swatch of colours in palette. See /"+palette.name+".png .."
    fig.savefig(palette.name+".png")
    

    return None

def load_palette(colourSet, palette):
    '''
    When supplied with the --palette flag and a palettable palette has been supplied, the user is asked to give a colour for each object of the image.
    This is pretty clunky but easy way of exploring colour arrangements that you like.
    Requires raw unput from command line.
    Outputs a dictionary of the selected colours in rgb 0.0 to 1.0
    '''

    import importlib

    mod = palette.split(".")[-1]
    package = ".".join(palette.split(".")[0:-1])

    palette_package = importlib.import_module(package)
    colour_object = getattr(palette_package, mod)

    show_palette(colour_object)

    output = {}
    
    standard_colours = {"white" : (1.0, 1.0, 1.0),
                        "lightgrey" : (0.75, 0.75, 0.75),
                        "grey": (0.5, 0.5, 0.5),
                        "darkgrey": (0.25, 0.25, 0.25),
                        "black": (0.0, 0.0, 0.0)
    }

    print "Pick a colour from the colour swatch. Use 0 (bottom colour) 1,2,3.. etc. Else choose from white, lightgrey, grey, darkgrey, black."

    for item in colourSet.keys():
        selected_colour = raw_input('Select colour for '+item+' :')
        if selected_colour.isdigit():
            output[item] = colour_object.mpl_colors[int(selected_colour)]
        elif selected_colour in standard_colours.keys():
            output[item] = standard_colours[selected_colour]
    return output

################################### Figure generator ####################################

def structure_layer(structure, do_ray):
    '''
    This is where all the PyMOL work is done.
    An image on a transparent canvas is outputted as png.


    Uses PyMOL to visualise the structure provided.
    make do_ray = 1 to ray trace the image (time consuming)
    make do_ray = 0 to skip retracing i.e. debugging or testing other parts of the code.
    '''
    # Building our first figure
    pymol.cmd.hide("everything", "all")

    # show HLA
    pymol.cmd.show("surface", "HLAA_obj")
    pymol.cmd.show("surface", "HLAB_obj")

    # give these surfaces some transparency
    pymol.cmd.set("transparency", 0.1, "HLAA_obj")
    pymol.cmd.set("transparency", 0.1, "HLAB_obj")

    # show HLA helices
    pymol.cmd.show("cartoon", "HLA_a1a2")

    # show peptide
    pymol.cmd.show("surface", "peptide_obj")
    pymol.cmd.set("transparency", 0.1, "peptide_obj")

    # show peptide as sticks
    pymol.cmd.show("sticks", "peptide_obj")
  
    # show TCR
    pymol.cmd.show("surface", "TCRA_obj")
    pymol.cmd.set("transparency", 0.0, "TCRA_obj")
    
    pymol.cmd.show("surface", "TCRB_obj")
    pymol.cmd.set("transparency", 0.0, "TCRB_obj")

    # set the view, save the scene and render the image (if do_ray=True)
    complex_view = ("\
    -0.194591478,    0.968047142,   -0.158180326,\
    -0.116970502,   -0.183013782,   -0.976128042,\
    -0.973888516,   -0.171443373,    0.148845106,\
     0.000085423,    0.000001311, -432.039001465,\
   -35.504920959,   24.705741882,  -28.177782059,\
   377.880249023,  486.194641113,   20.000000000")
    
    pymol.cmd.set_view(complex_view)
    pymol.cmd.scene("complex_image", "store")
    rayTime("complex_image.png", do_ray)

    return None

#########################################################################################
################################### Plaque generator ####################################


def author_on_lines(author_list, max_on_line):
    '''
    This is some really clunky code to try and wrap the author list onto a sensibile width.
    It sets the maximum number of author entries to max_on_line and generates a string with newlines.
    '''

    author_length = len(author_list)
    author_lines =  len(author_list) / max_on_line
    author_remainder = len(author_list) % max_on_line

    length_to_split = [max_on_line]*author_lines
    length_to_split.append(author_remainder)

    from itertools import islice

    # Using islice 
    Inputt = iter(author_list) 
    author_out = [list(islice(Inputt, elem)) 
              for elem in length_to_split] 

    output = ''
    for line in author_out:
        output += ', '.join(line)
        output += "\n"

    return output[:-1]

def title_on_lines(title, max_on_line):
    '''
    This is some really clunky code to try and wrap the title onto a sensibile width.
    It sets the maximum number of words in the title to max_on_line and generates a string with newlines.
    '''

    title_list = title.split(" ")
    title_length = len(title_list)
    title_lines =  len(title_list) / max_on_line
    title_remainder = len(title_list) % max_on_line

    length_to_split = [max_on_line]*title_lines
    length_to_split.append(title_remainder)

    from itertools import islice

    # Using islice 
    Inputt = iter(title_list) 
    title_out = [list(islice(Inputt, elem)) 
              for elem in length_to_split] 
    output = []
    for line in title_out:
        new_line = ''
        new_line = ' '.join(line)

        output.append(new_line)

    return output

def parse_pdb_info(pdb_id):
    '''
    Rips out some info from the mmcif header and returns a dict
    Used to form the plaque.
    This is currently only mmcif compatible and not pdb (because 6R0E is uploaded only in cif format.. no PDB)
    '''

    # parser = MMCIFparser()
    # structure = parser.get_structure('complex', '6r0e.cif')
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    mmcif_dict = MMCIF2Dict(pdb_id)

    output = {}

    # select out the values from the mmcif file which have some difficult to read names.
    # mmcif2dict parses a list from each entry. Therefore we need to index the first entry hence the [0]

    output["id"] = mmcif_dict["_entry.id"][0]
    output["title"] = mmcif_dict["_citation.title"][0]
    output["journal"] =  mmcif_dict["_citation.journal_abbrev"][0]
    output["year"] = mmcif_dict["_citation.year"][0]
    output["authors"] = mmcif_dict["_citation_author.name"]
    output["resolution"] = float(mmcif_dict["_reflns.d_resolution_high"][0])
    output["spacegroup"] = mmcif_dict["_symmetry.space_group_name_H-M"][0]

    output["a"] = float(mmcif_dict["_cell.length_a"][0])
    output["b"] = float(mmcif_dict["_cell.length_b"][0])
    output["c"] = float(mmcif_dict["_cell.length_c"][0])

    output["alpha"] = float(mmcif_dict["_cell.angle_alpha"][0])
    output["beta"] = float(mmcif_dict["_cell.angle_beta"][0])
    output["gamma"] = float(mmcif_dict["_cell.angle_gamma"][0])

    # Let's stick all the unit cell info into a new dict entry which is a string as will be printed on the plaque
    output["unitcell"] =  r'$a$ = ' + '%.1f' % output["a"] + ", " + r'$b$ = ' + '%.1f' % output["b"] + ", " + r'$c$ = ' + '%.1f' % output["c"] + \
                          ", " + r'$\alpha$ = ' + '%.1f' % output["alpha"]+ r'$\degree$'+", " + r'$\beta$ = ' + '%.1f' % output["beta"] + r'$\degree$'+ ", " + r'$\gamma$ = ' + '%.1f' % output["gamma"]+ r'$\degree$'

    # wrap the authors and title into a set line length
    output["authorlines"] = author_on_lines(mmcif_dict["_citation_author.name"], 8)
    output["titlelines"] = title_on_lines(mmcif_dict["_citation.title"][0], 9)

    return output

#########################################################################################
################################## Canvas generator #####################################
def canvas(image_name, info_dict):
    '''
    Creates the canvas for the artwork.
    This takes the form of a matplotlib plot with no axes and no space around the edges.
    We can then treat the plot like a canvas
    '''

    # ferameon=False required to remove the margin around the plot
    fig = plt.figure(figsize=(11.7,16.5), dpi = 300, facecolor="w", frameon=False)

    # puts the canvas onto a axis scale between 0 and 1 in both directions
    ax = fig.add_axes([0, 0, 1, 1])
    # turn off axes
    ax.axis('off')
    # This is reqiured to remove margin around the plot
    plt.autoscale(tight=True)

    # This creates a split background. The bottom is meant to represent the cell surface.
    ax.axhspan(0, 0.33,facecolor=colourSet['bottom'])
    ax.axhspan(0.33, 1.0 ,facecolor=colourSet['top'])

    # load the png image outputted by PyMOL
    im = plt.imread(image_name)

    # show the image on the canvas      
    fig.figimage(im, 600, 1000, zorder=10)

    # Generate the "plaque"
    # Forms a big string to be annotated onto the canvas
    big_string = r'$\bf{%s}$' % str(info_dict["id"])
    big_string += "\n" + "Resoution: "+'%.1f' % info_dict["resolution"] + r' $\AA$' + ". " +"Space group: "+info_dict["spacegroup"]+ ". "+ "Unit cell: "+info_dict["unitcell"]+"."+"\n"
    
    for line in info_dict["titlelines"]:
        line_new = line.replace(" ", "\\ ")
        big_string += "\n" + r'$\bf{%s}$' % line_new

    big_string += "\n" + info_dict["journal"]+" "+info_dict["year"]
    big_string += "\n\n" + info_dict["authorlines"]

    # Place the plaque
    ax.annotate(big_string, xy=(0.5, 0.075), xycoords='axes fraction', ha="center", bbox=dict(facecolor='white', edgecolor='black', pad=10.0, linewidth=5.0))

    # Add a signature link to this source code on github
    watermark = ax.annotate("Source code: github.com/brucemaclachlan/6R0E_artwork", xy=(0.95, 0.035), xycoords='axes fraction', ha="right", color='black', weight='bold', fontsize=9)
    watermark.set_alpha(0.2)

    import time
    timestr = time.strftime("%Y%m%d-%H%M%S")


    fig.savefig("canvas_"+timestr+".png", transparent=True, pad_inches = 0)
    plt.savefig("canvas_"+timestr+".pdf", pad_inches = 0)

    return None


#########################################################################################
############################       END FUNCTIONS      ###################################
#########################################################################################

#########################################################################################
################################       BODY       #######################################
#########################################################################################
# This is the "start" of the script.
    
# Handle commoand line arguments
args = parse_args()
print "Running analysis with the following inputs.. "
print args
ray = args.do_ray
view = args.do_view
palette = args.palette

# need to add this line as pymol.png wants 0 or 1 not bool
if ray == True:
    do_ray = 1
if ray == False:
    do_ray = 0

# If no palette is requested, use the default colours at the start of the script.
if palette != None:
    colourSet = load_palette(colourSet, palette)
else:
    None

# This creates a new session of PyMOL and runs my favourite viewing parameters for PyMOL
# These parameters can be seen in def intialisePymol() function above
initialisePymol()

# Fetch the mmcif file of the structure from RCSB
pymol.cmd.fetch('6r0e', name="complex", type='cif')

# sort the chains into objects
chains = chains_dict.keys()

print colourSet

# load the different chains, colour them and create objects out of each chain
for chain in chains:
    pymol.cmd.select(chain, selection="complex and chain "+chains_dict[chain])

    if colourSet[chain][0] == "#":
        set_new_colour(chain+"_colour", hex_to_fraction(colourSet[chain]))
    else:
        set_new_colour(chain+"_colour", colourSet[chain])
    pymol.cmd.color(chain+"_colour", chain)
    pymol.cmd.create(chain+"_obj", selection=chain)

pymol.cmd.select("HLA_a1a2", selection="complex and chain "+chains_dict["HLAA"]+" and resi 46-78 or complex and chain "+chains_dict["HLAB"]+" and resi 54-91")
pymol.cmd.create("HLA_a1a2_obj", selection="HLA_a1a2")

# generate the image of the structure using PyMOL
structure_layer("complex", do_ray)

# Parse the info from the mmcif file to create the info for the plaque
structure_info = parse_pdb_info('6r0e.cif')

# Generate the canvas, place the PyMOL image onto the canvas and draw the background and plaque
canvas("complex_image.png", structure_info)

###########
# This is the end. We finally want to save the PyMOL session file and quit pymol. 
# If we have asked to view the results in PyMOL, we will finish by opening up the session file.
###########
pymol.cmd.save("6R0E_artwork.pse")
pymol.cmd.quit()
if view == True:
    subprocess.call(["pymol", "6R0E_artwork.pse"])

#########################################################################################
#################################     END     ###########################################
#########################################################################################
