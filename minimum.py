""" Minimum working example of an SME script 
"""
import os.path

from pysme.gui import plot_plotly
from pysme import sme as SME
from pysme import util
from pysme.solve import solve
from pysme.synthesize import synthesize_spectrum

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.persistence import save_as_idl

if __name__ == "__main__":

    # Define the location of all your files
    # this will put everything into the example dir
    target = "sun"
    examples_dir = os.path.dirname(os.path.realpath(__file__))
    in_file = os.path.join(examples_dir, "sun_6440_grid.inp")
    out_file = os.path.join(examples_dir, f"{target}.sme")
    plot_file = os.path.join(examples_dir, f"{target}.html")
    log_file = os.path.join(examples_dir, f"{target}.log")

    # Start the logging to the file
    util.start_logging(log_file)

    # Load your existing SME structure or create your own
    sme = SME.SME_Structure.load(in_file)
    sme.abund = Abund.solar()
    # sme.linelist = ValdFile(os.path.join(examples_dir, "sun.lin"))
    sme.linelist = ValdFile('sun.lin.5')
    # Change parameters if your want
    sme.vsini = 0
    sme.vrad = 0.35
    sme.teff, sme.logg, sme.monh = 5000, 4.4, -0.1
    sme.h2broad = True
    sme.vrad_flag = "whole"
    sme.cscale_flag = "linear"
    sme.cscale_type = "mask"


    sme.atmo.source = "marcs2012.sav"
    sme.atmo.method = "grid"
    sme.atmo.geom = "PP"
    
    
    sme.nlte.set_nlte("Fe", "marcs2012_Fe2016.grd")
    sme.nlte.set_nlte("Si", "marcs2012_Si2016.grd")
    sme.nlte.set_nlte("Ca", 'marcs2012p_t1.0_Ca.grd')    
    
    # Define any fitparameters you want
    # For abundances use: 'abund {El}', where El is the element (e.g. 'abund Fe')
    # For linelist use: 'linelist {Nr} {p}', where Nr is the number in the
    # linelist and p is the line parameter (e.g. 'linelist 17 gflog')
    fitparameters = ['vmic','vmac']

    # Start SME solver
    # sme = synthesize_spectrum(sme)
    sme = solve(sme, fitparameters)

    print(sme.citation())

    # Save results
    sme.save(out_file)

    # Plot results
    fig = plot_plotly.FinalPlot(sme)
    fig.save(filename=plot_file)
