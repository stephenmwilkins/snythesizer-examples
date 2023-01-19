

import numpy as np
import matplotlib.pyplot as plt
from unyt import yr, Myr

from synthesizer.grid import get_available_lines, Grid
from synthesizer.line import get_line_id
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh
from synthesizer.galaxy.parametric import ParametricGalaxy as Galaxy


if __name__ == '__main__':

    grid_dir = '/Users/stephenwilkins/Dropbox/Research/data/synthesizer/grids'

    default_grid_name = 'fsps-v3.2_imf3:2.3_cloudy-v17.03_log10Uref-2'

    available_lines, wavelengths = get_available_lines(
        default_grid_name, grid_dir, include_wavelengths=True)

    print(available_lines)
    print(wavelengths)

    """ Predict the equivalent widths of various lines for different IMFs """

    # list of lines. Lines in nested lists (or tuples) ot strings containing commas denote doublets for which the combined line properties are calculated
    line_ids = ['HI4861', 'OIII4959', 'OIII5007', 'OIII4959,OIII5007']
    line_ids = available_lines  #  all lines

    # define the parameters of the star formation and metal enrichment histories
    sfh_p = {'duration': 10 * Myr}
    Z_p = {'log10Z': -3.0}  # can also use linear metallicity e.g. {'Z': 0.01}

    # define the functional form of the star formation and metal enrichment histories
    sfh = SFH.Constant(sfh_p)  # constant star formation
    Zh = ZH.deltaConstant(Z_p)  # constant metallicity

    # open test grid though without reading spectra BUT reading the required lines

    # this is NOT the best way of doing this. It would be better to make a glob of different grids and then read the high-mass slope from the grid itself. However, those attributes are not currently set.
    a3s = np.arange(1.5, 3.1, 0.1)
    # a3s = [2.3]
    grid_names = [
        f'fsps-v3.2_imf3:{a3:.1f}_cloudy-v17.03_log10Uref-2' for a3 in a3s]

    # dictionary holding the line EWs for different grids (IMFs in this context)
    line_ews = {line_id: np.zeros(len(a3s)) for line_id in line_ids}  #  uses line_id as key

    for i, grid_name in enumerate(grid_names):

        print(grid_name)

        grid = Grid(grid_name, grid_dir=grid_dir, read_spectra=False, read_lines=line_ids)

        #  get the 2D star formation and metal enrichment history for the given SPS grid. This is (age, Z).
        sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

        # create the Galaxy object
        galaxy = Galaxy(sfzh)

        # --- create the Lines dictionary which contains line objects
        lines = galaxy.get_intrinsic_line(grid, line_ids)

        # --- print summaries of each line
        for line_id, line in lines.items():
            line_ews[line_id][i] = line._ew  # the unit-less value
            # line_ews[line_id][i] = line.ew.value  # the value of the unyt quantity
            # line_ews[line_id][i] = line.ew.to('Angstrom').value  # explicit conversion to Angstrom. Actually unnecessary here as the implicit units are Angstrom

    print(line_ews['HI4861'])

    # EW vs high-mass slope plot

    fig = plt.figure(figsize=(3.5, 5))

    left = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    s = (wavelengths < 7000.)

    for line_id in np.array(line_ids)[s]:

        ew_limit = 10

        if line_ews[line_id][0] > ew_limit:
            ax.plot(a3s, line_ews[line_id])
            ax.text(1.45, line_ews[line_id][0],
                    line_id, ha='right', va='center', fontsize=6)

    ax.set_yscale('log')
    ax.set_xlim([1.0, 3.0])
    ax.set_xlabel('high-mass slope')
    ax.set_ylabel(r'$\rm EW/\AA$')

    fig.savefig('figs/line_ew_imf.pdf')

    # wavelength vs. different in EW

    fig = plt.figure(figsize=(5, 3.5))

    left = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    for line_id, wavelength in zip(line_ids, wavelengths):

        ew_limit = 0.1

        if line_ews[line_id][0] > ew_limit:

            diff = np.log10(line_ews[line_id][0]/line_ews[line_id][-1])

            ax.scatter(wavelength, diff)
            ax.text(wavelength, diff,
                    line_id, ha='right', va='center', fontsize=6)

    # ax.set_yscale('log')
    ax.set_xlim([1000, 7000])
    ax.set_xlabel(r'$\rm \lambda/\AA$')
    ax.set_ylabel(r'$\rm log_{10}(EW(\alpha_3=1.5)/EW(\alpha_3=3.0))$')

    fig.savefig('figs/line_ew_imf_diff.pdf')

    # EW vs difference in EW

    fig = plt.figure(figsize=(5, 4.5))

    left = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    for line_id, wavelength in zip(line_ids, wavelengths):

        diff = np.log10(line_ews[line_id][0]/line_ews[line_id][-1])

        if diff > 0.3:
            ax.scatter(line_ews[line_id][0], diff)
            ax.text(line_ews[line_id][0], diff,
                    line_id, ha='right', va='center', fontsize=6)

    ax.set_xscale('log')
    ax.set_xlim([0.1, 3000])
    ax.set_ylim([0.3, 0.6])
    ax.set_xlabel(r'$\rm EW(\alpha_3=1.5)$')
    ax.set_ylabel(r'$\rm log_{10}(EW(\alpha_3=1.5)/EW(\alpha_3=3.0))$')

    fig.savefig('figs/line_ew_imf_diff2.pdf')
