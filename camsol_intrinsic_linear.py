#!/usr/bin/env python
#  Copyright (c) 2020.
#  The Chancellor, Masters and Scholars of the University of Cambridge,
#  Dr Pietro Sormanni and Prof Michele Vendruscolo, all rights reserved.
#  This program is distributed with a license agreement that accompanies it
#  ("CamSol License agreement"),
#  by using this software you agree to the terms of this license
#  You should have received a copy of the license with this software.
#  IF YOU HAVEN'T RECEIVED THE LICENSE DO NOT USE THE SOFTWARE,
#  you can ask for a licence by emailing
#  Prof M. Vendruscolo (mv245@cam.ac.uk) or Dr P. Sormanni (ps589@cam.ac.uk)


import os
import sys

import numpy

import pandas as pd

# below is useful if running this file as a script from master folder and camsol is not installed as a package
if os.path.exists("camsol"):
    sys.path.append("camsol")
if os.path.exists("../camsol"):
    sys.path.append("../camsol")
try:
    extra = os.path.join(os.path.split(os.path.realpath(__file__))[0], "camsol")
    if os.path.exists(extra):
        sys.path.append(extra)
    extra2 = os.path.join(os.path.split(os.path.realpath(__file__))[0], "../camsol")
    if os.path.exists(extra):
        sys.path.append(extra2)
except Exception:
    pass


try:
    from . import (
        plotter,  # plotting module, not strictly required. Try import as using camsol as package/module
    )
except Exception:
    try:
        from camsol import (
            plotter,  # running camsol_intrinsic as script, (but camsol module installed)
        )
    except Exception:
        try:
            import plotter  # see if plotter module is in pythonpath
        except Exception:

            class tmpplotter:  # dummy class, raise exception only if asked to plot.
                default_figure_sizes = {}

                def profile(self, *args, **kwargs):
                    raise Exception("PLOTTER MODULE NOT AVAILABLE\n")

                def plot_seq_profile(self, *args, **kwargs):
                    raise Exception("PLOTTER MODULE NOT AVAILABLE\n")

            plotter = tmpplotter()


def smooth_profile(profile, smooth_per_side=3, weights=None):
    """
    # this function smooth a profile by replacing each entry at position i with the average of the entries 
    #  in i-smooth_per_side, i+smooth_per_side.
    if use_savgol_filter is True then
        window_length = (2*smooth_per_side)+(1-(2*smooth_per_side)%2) (makes odd)
        and polyorder=min([5,win-1]) mode='interp'
    otherwise use_savgol_filter can be a tuple and then
        window_length,polyorder=use_savgol_filter
    """
    smoothed = []
    if weights is not None:
        if len(weights) != len(profile):
            sys.stderr.write(
                "**WARNING** in smooth_profile() len(weights)!=len(profile) %d %d\n"
                % (len(weights), len(profile))
            )
        profile = numpy.array(profile)
        weights = numpy.array(weights)
        for i in range(0, len(profile)):
            if i - smooth_per_side < 0:
                s = 0
            else:
                s = i - smooth_per_side
            den = sum(weights[s : i + smooth_per_side + 1])
            if den == 0:
                smoothed += [
                    0
                ]  # SHOULD raise warning but within CamSol this is used for solvent-exposure so for non-exposed regions 0 is the expected results
            else:
                smoothed += [
                    float(
                        sum(
                            profile[s : i + smooth_per_side + 1]
                            * weights[s : i + smooth_per_side + 1]
                        )
                    )
                    / den
                ]
    else:
        for i in range(0, len(profile)):
            if i - smooth_per_side < 0:
                s = 0
            else:
                s = i - smooth_per_side
            # print i, s, profile[s:i+smooth_per_side+1],len(profile[s:i+smooth_per_side+1]),float(sum(profile[s:i+smooth_per_side+1]))/len(profile[s:i+smooth_per_side+1])
            smoothed += [
                float(sum(profile[s : i + smooth_per_side + 1]))
                / len(profile[s : i + smooth_per_side + 1])
            ]
    return numpy.array(smoothed)


def nan_smooth_profile(profile, smooth_per_side=3, weights=None):
    """
DOES NOT CONSERVE WINDOW LENGTH AROUND nan (not used anyway)
    # this function smooth a profile by replacing each entry at position i with the average of the entries 
    #  in i-smooth_per_side, i+smooth_per_side.
    if use_savgol_filter is True then
        window_length = (2*smooth_per_side)+(1-(2*smooth_per_side)%2) (makes odd)
        and polyorder=min([5,win-1]) mode='interp'
    otherwise use_savgol_filter can be a tuple and then
        window_length,polyorder=use_savgol_filter
    """
    smoothed = []
    if weights is not None:
        if len(weights) != len(profile):
            sys.stderr.write(
                "**WARNING** in smooth_profile() len(weights)!=len(profile) %d %d\n"
                % (len(weights), len(profile))
            )
        profile = numpy.array(profile)
        weights = numpy.array(weights)
        for i in range(0, len(profile)):
            if numpy.isnan(profile[i]):
                smoothed += [numpy.nan]
                continue
            if i - smooth_per_side < 0:
                s = 0
            else:
                s = i - smooth_per_side
            den = numpy.nansum(weights[s : i + smooth_per_side + 1])
            if den == 0:
                smoothed += [
                    0
                ]  # SHOULD raise warning but within CamSol this is used for solvent-exposure so for non-exposed regions 0 is the expected results
            else:
                smoothed += [
                    float(
                        numpy.nansum(
                            profile[s : i + smooth_per_side + 1]
                            * weights[s : i + smooth_per_side + 1]
                        )
                    )
                    / den
                ]
    else:
        for i in range(0, len(profile)):
            if i - smooth_per_side < 0:
                s = 0
            else:
                s = i - smooth_per_side
            # print i, s, profile[s:i+smooth_per_side+1],len(profile[s:i+smooth_per_side+1]),float(sum(profile[s:i+smooth_per_side+1]))/len(profile[s:i+smooth_per_side+1])
            if numpy.isnan(profile[i]):
                smoothed += [numpy.nan]
            else:
                smoothed += [
                    numpy.nansum(profile[s : i + smooth_per_side + 1])
                    / nanlen(profile[s : i + smooth_per_side + 1])
                ]
    return numpy.array(smoothed)


class Residue_properties:
    def __init__(
        self,
        pKa=0,
        chargeAbovepKa=0,
        chargeBelowpKa=0,
        hydrophobicAbovepKa=0,
        hydrophobicBelowpKa=0,
        betaSheet=0,
        alphaHelix=0,
    ):
        self.pKa = pKa
        self.chargeAbovepKa = chargeAbovepKa
        self.chargeBelowpKa = chargeBelowpKa
        self.hydrophobicAbovepKa = hydrophobicAbovepKa
        self.hydrophobicBelowpKa = hydrophobicBelowpKa
        self.betaSheet = betaSheet
        self.alphaHelix = alphaHelix

    def get_all_properties(self, pH):
        """
        returns all properties according to pH as a list  (4 properties)
        [self.hydrophobicBelowpKa, self.alphaHelix, self.betaSheet, self.chargeBelowpKa]
        """
        if pH <= self.pKa:
            return [
                self.hydrophobicBelowpKa,
                self.alphaHelix,
                self.betaSheet,
                self.chargeBelowpKa,
            ]
        else:
            return [
                self.hydrophobicAbovepKa,
                self.alphaHelix,
                self.betaSheet,
                self.chargeAbovepKa,
            ]

    def charge(self, pH):
        if pH <= self.pKa:
            return self.chargeBelowpKa
        else:
            return self.chargeAbovepKa

    def hydrophobicity(self, pH):
        if pH <= self.pKa:
            return self.hydrophobicBelowpKa
        else:
            return self.hydrophobicAbovepKa


s = {
    "A": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.70,
        alphaHelix=0.79,
        hydrophobicAbovepKa=-0.39,
        hydrophobicBelowpKa=-0.39,
    ),
    "C": Residue_properties(
        pKa=8.33,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=0.84,
        alphaHelix=0.49,
        hydrophobicAbovepKa=2.25,
        hydrophobicBelowpKa=-0.25,
    ),
    "D": Residue_properties(
        pKa=3.65,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=0.68,
        alphaHelix=0.65,
        hydrophobicAbovepKa=3.81,
        hydrophobicBelowpKa=1.31,
    ),
    "E": Residue_properties(
        pKa=4.25,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=0.73,
        alphaHelix=0.97,
        hydrophobicAbovepKa=2.91,
        hydrophobicBelowpKa=1.11,
    ),
    "F": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.98,
        alphaHelix=0.77,
        hydrophobicAbovepKa=-2.27,
        hydrophobicBelowpKa=-2.27,
    ),
    "G": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.228,
        alphaHelix=0.184,
        hydrophobicAbovepKa=0.0693178842655,
        hydrophobicBelowpKa=0.0693178842655,
    ),
    "H": Residue_properties(
        pKa=6.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=0.80,
        alphaHelix=0.80,
        hydrophobicAbovepKa=0.64,
        hydrophobicBelowpKa=1.0,
    ),
    "I": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=1.00,
        alphaHelix=0.93,
        hydrophobicAbovepKa=-1.82,
        hydrophobicBelowpKa=-1.82,
    ),
    "K": Residue_properties(
        pKa=10.52,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=0.73,
        alphaHelix=0.88,
        hydrophobicAbovepKa=2.77,
        hydrophobicBelowpKa=2.77,
    ),
    "L": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.838,
        alphaHelix=0.828,
        hydrophobicAbovepKa=-1.46760381926,
        hydrophobicBelowpKa=-1.46760381926,
    ),
    "M": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.82,
        alphaHelix=0.82,
        hydrophobicAbovepKa=-0.96,
        hydrophobicBelowpKa=-0.96,
    ),
    "N": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.60,
        alphaHelix=0.61,
        hydrophobicAbovepKa=1.91,
        hydrophobicBelowpKa=1.91,
    ),
    "P": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.238,
        alphaHelix=0.260,
        hydrophobicAbovepKa=-0.426364344181,
        hydrophobicBelowpKa=-0.426364344181,
    ),
    "Q": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.876,
        alphaHelix=0.834,
        hydrophobicAbovepKa=1.65394825958,
        hydrophobicBelowpKa=1.65394825958,
    ),
    "R": Residue_properties(
        pKa=12.48,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=0.79,
        alphaHelix=0.95,
        hydrophobicAbovepKa=3.50,
        hydrophobicBelowpKa=3.50,
    ),
    "S": Residue_properties(
        pKa=13.6,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.67,
        alphaHelix=0.67,
        hydrophobicAbovepKa=1.24,
        hydrophobicBelowpKa=1.24,
    ),
    "T": Residue_properties(
        pKa=13.6,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.745,
        alphaHelix=0.551,
        hydrophobicAbovepKa=0.958057550723,
        hydrophobicBelowpKa=0.958057550723,
    ),
    "V": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=1.05,
        alphaHelix=0.75,
        hydrophobicAbovepKa=-1.3,
        hydrophobicBelowpKa=-1.3,
    ),
    "W": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.86,
        alphaHelix=0.64,
        hydrophobicAbovepKa=-2.13,
        hydrophobicBelowpKa=-2.13,
    ),
    "Y": Residue_properties(
        pKa=10.06,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=0.89,
        alphaHelix=0.73,
        hydrophobicAbovepKa=-1.47,
        hydrophobicBelowpKa=-1.47,
    )
    #'-' :Residue_properties(pKa=20,chargeAbovepKa= numpy.nan,chargeBelowpKa=numpy.nan,betaSheet=numpy.nan ,alphaHelix=  numpy.nan ,hydrophobicAbovepKa= numpy.nan,hydrophobicBelowpKa= numpy.nan )
}


AminoAcids = {
    "A": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.24849524205,
        alphaHelix=4.36944785247,
        hydrophobicAbovepKa=-0.39,
        hydrophobicBelowpKa=-0.39,
    ),
    "C": Residue_properties(
        pKa=8.33,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=4.43081679884,
        alphaHelix=3.89182029811,
        hydrophobicAbovepKa=2.25,
        hydrophobicBelowpKa=-0.25,
    ),
    "D": Residue_properties(
        pKa=3.65,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=4.21950770518,
        alphaHelix=4.1743872699,
        hydrophobicAbovepKa=3.81,
        hydrophobicBelowpKa=1.31,
    ),
    "E": Residue_properties(
        pKa=4.25,
        chargeAbovepKa=-1.0,
        chargeBelowpKa=0.0,
        betaSheet=4.29045944115,
        alphaHelix=4.5747109785,
        hydrophobicAbovepKa=2.91,
        hydrophobicBelowpKa=1.11,
    ),
    "F": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.58496747867,
        alphaHelix=4.34380542185,
        hydrophobicAbovepKa=-2.27,
        hydrophobicBelowpKa=-2.27,
    ),
    "G": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=3.12676053596,
        alphaHelix=2.91235066461,
        hydrophobicAbovepKa=0.0693178842655,
        hydrophobicBelowpKa=0.0693178842655,
    ),
    "H": Residue_properties(
        pKa=6.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=4.38202663467,
        alphaHelix=4.38202663467,
        hydrophobicAbovepKa=0.64,
        hydrophobicBelowpKa=1.0,
    ),
    "I": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.60517018599,
        alphaHelix=4.53259949315,
        hydrophobicAbovepKa=-1.82,
        hydrophobicBelowpKa=-1.82,
    ),
    "K": Residue_properties(
        pKa=10.52,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=4.29045944115,
        alphaHelix=4.47733681448,
        hydrophobicAbovepKa=2.77,
        hydrophobicBelowpKa=2.77,
    ),
    "L": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.42843300749,
        alphaHelix=4.41642806139,
        hydrophobicAbovepKa=-1.46760381926,
        hydrophobicBelowpKa=-1.46760381926,
    ),
    "M": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.40671924726,
        alphaHelix=4.40671924726,
        hydrophobicAbovepKa=-0.96,
        hydrophobicBelowpKa=-0.96,
    ),
    "N": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.09434456222,
        alphaHelix=4.11087386417,
        hydrophobicAbovepKa=1.91,
        hydrophobicBelowpKa=1.91,
    ),
    "P": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=3.16968558068,
        alphaHelix=3.25809653802,
        hydrophobicAbovepKa=-0.426364344181,
        hydrophobicBelowpKa=-0.426364344181,
    ),
    "Q": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.47278099794,
        alphaHelix=4.42364830936,
        hydrophobicAbovepKa=1.65394825958,
        hydrophobicBelowpKa=1.65394825958,
    ),
    "R": Residue_properties(
        pKa=12.48,
        chargeAbovepKa=0.0,
        chargeBelowpKa=1.0,
        betaSheet=4.36944785247,
        alphaHelix=4.5538768916,
        hydrophobicAbovepKa=3.50,
        hydrophobicBelowpKa=3.50,
    ),
    "S": Residue_properties(
        pKa=13.6,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.20469261939,
        alphaHelix=4.20469261939,
        hydrophobicAbovepKa=1.24,
        hydrophobicBelowpKa=1.24,
    ),
    "T": Residue_properties(
        pKa=13.6,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.31079912539,
        alphaHelix=4.00914971616,
        hydrophobicAbovepKa=0.958057550723,
        hydrophobicBelowpKa=0.958057550723,
    ),
    "V": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.65396035016,
        alphaHelix=4.31748811354,
        hydrophobicAbovepKa=-1.3,
        hydrophobicBelowpKa=-1.3,
    ),
    "W": Residue_properties(
        pKa=7.0,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.45434729625,
        alphaHelix=4.15888308336,
        hydrophobicAbovepKa=-2.13,
        hydrophobicBelowpKa=-2.13,
    ),
    "Y": Residue_properties(
        pKa=10.06,
        chargeAbovepKa=0.0,
        chargeBelowpKa=0.0,
        betaSheet=4.48863636973,
        alphaHelix=4.29045944115,
        hydrophobicAbovepKa=-1.47,
        hydrophobicBelowpKa=-1.47,
    )
    #'-' :Residue_properties(pKa=20,chargeAbovepKa= numpy.nan,chargeBelowpKa=numpy.nan,betaSheet=numpy.nan ,alphaHelix=  numpy.nan ,hydrophobicAbovepKa= numpy.nan,hydrophobicBelowpKa= numpy.nan )
}


random_peptide = {
    "sequence": "MTGSSYDESQLIKGLVYKPNKAQNSNGNRVDLTKEKLIVESLSGWEGSYAVLRDMTGAAIACVYTQLGPLMLAIYYPEPKGFNITSECWGFVLMTIYGSFYIDVYNTKHKESEICDQLHNYRRIPVVFGPGECLQNPRSNPAVALDCSPMVKPRNLRADISLPANVDTRVATAAYKDNGYSGLLIFRLPEVLTPKAPSYTSVPYVMLNSTAKQKCFEECAKDKLVQLPRIKVIETHAEKSCFIIEFDKAGFMKGDKKPEPDYFMQLRQTFRRGLKMESEEDVASEDHVYSFRFGSVYTDSLEYYTLLNQTIITIEDFPNKITPDDLRNGEYYKTSIEGMKQNARKGEERHVLQMPINTKLDLHSGQEPKVEWNGASSERTPAIVFDITRRGEMSNAAISFIRLLLKGETINKLKKCSARYNICSEGTESMAQSKTEDDKRMSPDGEHCGLNLAGGEPQKLYFSFEKALTKKLGGPSFNAAFDKFLAAEQMTDEMTMETDNLSFYDGFTLPYCYLVILVSSKMVVETLKPSSQKINSTTMGKKKLNLAFLSLRILASGKINGQGAPLSLIIVAAKKDDKYPIITVVWSELDIVITVDISSYQPVEGLLKDIGMYDKTRGDKKFKALLQAEAYKSTYDIISKVNETGWPIQVRTKPYDENDKDRCVANRIKYRKAGLHEQYLDAHSPGVALDMFETKEASSFFYEPESDLARNQNLATVAYIQLNLISYLLPCQKLAYGDLGEFKIQFEDPQFATTKEMLTQQALTKLSVPEGAAAIEFSLIRDIYARQFKRVTNIPSALPPYRSVRQCSRPQDVPRVTEFIHMEKNRGNMRTLDMLEVSRALGDINEHAKRIRSSFYTAALDYRGKKFLGSLYDSGGRSGPVTSFKLKLLMAPKNNDSTPDMERTRIQDEMVDLSVNELKTCMKSSSALKPSEATDDGIGKRETKRSPGGGQKAWHTMIIEDQLERPIIDTNEDAWESKFDAVSRTDVLMISWFGDTNMEGLIKPILAQDFICMDGLHEARSTVLEDGITFSGNHRKVRTNSEGILPVVRAAWWAHDDSKAIGIKTLDSMAVIGGYYRWTDYYRSVLADEKQPAPGTELFVKKHAVGYRQSEDPGAFTATDKNYFTQLNSDRITLITYGVDDQGYKGKLGLFAKKFLWQFPISDEEKDQANHDGSRANPLLQDLWTCHPIRIVTKQDQIIELTLVMKVSGLEKPTLLRVELASLYHSNFCKASILSRNSYLEQDARPKDFVSPPEKNRSTEGTNTGYVLFYRANEAQAKIFLQHELRIMTEVQDRRLARNSVQSNLSGKIWEANPGPSSGFNGVELAYLRAERERQTYMIGLSGPLPIGLTAFDAAPAVENAKEGPLQKDPASLGIRNYIFSLGVLIPLGAPIKASGLAASSALEHAGLVTQPGKFYRHLLIQVREPITSLKKLDYVDELNDSTKILLLVAAHEVLKPVILLRRAEAAIGLKKGIDSAPEASGINEVPAVCLIEFLKPLRNRPPIDVQQPLELVSFQSIESFYQAKSKGGPLLIQVFDCVVSQISADATKGKEEYEWPVQRFPKNALTCSSNKAMGRTIVITLPRIRADAQEPAYWALPQELHVLRGFYTKGFDNSFCSWQRDARNYNMKHVDTVPVEITDAQTQHMVKISYLAYFDSTTIESKNSGRRTVEYQNNMSREEYCWVNKNVIAIPTRWKLVPGVAQCDHYALDGIIQHKTKLEKAAINYRKYTVYTQNRFARDQIFPFRELPYTEVTSRSTPPPPMKVSFFTSWRELTAFDLGRGPAELNRHYTNKIPKEIENEDQLSRLIENTGLKISRKYLLPFDQISLYRFDEVSRAGLEIHPAADTNTPFGLGLSSGEPLAVTVFVAAIELYALYQGRRNSLFNRIVLSPASMTPGSPIERVCNLNHVQLNAASDNQNFAWPTLKKVIGGSIGAKVDAVHLAIMTKSNSLDFCAQTMDEAKAIVARVQNWFPVRGTPDVLPDDVCISWKRYCDGCQLDLFSVGLARGPAVDLEAIPQDKLVYEQLDEEEFNMSKQVYAQTLRELQPDESTSGYASSLGRRARFSEQFNTHDADPQEVRENRVTRPLKLKIVTSLEIKSNEQGVSYKLQRSQIPFSSPDHRKYSQRPDNVSMPPNEHTVLVAMPVYRVITMGKKSLKFTTVSMLLGELAAVKNGKGSIVEWCDHVTQLEPGLQATQLKRAYTLLEREIDVRYSHLAGSVFIEDMTSGCSTLINLKTQLFITQRAHEADTSIASGAIYIEAHSCESGDDRLSHLYGGGLRWSLAPNNAKTAKLATAPFFFLKKTEYKFISRAQGKQEQRTMDVLGQNTDANSPNDERVTELSERSTESGDQMEGHPRLKFWEFLFTFEILCRSDRLQLLRIMSNTPGQGCTYLIVSFTAASTVVQAEYKSRCEERQVKPVYGLKMKTHGSNFSVELEQCSMPDTDKSLQISRQGHCLDIYNYEISNVPKLPLDGMLLILFPVSEPSLWQRHQGAVPIILAQNADIAGINYMGGVIDVIQLGTLDYARYIVEASTGLKSGIPGEKAIEPNLSTLAPGSLGFMMRRHFNNQRIVSLIQKMQVKVDLKYPCGKGGNLDDWYIFWASEGDYDSGGEKVVCKEPHGGVFQDLSGADESLFNKVCDVPAHNNLTSAVLSLGALIHSFNGSAKGETQNIKTAQRCPAKQEAFVTRVNKRMSPVVTHEKLRFVSHCSGPDASSKIAELFMSVKDIPKNLMEIATVMESLAVVQYVKTAPDFCKKYKEAHVAMQLKFIVVDHVVSQMFTKAGEYQQTREGGIAYQGVDPLAACEDELERLNRSGKAGAATRTDGADLTQQERLPSKDNSVKVRKNCRGGFPAIVYRGTDSRGNINKGAHFVKIIIQCVGPEGQLMWEAPILDNAKKATNESLEVEYKKESELQFMAPLGEDGKVDDSNLRSYWTRAKAADLERALNTVGHRLKKYLGDDSIIDLEQAYTKEAMRPNPISEILVYIKATARIPRYVIRLRVSQKAHLKTHEIQQIPSLWWMKELIRTTAMIAESYDMIWVREGGADTDFSEILAAGQDDGYYMFSFQFETMIISQNWLWVQAFITLAIHNTADAPVRYPRVILLEAFDLQYSMVFYLVIPVITQANSDRAKKWHAPARMNDGYAVRFGYSDKASYICYETSFGTSAIGKALQAGKDILLQLGTRSNDIRILLAVAESLDKINRNRKDENGGVTGRDRNFNEEKLQLHLLDDVRAIIWKVVAPNDTQINGSTDESDDYVIAFEPTSDLVTNVGQAFDVGELTGSIMSTLRGAPIRSVVLLAQREIVYAETKGPFLSQLYPEAALKLNYLPEEIEDEQTDPNGYLVDSSEGNAENAFYVLPNRTIKARNKQASLGRDATSRVKDSQDDMFAPMKVAWRKLGSLVQAYRGGSGVQVPSEQVAIAGVKGPAGDPIEDSYKLPVAELVYSSRPCTCERQTIFYFALVVTKKDDDMHRSDEIGDVEVSLTREEVRSMPKTCQILLHKIKFPPPAVGMSLLEQGILELLTNVGEKTKVQLFINPSEEFVGATNEWQALVYEKDYGHLQRGRDIGAYNALFGSHQEMRLDIVEWNKQPSRVHSRVTKGEVTVHDTKQTFDDPKKPSGGMKKVEIPTINNNPHEWSMTLEYEYYPNRPDEKMYITASLKVLTIIQEKNTIQCEQLRLIDVNHTQLEVHVESKSRLTLDRLLGNQCSVVKGFHEIMDILSNYMNRTTYELERPFGNVTIKAPARWGSRENKVENSINWASKMERKNFMEATNYGMFLGPKDPGYDLLLVFSQSDAAYVCLKRGSPFILSRTRFCNADELDARAAGKDGKIQVTGCVERGVEHTIRVKSTNVKVPTPIVALLCIRKYPPTGFAVNLHTDSPLRKSRPPELLEMVGLVPMKSCASGKSYRDNVKNTSFGLLGFKYSITLSPRIANTPYFFKNYLSVEVWALVTAQTVDNIQVRRPYQLAGVQKTEELGGDKKQAKSTMNLKPTCNKYETEEKTF",
    "mean": None,
    "stdev": None,
    "pH": None,
}


memory_last_pattern = {"last_found": 0}


def hydrophobicHydrophilicPatternOfResidue(
    sequence,
    aa_index,
    pH,
    use_memory_last_pattern_value=True,
    amino_acid_properties=AminoAcids,
    tolerence_for_hydrophobicity=-0.5,
    tolerence_for_hydrophilicity=0.5,
):
    """
    calculate hydrophobic hydrophilic patterns of length 5
    return 1 if a pattern is found
    """
    if aa_index < 2 or aa_index >= len(sequence) - 2:
        memory_last_pattern[
            "last_found"
        ] = 0  # in case it is a new sequence it resets the memory
        return 0
    if use_memory_last_pattern_value and memory_last_pattern["last_found"] > 0.2:
        memory_last_pattern["last_found"] = 0
        return 0
    if "P" in sequence[aa_index - 2 : aa_index + 3]:
        return 0  # special provision for prolines - never found in patterns
    if (
        amino_acid_properties[sequence[aa_index - 2]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[sequence[aa_index - 1]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[sequence[aa_index]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[sequence[aa_index + 1]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[sequence[aa_index + 2]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
    ):
        memory_last_pattern["last_found"] = 1
        return 1
    if (
        amino_acid_properties[sequence[aa_index - 2]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[sequence[aa_index - 1]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[sequence[aa_index]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[sequence[aa_index + 1]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[sequence[aa_index + 2]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
    ):
        memory_last_pattern["last_found"] = 1
        return 1
    memory_last_pattern["last_found"] = 0
    return 0


def gapless_stretch(sequence, aa_index, Naa_per_side=2):
    """
    return the part of the sequecne from aa_index-Naa_per_side to aa_index+Naa_per_side included but ingoring gaps
     also returns the index on the returned stretch that was aa_index in the input sequence
    """
    j, n = aa_index - 1, aa_index + 1
    before, after = "", ""
    while j >= 0:
        if sequence[j] != "-":
            before += sequence[j]
            if len(before) >= Naa_per_side:
                break
        j -= 1
    while n < len(sequence):
        if sequence[n] != "-":
            after += sequence[n]
            if len(after) >= Naa_per_side:
                break
        n += 1
    return before[::-1] + sequence[aa_index] + after, len(before)


def nan_hydrophobicHydrophilicPatternOfResidue(
    sequence,
    aa_index,
    pH,
    use_memory_last_pattern_value=True,
    amino_acid_properties=AminoAcids,
    tolerence_for_hydrophobicity=-0.5,
    tolerence_for_hydrophilicity=0.5,
):
    """
    calculate hydrophobic hydrophilic patterns of length 5
    return 1 if a pattern is found
    """
    if sequence[aa_index] == "-":
        return numpy.nan
    if aa_index < 2 or aa_index >= len(sequence) - 2:
        memory_last_pattern[
            "last_found"
        ] = 0  # in case it is a new sequence it resets the memory
        return 0
    if use_memory_last_pattern_value and memory_last_pattern["last_found"] > 0.2:
        memory_last_pattern["last_found"] = 0
        return 0
    stretch, mid_index = gapless_stretch(sequence, aa_index, Naa_per_side=2)
    if "P" in stretch:
        return 0  # special provision for prolines - never found in patterns
    if (
        amino_acid_properties[stretch[mid_index - 2]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[stretch[mid_index - 1]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[stretch[mid_index]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[stretch[mid_index + 1]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[stretch[mid_index + 2]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
    ):
        memory_last_pattern["last_found"] = 1
        return 1
    if (
        amino_acid_properties[stretch[mid_index - 2]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[stretch[mid_index - 1]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[stretch[mid_index]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
        and amino_acid_properties[stretch[mid_index + 1]].hydrophobicity(pH)
        > tolerence_for_hydrophilicity
        and amino_acid_properties[stretch[mid_index + 2]].hydrophobicity(pH)
        < tolerence_for_hydrophobicity
    ):
        memory_last_pattern["last_found"] = 1
        return 1
    memory_last_pattern["last_found"] = 0
    return 0


def gate_function(dist, aa, pH, amino_acid_properties=AminoAcids):
    if dist < 2:
        return amino_acid_properties[aa].charge(
            pH
        )  # tabulated form of exp( - x^4 /200 )
    elif dist == 2:
        return amino_acid_properties[aa].charge(pH) * 0.92
    elif dist == 3:
        return amino_acid_properties[aa].charge(pH) * 0.667
    elif dist == 4:
        return amino_acid_properties[aa].charge(pH) * 0.278
    elif dist == 5:
        return amino_acid_properties[aa].charge(pH) * 0.144
    elif dist == 6:
        return amino_acid_properties[aa].charge(pH) * 0.044
    return 0


def gatekeepersOfResidue(sequence, aa_index, pH, amino_acid_properties=AminoAcids):
    """
    calculate gatekeeper value for aa_index
     exploits a tabulated form of the function exp( - x^4 /200 ) to smooth with a distance threshold the contribution of gatekeepers
    """
    gate = 0
    for j in range(aa_index - 6, aa_index + 7):
        if j < 0:
            continue
        if j >= len(sequence):
            break
        d = abs(aa_index - j)
        gate += gate_function(
            d, sequence[j], pH, amino_acid_properties=amino_acid_properties
        )
    return abs(gate)


def nan_gatekeepersOfResidue(sequence, aa_index, pH, amino_acid_properties=AminoAcids):
    """
    calculate gatekeeper value for aa_index
     exploits a tabulated form of the function exp( - x^4 /200 ) to smooth with a distance threshold the contribution of gatekeepers
    """
    if sequence[aa_index] == "-":
        return numpy.nan
    gate = 0
    left, right, j, i = 1, 0, aa_index - 1, aa_index
    while j >= 0:
        if sequence[j] != "-":
            gate += gate_function(
                left, sequence[j], pH, amino_acid_properties=amino_acid_properties
            )
            left += 1
            if left > 6:
                break
        j -= 1
    while i < len(sequence):
        if sequence[i] != "-":
            gate += gate_function(
                right, sequence[i], pH, amino_acid_properties=amino_acid_properties
            )
            right += 1
            if right > 6:
                break
        i += 1
    return abs(gate)


def nanlen(profile):
    return numpy.count_nonzero(~numpy.isnan(profile))


def get_scoreOLD(profile, low_th=-0.7, up_th=0.7):
    """
    as published, this score is very good for a very small numbers of mutations (1-2) localized in (or in the vicinity of) aggregation promoting regions.
    For a more general solubility ranking use get_scoreNEW()
    """
    return (profile[(profile < low_th) | (profile > up_th)]).sum() / len(profile)


#
def nan_get_scoreOLD(profile, low_th=-0.7, up_th=0.7):
    """
    as published, this score is very good for a very small numbers of mutations (1-2) localized in (or in the vicinity of) aggregation promoting regions.
    For a more general solubility ranking use get_scoreNEW()
    """
    return (profile[(profile < low_th) | (profile > up_th)]).sum() / nanlen(profile)


# the following where determined through an independent optimization of the neural-network based approach on historic measurements of aggregation rates.
new_new = [
    0.604287194505,
    -0.701442827953,
    0.963037609793,
    0.409513991578,
    0.686110422566,
    0.878544278758,
    0.0,
]
new_C = [
    0.602768701070882,
    -0.7054521794065463,
    0.9513753487882796,
    0.4004571570064588,
    0.32485962536814583,
    0.5181584799174997,
    0.0,
]
#         0                 1                     2                   3                    4                   5 (4 and 5 are just length normalizations so they don't matter when comparing protein variants.)
def get_scoreNEW(profile, coeffs_zscore=new_C):
    """
    new way of calculating the solubility score. This score is more precise in ranking according to their solubility 
     protein variants with a large number of mutations, wherever these are localized. 
     It performs relatively well also in ranking completely unrelated proteins.
    Note well - the profile of input should be -1* the CamSol profile (because of the way coefficients had been originally optimised)
    """
    sc = (
        (
            (profile[profile > coeffs_zscore[0]] - coeffs_zscore[0]) * coeffs_zscore[2]
        ).sum()
        + (
            (profile[profile < coeffs_zscore[1]] - coeffs_zscore[1]) * coeffs_zscore[3]
        ).sum()
    ) / (coeffs_zscore[4] * len(profile) ** coeffs_zscore[5])
    return (
        -1.0 * (sc - 8.0) / 4.4
    )  # prev(sc-0.34521)/0.2193 # mean and stdev of random peptides between 10 and 500 aa long


def nan_get_scoreNEW(profile, coeffs_zscore=new_C):
    """
    new way of calculating the solubility score. This score is more precise in ranking according to their solubility 
     protein variants with a large number of mutations, wherever these are localized. 
     It performs relatively well also in ranking completely unrelated proteins.
    """
    sc = (
        numpy.nansum(
            (profile[profile > coeffs_zscore[0]] - coeffs_zscore[0]) * coeffs_zscore[2]
        )
        + numpy.nansum(
            (profile[profile < coeffs_zscore[1]] - coeffs_zscore[1]) * coeffs_zscore[3]
        )
    ) / (coeffs_zscore[4] * nanlen(profile) ** coeffs_zscore[5])
    return (
        -1.0 * (sc - 8.0) / 4.4
    )  # prev(sc-0.34521)/0.2193 # mean and stdev of random peptides between 10 and 500 aa long


#                   Hydrophob      , alpha-hel         , beta-str          ,charges             ,pattern(0.7 is patfactor),  gatekeepers
profile_coeffs = [
    -0.7182318543722275,
    -3.122306452144164,
    4.7438896491092715,
    -0.4555325379008295,
    0.7 * 2.81642,
    -0.10063112005578956,
]


def camsol_intrinsic(
    sequence,
    weights=None,
    pH=7.0,
    full_return=False,
    profile_coeffs=profile_coeffs,
    amino_acid_properties=AminoAcids,
    normalise_for_peptide=True,
    use_gatekeepers=True,
    use_old_score=False,
    use_old_table=False,
):
    """
    this function runs camsol_intrinsic on a given amino acid sequence (only 20 standard amino acids capital letters)
    return intrinsic_profile, intrinsic_solublity_score
    if full_return:
        return intrinsic_profile, intrinsic_solublity_score, old_score_pm07, all_properties,residue_only_profile,gatekeepers,old_intrinsic_score_from_global_properties,global_hydrophobicity, global_alpha_propensity, global_beta_propensity, global_net_charge,N_patterns
    use_old_table is equivalent to use_old_score and it is there for back-compatibility
    use_old_score computes the intrinsic solubility score by using regions above 0.7 and below -0.7 on the intrinsic solubility profile
        and it may be more accurate in classifying single point mutations.  
    """
    if use_old_table:
        use_old_score = True  # back-comptatibility of keyword - don't use use_old_table for new applications.
    if len(sequence) < 7:
        raise Exception(
            "Sequence %s too short for camsol_intrinsic (len=%d minimum=7)\n"
            % (sequence, len(sequence))
        )
    all_properties = numpy.array(
        [
            amino_acid_properties[aa].get_all_properties(pH)
            + [
                hydrophobicHydrophilicPatternOfResidue(
                    sequence, j, pH, amino_acid_properties=amino_acid_properties
                )
            ]
            for j, aa in enumerate(sequence)
        ]
    )
    all_properties[..., 3] = numpy.abs(
        all_properties[..., 3]
    )  # absolute value of net charge
    residue_only_profile = (all_properties * profile_coeffs[:5]).sum(axis=1)
    smoothed_profile = smooth_profile(
        residue_only_profile, smooth_per_side=3, weights=weights
    )
    # calculate gatekeepers and add to profile
    if use_gatekeepers:
        gatekeepers = numpy.array(
            [
                gatekeepersOfResidue(
                    sequence, j, pH, amino_acid_properties=amino_acid_properties
                )
                for j in range(len(sequence))
            ]
        )
        smoothed_profile += profile_coeffs[5] * gatekeepers
    # do the random peptide normalisation of the profile
    if normalise_for_peptide:
        if (
            pH != random_peptide["pH"]
        ):  # normalise only once per import or again if the pH is changed from previous run
            peptide_profile, _ = camsol_intrinsic(
                random_peptide["sequence"],
                pH=pH,
                profile_coeffs=profile_coeffs,
                amino_acid_properties=amino_acid_properties,
                normalise_for_peptide=False,
                full_return=False,
            )
            random_peptide["mean"] = numpy.mean(-1 * peptide_profile)
            random_peptide["stdev"] = numpy.std(-1 * peptide_profile)
            random_peptide["pH"] = pH
        smoothed_profile = (smoothed_profile - random_peptide["mean"]) / random_peptide[
            "stdev"
        ]
    if full_return:
        # can calculate the intrinsic aggregation propensity from global properties (like in the time of old Zyggregator - but usually does not correlate well with experiments)
        (
            global_hydrophobicity,
            global_alpha_propensity,
            global_beta_propensity,
            global_net_charge,
            N_patterns,
        ) = all_properties.sum(axis=0)
        old_intrinsic_score_from_global_properties = (
            (
                profile_coeffs[0] * global_hydrophobicity
                + profile_coeffs[1] * global_alpha_propensity
                + profile_coeffs[2] * global_beta_propensity
            )
            / (10000.0 * numpy.tanh(len(sequence) / 10000.0))
            + profile_coeffs[3] * abs(global_net_charge)
            + profile_coeffs[4] * N_patterns
        )
        return (
            -1 * smoothed_profile,
            get_scoreNEW(smoothed_profile),
            -1.0 * get_scoreOLD(smoothed_profile),
            all_properties,
            residue_only_profile,
            gatekeepers,
            old_intrinsic_score_from_global_properties,
            global_hydrophobicity,
            global_alpha_propensity,
            global_beta_propensity,
            global_net_charge,
            N_patterns,
        )
    if use_old_score:
        score = -1.0 * get_scoreOLD(smoothed_profile)
    else:
        score = get_scoreNEW(smoothed_profile)
    return -1.0 * smoothed_profile, score


def contiguous_split(list_of_int, reverse=False):
    """
    split a list of int into contiguous bits
    [1,2,3,4,6,7,8,11,13,15,16] -> [[1, 2, 3, 4], [6, 7, 8], [11], [13], [15, 16]]
    """
    if reverse:
        a = -1
    else:
        a = 1
    conts = [[]]
    for j, i in enumerate(list_of_int[:-1]):
        if conts[-1] == []:
            conts[-1] += [i]
        if a * (list_of_int[j + 1] - i) == 1:
            conts[-1] += [list_of_int[j + 1]]
        else:
            conts += [[]]
    if conts[-1] == []:
        conts[-1] += [list_of_int[-1]]
    return conts


def get_APR(camsol_profile, low_threshold=-0.9, min_length=3):
    """
    get the aggregation promoting regions from a CamSol profiles
    :param camsol_profile:
    :param low_threshold:
    :param min_length:
    :return: APR_indices,APR_scores, sum(APR_scores)
    """
    if len(camsol_profile[camsol_profile <= -0.9]) == 0:
        return [], [], 0
    # gets the aggregation prone regions APRs from a profile
    APR_indices = [
        numpy.array(reg)
        for reg in contiguous_split(numpy.where(camsol_profile <= low_threshold)[0])
        if len(reg) >= min_length
    ]
    APR_scores = [
        camsol_profile[reg].sum() for reg in APR_indices
    ]  # all negative, the most negative the most aggregation prone
    return APR_indices, APR_scores, sum(APR_scores)


def _get_s2D_disorder_profile(output_matrix, population_low_threshold=0.45):
    """
    auxiliary function of s2D_corrected_CamSol
    """
    if output_matrix.shape[1] == 4:
        inds = numpy.array([0, 1, 3])  # helix, beta and PPII
    else:
        inds = inds = numpy.array([0, 1])
    ss_no_coil = output_matrix[:, inds]
    disorder = numpy.ones(len(ss_no_coil))
    ii = numpy.where(numpy.max(ss_no_coil, axis=1) >= population_low_threshold)
    disorder[ii] -= numpy.max(ss_no_coil, axis=1)[ii]
    return disorder


def _get_s2D_disorder_Helix_profile(output_matrix, population_low_threshold=0.45):
    """
    auxiliary function of s2D_corrected_CamSol
    """
    if output_matrix.shape[1] == 4:
        inds = numpy.array([1, 3])  # beta and PPII
    else:
        inds = inds = numpy.array([1])
    ss_no_coil = output_matrix[:, numpy.array(inds)]
    disorder = (
        numpy.ones(len(ss_no_coil)) - output_matrix[:, 0]
    )  # subtract helix as this is (almost) always aggregation-protecting
    ii = numpy.where(
        numpy.max(ss_no_coil, axis=1) >= population_low_threshold
    )  # further subtract most populated structure if not helix and if greater than population_low_threshold
    disorder[ii] -= numpy.max(ss_no_coil, axis=1)[ii]
    return disorder


def s2D_corrected_CamSol(
    s2D_output,
    sequences=None,
    seq_names=None,
    return_sequences=False,
    **kwargs_camsol_intrinsic
):
    """
    built to get one or more input sequences,
    can read a s2D_output file or instead one can give the s2D_output_matrix in its palce and corresponding sequences
    returns a list of lists in the form:
    [ [intrinsic_profile,intrinsic_score,s2D_corrected_profile,s2D_corrected_score, seq_name ]
       where seq_name may be None if not given as input (or not read from file).
    if return_sequences it also returns the sequence as the last element (useful if run on a s2D_output file)
    it uses _get_s2D_disorder_Helix_profile to get exposure weights
     the s2D methos is from Sormanni et al. J Mol Biol 2015 and it is 
     freely available for download at www-mvsoftware.ch.cam.ac.uk
    """
    results = []
    if type(s2D_output) is str:
        try:
            from . import s2D_class
        except Exception:
            import s2D_class  # raise error if not found in either way
        seq_names, sequences, s2D_output, ss_kind_str = s2D_class.read_output(
            s2D_output
        )
    elif (
        type(s2D_output) is not list
    ):  # One list per sequence # and not ( hasattr(s2D_output,'shape') and len(s2D_output.shape)>1 ): # if it is not a list and not an (at least) bidimensional numpy array
        s2D_output = [s2D_output]
    if type(sequences) is not list and type(sequences) is not tuple:
        sequences = [sequences]
    if type(seq_names) is not list and type(seq_names) is not tuple:
        if seq_names is None:
            seq_names = [""] * len(sequences)
        else:
            seq_names = [seq_names]
    for j, seq in enumerate(sequences):
        intrinsic_profile, intrinsic_score = camsol_intrinsic(
            seq, **kwargs_camsol_intrinsic
        )
        disorder = _get_s2D_disorder_Helix_profile(s2D_output[j])
        s2D_corrected_profile, _ = camsol_intrinsic(seq, weights=disorder)
        s2D_corrected_profile *= disorder  # DOUBLE CORRECT like in zsurf!!
        s2D_corrected_score = get_scoreNEW(-1 * s2D_corrected_profile)
        if return_sequences:
            results += [
                [
                    intrinsic_profile,
                    intrinsic_score,
                    s2D_corrected_profile,
                    s2D_corrected_score,
                    seq_names[j],
                    seq,
                ]
            ]
        else:
            results += [
                [
                    intrinsic_profile,
                    intrinsic_score,
                    s2D_corrected_profile,
                    s2D_corrected_score,
                    seq_names[j],
                ]
            ]
    return results


def seq_profile_to_alignment(
    profile, seq_from_alignment, gap_profile_value=numpy.nan, gap_symbol="-"
):
    """
    given a sequence profile (list of floats) and a sequence as read from an alignemt file
     it adds in the profile list entries equal to gap_profile_value
     at positions corresponding to the gaps in the sequences '-' 
    """
    if len(profile) == 2 and hasattr(
        profile[0], "__len__"
    ):  # most likely a CI from errors
        return [
            seq_profile_to_alignment(
                pr,
                seq_from_alignment,
                gap_profile_value=gap_profile_value,
                gap_symbol=gap_symbol,
            )
            for pr in profile
        ]
    al_profile = []
    i = 0
    for res in seq_from_alignment:
        if res != gap_symbol:
            al_profile += [float(profile[i])]
            i += 1
        else:
            al_profile += [gap_profile_value]
    if len(profile) != i:
        sys.stderr.write(
            "WARNING in seq_profile_to_alignment() len(profile)!= number of non-gap residues in seq_from_alignment (%d!=%d)\n"
            % (len(profile), i)
        )
    return numpy.array(al_profile)


def camsol_intrinsic_gaps(
    sequence_with_gaps,
    weights=None,
    pH=7.0,
    gap_symbol="-",
    return_also_gapless_profile=False,
    full_return=False,
    profile_coeffs=profile_coeffs,
    amino_acid_properties=AminoAcids,
    normalise_for_peptide=True,
    use_gatekeepers=True,
    use_old_score=False,
    use_old_table=False,
):
    """
    this function runs camsol_intrinsic on a given amino acid sequence (only 20 standard amino acids capital letters)
    return intrinsic_profile, intrinsic_solublity_score
    if full_return:
        return intrinsic_profile, intrinsic_solublity_score, old_score_pm07, all_properties,residue_only_profile,gatekeepers,old_intrinsic_score_from_global_properties,global_hydrophobicity, global_alpha_propensity, global_beta_propensity, global_net_charge,N_patterns
    use_old_table is equivalent to use_old_score and it is there for back-compatibility
    use_old_score computes the intrinsic solubility score by using regions above 0.7 and below -0.7 on the intrinsic solubility profile
        and it may be more accurate in classifying single point mutations.  
    """
    if use_old_table:
        use_old_score = True  # back-comptatibility of keyword

    return_aligned_profiles = False
    if gap_symbol in sequence_with_gaps:
        return_aligned_profiles = True
        sequence = numpy.array(list(sequence_with_gaps))
        sequence = sequence[sequence != gap_symbol]
    else:
        sequence = sequence_with_gaps

    if len(sequence) < 7:
        raise Exception(
            "Sequence (gapless) too short for camsol_intrinsic (len=%d minimum=7)\n"
            % (len(sequence))
        )

    all_properties = numpy.array(
        [
            amino_acid_properties[aa].get_all_properties(pH)
            + [
                hydrophobicHydrophilicPatternOfResidue(
                    sequence, j, pH, amino_acid_properties=amino_acid_properties
                )
            ]
            for j, aa in enumerate(sequence)
        ]
    )
    if full_return:
        # can calculate the intrinsic aggregation propensity from global properties (in a manner similar to that of old Zyggregator - but usually does not correlate well with experimental measurements)
        (
            global_hydrophobicity,
            global_alpha_propensity,
            global_beta_propensity,
            global_net_charge,
            N_patterns,
        ) = numpy.sum(all_properties, axis=0)
        old_intrinsic_score_from_global_properties = (
            (
                profile_coeffs[0] * global_hydrophobicity
                + profile_coeffs[1] * global_alpha_propensity
                + profile_coeffs[2] * global_beta_propensity
            )
            / (10000.0 * numpy.tanh(len(sequence) / 10000.0))
            + profile_coeffs[3] * abs(global_net_charge)
            + profile_coeffs[4] * N_patterns
        )

    # residue_only_profile = (all_properties*profile_coeffs[:5]).sum(axis=1)
    all_properties[..., 3] = numpy.abs(
        all_properties[..., 3]
    )  # absolute value of net charge
    residue_only_profile = numpy.sum(all_properties * profile_coeffs[:5], axis=1)
    smoothed_profile = smooth_profile(
        residue_only_profile, smooth_per_side=3, weights=weights
    )
    print(
        (
            "len(residue_only_profile)",
            len(residue_only_profile),
            "len(smoothed_profile)",
            len(smoothed_profile),
            "all_properties.shape",
            all_properties.shape,
        )
    )
    # calculate gatekeepers and add to profile
    if use_gatekeepers:
        gatekeepers = numpy.array(
            [
                gatekeepersOfResidue(
                    sequence, j, pH, amino_acid_properties=amino_acid_properties
                )
                for j in range(len(sequence))
            ]
        )
        smoothed_profile += profile_coeffs[5] * gatekeepers
    # do the random peptide normalisation of the profile
    if normalise_for_peptide:
        if pH != random_peptide["pH"]:  # changed from last time saved
            peptide_profile, _ = camsol_intrinsic(
                random_peptide["sequence"],
                pH=pH,
                profile_coeffs=profile_coeffs,
                amino_acid_properties=amino_acid_properties,
                normalise_for_peptide=False,
                full_return=False,
            )
            random_peptide["mean"] = numpy.mean(-1 * peptide_profile)
            random_peptide["stdev"] = numpy.std(-1 * peptide_profile)
            random_peptide["pH"] = pH
        smoothed_profile = (smoothed_profile - random_peptide["mean"]) / random_peptide[
            "stdev"
        ]
    if full_return or use_old_score:
        oldscore = -1.0 * get_scoreOLD(smoothed_profile)
    score = get_scoreNEW(smoothed_profile)
    if return_aligned_profiles:
        smoothed_profile_al = seq_profile_to_alignment(
            smoothed_profile,
            sequence_with_gaps,
            gap_profile_value=numpy.nan,
            gap_symbol=gap_symbol,
        )
    else:
        smoothed_profile_al = smoothed_profile
    if full_return:
        return (
            -1 * smoothed_profile_al,
            score,
            oldscore,
            all_properties,
            residue_only_profile,
            gatekeepers,
            old_intrinsic_score_from_global_properties,
            global_hydrophobicity,
            global_alpha_propensity,
            global_beta_propensity,
            global_net_charge,
            N_patterns,
        )
    if use_old_score:
        score = oldscore
    if return_also_gapless_profile:
        return -1.0 * smoothed_profile_al, score, smoothed_profile
    return -1.0 * smoothed_profile_al, score


def camsol_intrinsic_multichain(sequences, score_function=get_scoreNEW, **kwargs):
    """
    runs CamSol intrinsic on a multi-chain complex (e.g. a full length antibody) 
    the profiles are separate and identical to those that would be obtained by running camsol_intrinsic
    on each chain individually. The intrinsic solubility score however is calculated by merging all profiles of the various
    chains first.
    return score, profiles,scores
    """
    if hasattr(sequences, "keys"):
        seqs = [sequences[k] for k in sequences]
    else:
        seqs = sequences
    profiles, scores = [], []
    merg = numpy.empty(0)
    for seq in seqs:
        p, s = camsol_intrinsic(seq, **kwargs)
        profiles += [p]
        scores += [s]
        merg = numpy.hstack((merg, profiles[-1]))
        # if larger_values_more_soluble :profiles[-1]*=-1.
    score = score_function(-1.0 * merg)
    # if larger_values_more_soluble : score*= -1.
    return score, profiles, scores


def camsol_intrinsic_multichain_gaps(sequences, score_function=get_scoreNEW, **kwargs):
    """
    runs CamSol intrinsic on a multi-chain complex (e.g. a full length antibody) 
    the profiles are separate and identical to those that would be obtained by running camsol_intrinsic
    on each chain individually. The intrinsic solubility score however is calculated by merging all profiles of the various
    chains first.
    return score, profiles,scores
    """
    if hasattr(sequences, "__keys__"):
        seqs = [sequences[k] for k in sequences]
    else:
        seqs = sequences
    profiles, scores = [], []
    merg = numpy.empty(0)
    for seq in seqs:
        p, s, gp = camsol_intrinsic_gaps(
            seq, return_also_gapless_profile=True, **kwargs
        )
        profiles += [p]
        scores += [s]
        merg = numpy.hstack(
            (merg, gp)
        )  # merge profiles without nan (the gapless) but return the aligned
        # if larger_values_more_soluble :profiles[-1]*=-1.
    score = score_function(-1.0 * merg)
    # if larger_values_more_soluble : score*= -1.
    return score, profiles, scores


def is_protein(seq, aa_list=AminoAcids):
    """
    check if a sequence is an amino acid sequence
    return True,None
    or return False, non_prot_res
    """
    for aa in seq:
        if aa not in aa_list:
            return False, aa
    return True, None


class Seq:
    # mimic of the Biopython Seq class within SeqRecord object, but in this way there is no need to have biopython installed
    def __init__(self, seq="", seq_id="", seq_name="", description=""):
        self.seq = seq
        self.id = seq_id
        self.name = seq_name
        self.description = description

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __repr__(self):
        restr = "Seq:%s:" % (self.name)
        if len(self.seq) > 20:
            restr += self.seq[:15] + "..." + self.seq[-2:]
        else:
            restr += self.seq
        return restr

    def __getslice__(self, i, j):
        return self.seq[i:j]

    def __getitem__(self, y):
        return self.seq[y]

    def __add__(self, y):
        return Seq(
            seq=self.seq + str(y),
            seq_id=self.id,
            seq_name=self.name,
            description=self.description,
        )


def parse_input(f):
    records = []
    i = 1
    fid = None
    if os.path.isfile(f):  # threat it as a fasta formatted file
        fid = f.split("/")[-1].split(".")[0]
        with open(f) as fff:
            content = fff.read()
        if ">" not in content:
            sys.stderr.write(
                "\n**WARNING** input file %s does not appear to be in fasta format - will very likely generate errors\n\n"
                % (f)
            )
        content = content.split(">")  # assume fasta file and quickly parse it
        for seq_bit in content:
            s = seq_bit.split("\n")
            if len(s) <= 1:
                continue
            records += [
                Seq(
                    seq_name=s[0].strip(),
                    seq_id=s[0].strip(),
                    seq="".join(s[1:]).replace(" ", ""),
                )
            ]
            # print 'r:',records[-1].id,records[-1]
        # print '%d sequences in file %s' % (len(records),f)
    else:
        s = ""
        if len(f) < 7:
            sys.stderr.write(
                "**ERROR** input sequence too short (len=%d) Skipping!! at least 7 amino acids required\n"
                % (len(f))
            )
            sys.exit(1)
        for letter in f:  # threat it as a sequence
            if letter.upper() in [
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
            ]:
                s += letter.upper()
            else:
                sys.stderr.write(
                    "NOT STANDARD RESIDUE, deleting it: letter '%s' in sequence %s\n"
                    % (letter, f)
                )
        records.append(Seq(s, "seq_" + str(i), "seq_" + str(i)))
        i += 1
    return records, fid


def main(
    form_input,
    outfilename=None,
    pH=7.0,
    max_sequences=None,
    plot_single_sequence=True,
    set_publish=False,
    glob_path="",
    run_id="",
    write_tmp_file=True,
):
    if "module" not in str(type(plotter)).lower():
        sys.stderr.write("Cannot plot! [%s]\n" % (str(type(plotter))))
        plot_single_sequence = False
    records, fid = parse_input(form_input)
    if fid == None:
        one_seq = True
    else:
        one_seq = False
    if len(glob_path) > 0 and glob_path[-1] != "/":
        glob_path += "/"
    if outfilename is None:
        outfilename = glob_path + "CamSol_intrinsic" + run_id + ".txt"
    if one_seq:
        profile, score = camsol_intrinsic(
            str(records[0].seq), pH=pH
        )  # use linear module
        outstr = "The protein variant intrinsic solubility score is %lf\n" % (score)

        out = open(outfilename, "w")
        out.write("Pos\taa\taa intrinsic solubility\n")
        for j, sc in enumerate(profile):
            out.write("%-2d\t%s\t%lf\n" % (j + 1, records[0].seq[j], sc))
        out.close()
        if plot_single_sequence:
            plotter.default_figure_sizes["dpi"] = 100
            if set_publish:
                plotter.set_publish(small_figure=False)
            plotter.plot_seq_profile(
                str(records[0].seq),
                profile,
                ylabel="Intrinsic residue solubility",
                zygg_like_lines=True,
                save=outfilename[:-4] + ".png",
                linewidth=2.0,
                add_colorline=(plotter.camsol_to_color_val, plotter.camsol_colormap),
            )
    else:
        if max_sequences is not None and max_sequences < len(records):
            sys.stderr.write(
                "\n**ERROR** because of limited computational resources the maximum number of input sequences is %d. Found %d\n\n"
                % (max_sequences, len(records))
            )
            sys.exit(1)
        out = open(outfilename, "w")
        out.write("Name\tprotein variant score\tintrinsic solubility profile\n")
        for rec in records:
            if len(rec.seq) < 7:
                sys.stderr.write(
                    "Warning %s too short (len=%d) Skipping!!\n"
                    % (str(rec.id), len(rec.seq))
                )
                continue
            ok, probl = is_protein(rec.seq)
            if not ok:
                sys.stderr.write(
                    "Warning %s contains non-standard AA '%s'.. Skipping!!\n"
                    % (str(rec.id), probl)
                )
                continue

            profile, score = camsol_intrinsic(str(rec.seq), pH=pH)  # use linear module
            profile = ";".join(map(str, profile))
            out.write("%s\t%lf\t%s\n" % (str(rec.id), score, profile))
            if (
                plot_single_sequence and len(records) < 5
            ):  # plot up to 4 sequences (e.g. light and heavy chain of an antibody, or chains of a complex), if more don't plot anything
                plotter.default_figure_sizes["dpi"] = 100
                if set_publish:
                    plotter.set_publish(small_figure=False)
                nam = (
                    outfilename[:-4]
                    + (
                        r.id.replace(" ", "_")
                        .replace("/", "_")
                        .replace("\\", "")
                        .replace("|", "_")
                        .replace(":", "-")
                        .replace("<", "-"),
                        replace(">", "-")
                        .replace("*", "-")
                        .replace("?", "-")
                        .replace('"', "")
                        .replace("'", "")
                        .replace("`", "")
                        .strip(),
                    )
                    + ".png"
                )
                plotter.plot_seq_profile(
                    str(records[0].seq),
                    profile,
                    ylabel="Intrinsic residue solubility",
                    zygg_like_lines=True,
                    save=nam,
                    title=r.id,
                    linewidth=2.0,
                    add_colorline=(
                        plotter.camsol_to_color_val,
                        plotter.camsol_colormap,
                    ),
                )
        out.close()
        df = pd.read_csv(outfilename, sep="\t")
        df = df[["Name", "protein variant score"]].rename(
            columns={"Name": "binder_id", "protein variant score": "camsol_score"}
            )
        df.to_csv(outfilename.replace(".txt", ".csv"), index=False)
        outstr = "CamSol intrinsic run on %d sequences completed!\n" % (len(records))
        
    if write_tmp_file and (glob_path != "" or run_id != ""):
        out = open(glob_path + "tmp" + run_id, "w")
        out.write(outstr)
        out.close()
    print(outstr)
    return outstr, outfilename


PAPER = None
if __name__ == "__main__":
    import argparse

    if len(sys.argv) < 2:
        sys.stderr.write("\n**ERROR**\nGIVEN NO INPUT sequence or fasta_file\n\n")
    parser = argparse.ArgumentParser(
        description="Runs the CamSol intrinsic program, which carries out sequence-based prediction of solubility.",
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        metavar="input sequence or input fasta file containing one or multiple sequences",
        type=str,
        nargs="+",
        help="input sequence (1-letter capital, 20 standard amino acids only) or input fasta file containing one or multiple sequences. A second positional arguments can be given and it will be processed as output file name (also with -out).",
    )
    parser.add_argument(
        "-out",
        "--outfilename",
        dest="outfilename",
        type=str,
        default=None,
        help=" Output file name that will store the intrinsic solubility profile",
    )
    parser.add_argument(
        "-pH",
        "--pH",
        dest="pH",
        type=float,
        default=7.0,
        help=" pH value to run the calculations (biophysical properties of amino acids are extracted from table and are different below and above their pKa)",
    )
    parser.add_argument(
        "-plot",
        "--plot",
        dest="plot_single_sequence",
        action="store_true",
        default=False,
        help=" whether to plot the intrinsic profile for the input sequence (not available for multiple sequences)",
    )
    parser.add_argument(
        "-ms",
        "--max_sequences",
        dest="max_sequences",
        type=int,
        default=1e8,
        help=" maximum number of sequences that can be processed in an individual fasta file (set to small number on web server)",
    )
    # webserver:
    parser.add_argument(
        "-path",
        "--path",
        dest="glob_path",
        type=str,
        default="",
        help=" used in web server mode to indicate a path with write/read permission where output should be stored",
    )
    parser.add_argument(
        "-run_id",
        "--run_id",
        dest="run_id",
        type=str,
        default="",
        help=" used in web server mode to indicate a uniqe id for this run, appended to ouptut file names to avoid overwritings",
    )

    arguments = parser.parse_args()
    if (
        len(arguments.input_file) > 1
        and arguments.outfilename == "CamSol_intrinsic_reults.txt"
    ):
        arguments.outfilename = arguments.input_file[
            1
        ]  # assumes second argument is output file
    arg_dict = dict(
        (k, v) for k, v in list(vars(arguments).items()) if k != "input_file"
    )
    outstr, outfilename = main(arguments.input_file[0], **arg_dict)

    if PAPER is not None:
        sys.stdout.write(
            "\nIf you find CamSol helpful for your research please cite:\n%s\n"
            % (PAPER)
        )
    sys.exit(0)
