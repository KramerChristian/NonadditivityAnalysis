# Nonadditivity Analysis
#
# Copyright (c) 2019, F. Hoffmann-La Roche Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

########################
# BEFORE RUNNING THIS CODE
#
# This code uses other external codes. To make it run properly, you have to adjust
# the paths below to the corresponding installations on your system.

import rdkit
import os.path
from rdkit import RDConfig

salt_defns = os.path.join(RDConfig.RDDataDir, "Salts.txt")  # replace if you have more specific definitions

font_path = "arial.ttf"   # Only used in draw_pics, currently not supported

########################

import subprocess
import math
import pickle
import argparse
import random

from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import SaltRemover
from rdkit.Chem import Descriptors

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

from scipy import stats
import numpy as np


def numpy_std(values):
    return np.std(values)


def mad_std(values):
    """
    Estimate the standard deviation based on the median absolute deviation (MAD)
    In Gaussian distributions, SD = 1.4826 * MAD

    This is a fast and simple robust SD estimator. However, it may have problems
    for small datasets where the average is not well estimated.
    """

    median = np.median(values)
    ad = np.abs(values-median)
    mad = np.median(ad)

    return 1.4826 * mad


def Sn_MedMed_std(values):
    """
    Estimate the standard deviation based on the median of median pairwise differences
    as defined by Rousseeuw and Croux [1]:

    Sn = 1.1926 * med_i( med_j(abs(xi - xj)))

    [1] Rousseeuw, Peter J.; Croux, Christophe (December 1993),
    "Alternatives to the Median Absolute Deviation",
    Journal of the American Statistical Association, American Statistical Association,
    88 (424): 1273–1283, doi:10.2307/2291267, JSTOR 2291267
    """

    pairwise_medians = np.empty(len(values))

    for idx, i_value in enumerate(values):
        # Do I need to remove the case where i==j ? Tests indicate yes
        j_values = np.delete(values, [idx])
        abs_pairwise_diffs = np.abs(i_value - j_values)
        pairwise_medians[idx] = np.median(abs_pairwise_diffs)

    Sn = 1.1926 * np.median(pairwise_medians)

    return Sn


def build_neighbor_dictionary(mmps, no_chiral=False):
    """
    Generate neighbor dictionary from calculated MMPs
    """
    print("Analyzing neighborhoods")

    neighs = {}
    for line in mmps:
        smiles_lhs, smiles_rhs, id_lhs, id_rhs, transf, const = line.split("\t")
        if no_chiral and "@" in transf:
            continue
        var_lhs, var_rhs = transf.split(">>")
        # Skip pair if the transformation has more than one anchoring point
        # and the topological distance changes between those two (no reason to assume additivity then)
        if "[*:2]" in var_lhs:
            a = Chem.MolFromSmarts(var_lhs)
            b = Chem.MolFromSmarts(var_rhs)
            a_idx1 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:1]")
            a_idx2 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:2]")
            b_idx1 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:1]")
            b_idx2 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:2]")
            if not Chem.GetDistanceMatrix(a)[a_idx1, a_idx2] == Chem.GetDistanceMatrix(b)[b_idx1, b_idx2]:
                continue
            if "[*:3]" in var_lhs:
                a_idx3 = [atom.GetSmarts() for atom in a.GetAtoms()].index("[*:3]")
                b_idx3 = [atom.GetSmarts() for atom in b.GetAtoms()].index("[*:3]")
                if not Chem.GetDistanceMatrix(a)[a_idx1, a_idx3] == Chem.GetDistanceMatrix(b)[b_idx1, b_idx3]:
                    continue
                if not Chem.GetDistanceMatrix(a)[a_idx2, a_idx3] == Chem.GetDistanceMatrix(b)[b_idx2, b_idx3]:
                    continue
        # Add to neighbor dictionary
        if id_lhs in neighs.keys():
            if id_rhs not in [i[0] for i in neighs[id_lhs]]:
                neighs[id_lhs].append((id_rhs, transf))
            else:
                id_rhs_idx = [i[0] for i in neighs[id_lhs]].index(id_rhs)
                old_transf_len = len(neighs[id_lhs][id_rhs_idx][1])
                if len(transf) < old_transf_len:
                    neighs[id_lhs][id_rhs_idx] = (id_rhs, transf)
        else:
            neighs[id_lhs] = [(id_rhs, transf)]

    return neighs


def get_circles(neighs):
    """
    Assemble circle information
    """
    print("Assembling circles")

    still_in = {key: True for key in neighs.keys()}
    neighbors = {key: [i[0] for i in dat] for key, dat in neighs.items()}

    trans = {}    
    for key, dat in neighs.items():
        for partner, T1 in dat:
            trans[(key, partner)] = T1

    trans_pairs = {}
    for pair, T1 in trans.items():
        trans_pairs[T1] = trans_pairs.get(T1, [])
        trans_pairs[T1].append(pair)

    circles = []
    for C1 in iter(neighs.keys()):
        still_in_2 = still_in.copy()
        for C2 in neighbors[C1]:
            if not still_in[C2]: continue
            for C3, C4 in trans_pairs[trans[(C1, C2)]]:
                if not still_in_2[C3]: continue
                if not still_in_2[C4]: continue
                if C3 not in neighbors[C1]: continue
                if C4 not in neighbors[C2]: continue
                if C2 == C3: continue
                if C1 == C4: continue
                if trans[(C1, C3)] == trans[(C2, C4)]: circles.append((C1, C2, C4, C3))
            still_in_2[C2] = False
        still_in[C1] = False

    return circles


def write_circles_to_output(circles, meas, neighs, outfile, props, units, images=False,
                            include_censored=False, update=False, series_column=None):
    """
    Get Diffs and write to output
    """
    print("Writing Output.")

    #########
    # Assemble header

    header = "Compound1\tCompound2\tCompound3\tCompound4\t" \
             + "SMILES1\tSMILES2\tSMILES3\tSMILES4\tSeries\t" \
             + "Transformation1\tTransformation2\tProperty\t"

    for i in range(4):
        header = header + "Prop_Cpd" + str(i + 1) + "\t"

    header = header + "Nonadditivity\t"
    if images:
        header = header + "Circle_image\t"
        if not os.path.isdir("images"):
            os.mkdir("images")
        cpds = [Chem.MolFromSmiles(i["smiles"]) for i in meas.values()]
        mcss_tot = rdFMCS.FindMCS(cpds, completeRingsOnly=True, timeout=60)
        if mcss_tot.numAtoms < 7:
            mcss_tot = None
        else:
            mcss_tot = Chem.MolFromSmarts(mcss_tot.smartsString)
            AllChem.Compute2DCoords(mcss_tot)

    header = header + "Circle_ID\tTheo_Quantile\n"

    nonadd_percompound_file = outfile[:outfile.index(".")] + "_perCompound.txt"
    napc_header = "Compound_ID\tSMILES\tSeries\tProperty\tOperator\tMeasured\t"
    if series_column:
        napc_header = napc_header + "Nonadd_pure\tnOccurence_pure\tNonadd_SD_pure\t"
        napc_header = napc_header + "Nonadd_mixed\tnOccurence_mixed\tNonadd_SD_mixed\n"
    else:
        napc_header = napc_header + "Nonadd_pC\tnOccurence\tNonadd_SD\n"

    circle_to_cpd_file = outfile[:outfile.index(".")] + "_c2c.txt"
    c2c_header = "Circle_ID\tCompound_ID\n"

    with open(outfile, "w") as f, open(nonadd_percompound_file, "w") as g, open(circle_to_cpd_file, "w") as h:
        f.write(header)
        g.write(napc_header)
        h.write(c2c_header)

        print("Estimated Experimental Uncertainty")

        #########
        # Assemble and fill output lines
        for target_idx, target in enumerate(props):
            outlines = []
            add_diffs = []
            nonadd_percompound = {}
            c2c_lines = []
            series_id = []

            print("for property: ", target)

            for circle_idx, cpds in enumerate(circles):
                ############
                # 1st: Randomly reorder circles
                min_idx = random.choice(range(4))     # This does a random selection of the starting compound

                circles[circle_idx] = (circles[circle_idx] + circles[circle_idx])[min_idx:min_idx + 4]  # Reorder Ids

                smis = [meas[cpd]["smiles"] for cpd in cpds]  # Reorder SMILES
                smis = (smis + smis)[min_idx:min_idx + 4]

                acts = [meas[cpd]["pAct"][target_idx] for cpd in cpds]  # Reorder Activities
                if "NA" in acts:
                    continue

                cens_signs = [meas[cpd]["qualifiers"][target_idx] for cpd in cpds]
                if include_censored:
                    if len([i for i in cens_signs if not i == ""]) > 1:
                        continue

                else:
                    if not set(cens_signs) == {""}:
                        continue

                acts = (acts + acts)[min_idx:min_idx + 4]

                ############
                # 2nd: Write Info to output lines
                line = "\t".join(circles[circle_idx]) + "\t"  # Cpd IDs
                line = line + "\t".join(smis) + "\t"  # SMILES

                if series_column:
                    if len(set([meas[cpds[0]]["series"],
                                meas[cpds[1]]["series"],
                                meas[cpds[2]]["series"],
                                meas[cpds[3]]["series"]])) == 1:
                        line = line + meas[cpds[0]]["series"] + "\t"
                        series_id.append(meas[cpds[0]]["series"])
                        pure_series = True
                    else:
                        line = line + " \t"
                        series_id.append("")
                        pure_series = False
                else:
                    line = line + " \t"

                transf_idx1 = [j[0] for j in neighs[circles[circle_idx][0]]].index(circles[circle_idx][1])
                line = line + neighs[circles[circle_idx][0]][transf_idx1][1] + "\t"
                transf_idx2 = [j[0] for j in neighs[circles[circle_idx][1]]].index(circles[circle_idx][2])
                line = line + neighs[circles[circle_idx][1]][transf_idx2][1] + "\t"

                line = line + target + "\t"

                line = line + "{:.3g}".format(acts[0]) + "\t"
                line = line + "{:.3g}".format(acts[1]) + "\t"
                line = line + "{:.3g}".format(acts[2]) + "\t"
                line = line + "{:.3g}".format(acts[3]) + "\t"

                nonadd = (acts[2] - acts[3]) - (acts[1] - acts[0])
                add_diffs.append(nonadd)
                line = line + "{:.3g}".format(nonadd) + "\t"

                for idx, cpd in enumerate(cpds):
                    sign_switch = 1
                    if min_idx in [1, 3]:
                        sign_switch = -1
                    if idx in [1, 3]:
                        sign_switch = -1 * sign_switch
                    if series_column:
                        if cpd in nonadd_percompound:
                            if pure_series:
                                nonadd_percompound[cpd].append((sign_switch * nonadd, "pure"))
                            else:
                                nonadd_percompound[cpd].append((sign_switch * nonadd, "mixed"))
                        else:
                            if pure_series:
                                nonadd_percompound[cpd] = [(sign_switch * nonadd, "pure")]
                            else:
                                nonadd_percompound[cpd] = [(sign_switch * nonadd, "mixed")]
                    else:
                        if cpd in nonadd_percompound:
                            nonadd_percompound[cpd].append(sign_switch * nonadd)
                        else:
                            nonadd_percompound[cpd] = [sign_switch * nonadd]

                if images:
                    image_file = "_".join(["images/AddCyc", target] + list(circles[circle_idx]) + [".png"])
                    line = line + image_file + "\t"
                    if update and os.path.exists(image_file):
                        continue
                    else:
                        draw_image(ids=circles[circle_idx],
                                   smiles=smis,
                                   tsmarts=[neighs[circles[circle_idx][0]][transf_idx1][1],
                                            neighs[circles[circle_idx][1]][transf_idx2][1]],
                                   pActs=["{:.3g}".format(meas[k]["pAct"][target_idx]) for k in circles[circle_idx]],
                                   Acts=["" if units[target_idx] == "noconv"
                                         else "{:.3g}".format(meas[k]["Act"][target_idx])
                                         for k in circles[circle_idx]],
                                   qualifiers=[meas[k]["qualifiers"][target_idx] for k in circles[circle_idx]],
                                   nonadd="{:.2g}".format(nonadd),
                                   target=target,
                                   mcss_tot=mcss_tot,
                                   image_file=image_file)

                circle_id = "_".join(list(circles[circle_idx])+[target])
                line = line + circle_id + "\t"
                outlines.append(line)
                for i in range(4):
                    c2c_lines.append(circle_id+"\t"+circles[circle_idx][i]+"\n")

            if len(outlines) < 3:
                continue

            ###########
            # Calculate Theoretical Quantiles and sort accordingly
            if series_column:
                theo_quantiles = []
                for series in set(series_id):
                    a = sorted([(add_diff, x[0])
                                for add_diff, x in zip(add_diffs, enumerate(series_id))
                                if x[1] == series])
                    outline_idxs = [j[1] for j in a]
                    if len(outline_idxs) > 2:
                        quantiles = list(stats.probplot(outline_idxs)[0][0])
                    elif len(outline_idxs) == 2:                       # To avoid a scipy warning
                        quantiles = [-0.54495214,  0.54495214]
                    else:
                        quantiles = [0.]
                    theo_quantiles = theo_quantiles + list(zip(outline_idxs, quantiles))

                    # Write estimated SD per series to STDOUT.
                    # SD is multiplied by 1/sqrt(4) = 0.5 to obtain uncertainty of individual measurement
                    series_add_diffs = [j[0] for j in a]
                    if len(series_add_diffs) > 3:
                        print("and series ", series, " based on ", str(len(series_add_diffs)), " cycles.")
                        print("{0:.2f}".format(0.5 * numpy_std(series_add_diffs)), " from normal SD")
                        print("{0:.2f}".format(0.5 * mad_std(series_add_diffs)), " from MAD")
                        print("{0:.2f}".format(0.5 * Sn_MedMed_std(series_add_diffs)), " from Median of Medians")

                theo_quantiles = [j[1] for j in sorted(theo_quantiles)]

            else:
                a = sorted(zip(add_diffs, range(len(outlines))))
                outline_idxs = [j[1] for j in a]
                theo_quantiles = stats.probplot(sorted(add_diffs))[0][0]
                theo_quantiles = [j[1] for j in sorted(zip(outline_idxs, list(theo_quantiles)))]

                series_add_diffs = [j[0] for j in a]
                if len(series_add_diffs) > 3:
                    print("based on ", str(len(series_add_diffs)), " cycles.")
                    print("{0:.2f}".format(0.5 * numpy_std(series_add_diffs)), " from normal SD")
                    print("{0:.2f}".format(0.5 * mad_std(series_add_diffs)), " from MAD")
                    print("{0:.2f}".format(0.5 * Sn_MedMed_std(series_add_diffs)), " from Median of Medians")

            print()

            for outline, theo_quantile in zip(outlines, theo_quantiles):
                f.write(outline + "{:4.3}".format(theo_quantile) + "\n")

            for ID, values in nonadd_percompound.items():
                if series_column:
                    if len(values) == 0:
                        continue
                    series = meas[ID]["series"]
                    outline = [ID,
                               meas[ID]["smiles"],
                               series,
                               target,
                               meas[ID]["qualifiers"][target_idx],
                               "{:.2g}".format(meas[ID]["pAct"][target_idx])]
                    values_pure = [value[0] for value in values if value[1] == "pure"]
                    if len(values_pure) > 0:
                        nonadd_cnt_pure = len(values_pure)
                        nonadd_pure = sum(values_pure) / nonadd_cnt_pure
                        nonadd_sd_pure = math.sqrt(sum([(i-nonadd_pure)**2 for i in values_pure]) / nonadd_cnt_pure)
                        outline = outline + ["{:.2g}".format(nonadd_pure),
                                             str(nonadd_cnt_pure),
                                             "{:.2g}".format(nonadd_sd_pure)]
                    else:
                        outline = outline + ["", "", ""]

                    values_mixed = [value[0] for value in values if value[1] == "mixed"]
                    if len(values_mixed) > 0:
                        nonadd_cnt_mixed = len(values_mixed)
                        nonadd_mixed = sum(values_mixed) / nonadd_cnt_mixed
                        nonadd_sd_mixed = math.sqrt(sum([(i-nonadd_mixed)**2 for i in values_mixed]) / nonadd_cnt_mixed)
                        outline = outline + ["{:.2g}".format(nonadd_mixed),
                                             str(nonadd_cnt_mixed),
                                             "{:.2g}".format(nonadd_sd_mixed)]
                    else:
                        outline = outline + ["", "", ""]

                    g.write("\t".join(outline)+"\n")

                else:
                    if len(values) == 0:
                        continue
                    nonadd_cnt = len(values)
                    nonadd = sum(values) / nonadd_cnt
                    nonadd_sd = math.sqrt(sum([(i-nonadd)**2 for i in values]) / nonadd_cnt)
                    outline = "\t".join([ID,
                                         meas[ID]["smiles"],
                                         "",
                                         target,
                                         meas[ID]["qualifiers"][target_idx],
                                         "{:.2g}".format(meas[ID]["pAct"][target_idx]),
                                         "{:.2g}".format(nonadd),
                                         str(nonadd_cnt),
                                         "{:.2g}".format(nonadd_sd)])
                    g.write(outline+"\n")

            for line in c2c_lines:
                h.write(line)

    return


def draw_image(ids, smiles, tsmarts, pActs, Acts, qualifiers, nonadd, target, mcss_tot, image_file):
    """
    Draw Nonadditivity Circle to Image file
    """
    
    cpds = [Chem.MolFromSmiles(i) for i in smiles]

    #########
    # Compute Coordinates of local MCSS, aligned with global MCSS
    mcss_loc = Chem.MolFromSmarts(rdFMCS.FindMCS(cpds, completeRingsOnly=True, timeout=10).smartsString)
    
    if mcss_tot:
        mcss_tot_coords = [mcss_tot.GetConformer().GetAtomPosition(x) for x in range(mcss_tot.GetNumAtoms())]
        coords2D_tot = [Geometry.Point2D(pt.x, pt.y) for pt in mcss_tot_coords]
    
        mcss_match = mcss_loc.GetSubstructMatch(mcss_tot)
    
        coordDict = {}
        for i, coord in enumerate(coords2D_tot):
            coordDict[mcss_match[i]] = coord
        
        AllChem.Compute2DCoords(mcss_loc, coordMap=coordDict)
    else:
        AllChem.Compute2DCoords(mcss_loc)
    
    #########
    # Align circle compounds to local MCSS
    
    matchVs = [x.GetSubstructMatch(mcss_loc) for x in cpds]

    # compute reference coordinates:
    mcss_loc_coords = [mcss_loc.GetConformer().GetAtomPosition(x) for x in range(mcss_loc.GetNumAtoms())]
    coords2D_loc = [Geometry.Point2D(pt.x, pt.y) for pt in mcss_loc_coords]

    # generate coords for the other molecules using the common substructure
    for molIdx in range(4):
        mol = cpds[molIdx]
        coordDict = {}
        for i, coord in enumerate(coords2D_loc):
            coordDict[matchVs[molIdx][i]] = coord
        AllChem.Compute2DCoords(mol, coordMap=coordDict)

    ##########
    # Assemble Image

    qualifiers_inv = ["" for i in range(4)]
    for i in range(4):
        if qualifiers[i] == ">":
            qualifiers_inv[i] = "<"
        elif qualifiers[i] == "<":
            qualifiers_inv[i] = ">"
        else:
            continue

    new_im = Image.new("RGB", size=(650, 670), color=(255, 255, 255, 0))
    if Acts[0] != "":
        new_im.paste(Draw.MolToImage(cpds[0],
                                     size=(300, 300),
                                     legend=ids[0] + "        " + qualifiers_inv[0] + Acts[0] + " ("
                                     + qualifiers[0] + pActs[0] + ")"),
                     (0, 0))
        new_im.paste(Draw.MolToImage(cpds[1],
                                     size=(300, 300),
                                     legend=ids[1] + "        " + qualifiers_inv[1] + Acts[1] + " ("
                                     + qualifiers[1] + pActs[1] + ")"),
                     (350, 0))
        new_im.paste(Draw.MolToImage(cpds[2],
                                     size=(300, 300),
                                     legend=ids[2] + "        " + qualifiers_inv[2] + Acts[2] + " ("
                                     + qualifiers[2] + pActs[2]+")"),
                     (350, 350))
        new_im.paste(Draw.MolToImage(cpds[3],
                                     size=(300, 300),
                                     legend=ids[3] + "        " + qualifiers_inv[3] + Acts[3] + " ("
                                     + qualifiers[3] + pActs[3] + ")"),
                     (0, 350))

        draw = ImageDraw.Draw(new_im)
        font = ImageFont.truetype(font_path, 14)
        draw.text((260, 330), "Nonadditivity: " + nonadd, fill=(0, 0, 0, 0), font=font)

        font = ImageFont.truetype(font_path, 10)
        draw.text((10, 650), "[uM]  (-log10[M])  Activity in Assay: " + target, fill=(0, 0, 0, 0), font=font)
    else:
        new_im.paste(Draw.MolToImage(cpds[0],
                                     size=(300, 300),
                                     legend=ids[0]+"        "+qualifiers[0]+pActs[0]),
                     (0, 0))
        new_im.paste(Draw.MolToImage(cpds[1],
                                     size=(300, 300),
                                     legend=ids[1]+"        "+qualifiers[1]+pActs[1]),
                     (350, 0))
        new_im.paste(Draw.MolToImage(cpds[2],
                                     size=(300, 300),
                                     legend=ids[2]+"        "+qualifiers[2]+pActs[2]),
                     (350, 350))
        new_im.paste(Draw.MolToImage(cpds[3],
                                     size=(300, 300),
                                     legend=ids[3]+"        "+qualifiers[3]+pActs[3]),
                     (0, 350))

        draw = ImageDraw.Draw(new_im)
        font = ImageFont.truetype(font_path, 14)
        draw.text((260, 330), "Nonadditivity: " + nonadd, fill=(0, 0, 0, 0), font=font)

        font = ImageFont.truetype(font_path, 10)
        draw.text((10, 650), "Activity in Assay: " + target, fill=(0, 0, 0, 0), font=font)

    # Draw Arrows
    draw.line((300, 150, 350, 150), fill=0, width=2)
    draw.line((340, 145, 350, 150), fill=0, width=2)
    draw.line((340, 155, 350, 150), fill=0, width=2)

    draw.line((300, 500, 350, 500), fill=0, width=2)
    draw.line((340, 495, 350, 500), fill=0, width=2)
    draw.line((340, 505, 350, 500), fill=0, width=2)

    draw.line((150, 300, 150, 350), fill=0, width=2)
    draw.line((145, 340, 150, 350), fill=0, width=2)
    draw.line((155, 340, 150, 350), fill=0, width=2)

    draw.line((500, 300, 500, 350), fill=0, width=2)
    draw.line((495, 340, 500, 350), fill=0, width=2)
    draw.line((505, 340, 500, 350), fill=0, width=2)

    # Add Reaction Parts
    b = Chem.MolFromSmiles(tsmarts[0][:tsmarts[0].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 90))

    b = Chem.MolFromSmiles(tsmarts[0][tsmarts[0].index(">")+2:])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 160))

    b = Chem.MolFromSmiles(tsmarts[0][:tsmarts[0].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 440))

    b = Chem.MolFromSmiles(tsmarts[0][tsmarts[0].index(">")+2:])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 510))

    b = Chem.MolFromSmiles(tsmarts[1][:tsmarts[1].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (80, 300))

    b = Chem.MolFromSmiles(tsmarts[1][tsmarts[1].index(">")+2:])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (170, 300))

    b = Chem.MolFromSmiles(tsmarts[1][:tsmarts[1].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (430, 300))

    b = Chem.MolFromSmiles(tsmarts[1][tsmarts[1].index(">")+2:])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (520, 300))

    new_im.save(image_file)

    return


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def write_smiles_id_file(meas, infile, max_heavy):
    """
    Format dataset for MMP analysis and write to temp file
    """

    smifile = infile[:infile.index(".")]+"_ligands.smi"
    f = open(smifile, 'w')

    for entry in meas.keys():
        try:
            mol = Chem.MolFromSmiles(meas[entry]["smiles"])
            if mol.GetNumHeavyAtoms() > max_heavy:
                continue
            f.write(meas[entry]["smiles"] + "\t" + entry + '\n')
        except:
            print("Skipping compound", entry)
            continue

    f.close()

    return


def calc_raw_MMPs(infile, update):
    """
Generate MMP Indexing and Matching using mmpdb
    """
    smifile = infile[:infile.index(".")] + "_ligands.smi"
    fragfile = infile[:infile.index(".")] + ".fragments"
    mmp_outfile = infile[:infile.index(".")] + "_mmp_raw.csv"

    # TODO switch system calls to just importing the python code and using it directly
    if update:
        print("Updating MMP Fragments for " + infile)
        sp_call = 'python -m mmpdblib fragment ' \
                  + ' --num-jobs 20 --delimiter tab --cache ' \
                  + fragfile + ' --output ' + fragfile + ' ' + smifile
    else:
        print("Generating MMP Fragments for " + infile)
        sp_call = 'python -m mmpdblib fragment' \
                  + ' --num-jobs 20 -i smi --delimiter tab ' \
                  + ' --max-rotatable-bonds 20 --output ' + fragfile + ' ' + smifile

    subprocess.call(sp_call, shell=True)        # Fragmentation

    print("Indexing MMP Fragments for " + infile)
    sp_call = 'python -m mmpdblib index'
    sp_call = sp_call + ' --out csv --symmetric --output ' + mmp_outfile + ' ' + fragfile
    subprocess.call(sp_call, shell=True)        # Fragment_matching & Indexing

    return


def read_raw_mmps(infile):
    """
Read raw precalculated MMPs
    """
    mmp_outfile = infile[:infile.index(".")]+"_mmp_raw.csv"
    f = open(mmp_outfile, 'r')
    mmps = f.readlines()
    f.close()

    return mmps


def build_ligand_dictionary_from_infile(infile, props, units, delimiter, series_column):
    """
    Read input file and assemble dictionaries
    """

    error_files = infile[:infile.index(".")] + "_problem_smiles.smi"

    if delimiter == "tab":
        delimiter = "\t"
    elif delimiter == "space":
        delimiter = " "
    elif delimiter == "semicolon":
        delimiter = ";"
    else:
        delimiter = ","

    with open(infile, 'r') as f, open(error_files, 'w') as g:
        ########
        # Process header
        header = [i.strip('"') for i in f.readline().rstrip('\n').split(delimiter)]

        ########
        # Figure out Column ID of SMILES and ID column
        id_col = [i for i, name in enumerate(header) if "SRN" in name or "ID" in name]
        id_col = (0 if len(id_col) == 0 else id_col[0])

        smi_col = [i for i, name in enumerate(header) if "smiles" in name.lower()]
        smi_col = (1 if len(smi_col) == 0 else smi_col[0])

        if series_column:
            ser_col = header.index(series_column)

        ########
        # Figure out target column ids
        if not props:
            act_col = [2]
            props = [header[2]]
        else:
            try:
                act_col = [header.index(i) for i in props]
            except:
                print("Could not find all given Activity columns in file header. Exiting")
                exit(0)

        ########
        # Figure out conversion of target columns
        # Valid Flags for not converting activity data to pActivity: pIC50, pEC50, pKi, pKd, noconv
        log_flags = ["pIC50", "pEC50", "pKi", "pKd", "pCC50", "pIC20", "pID50", "noconv"]
        col_convert = [False if any(log_flag.lower() in target.lower() for log_flag in log_flags)
                       else True for target in props]
        log10 = [False for target in props]

        #########
        # Write Identified Columns to STDOUT
        print("Identifier Column found: " + header[id_col])
        print("Smiles column found: " + header[smi_col])
        for i in range(len(props)):
            if len(units) > 0:
                if units[i] == "noconv":
                    col_convert[i] = False
                elif units[i] == "log10":
                    col_convert[i] = False
                    log10[i] = True
            if col_convert[i]:
                print("Activity column #" + str(i + 1) + ": " + props[i] + " will be converted to -log10(" + props[
                    i] + ")")
            elif log10[i]:
                print("Activity column #" + str(i + 1) + ": " + props[i] + " will be converted to log10(" + props[
                    i] + ")")
            else:
                print("Activity column #" + str(i + 1) + ": " + props[i])

        if series_column:
            print("Series Column found: " + header[ser_col])

        if id_col == smi_col or id_col in act_col or smi_col in act_col:
            print("Was not able to cleanly distinguish ID, SMILES, and activity columns.")
            print("Please assign unambiguous names (no overlap in 'SMILES', 'ID', 'SRN'.")
            print("Exiting.")
            exit(0)

            ########
            # Assemble data
        remover = SaltRemover.SaltRemover(defnFilename=salt_defns)
        unit_conv = {"M": 1.0,
                     "mM": 1E-3,
                     "uM": 1E-6,
                     "nM": 1E-9,
                     "pM": 1E-12,
                     "noconv": 1.0}
        meas = dict()
        smiles_registered = dict()
 
        for line in f:
            line = [i.strip('"') for i in line.rstrip('\n').split(delimiter)]
            if line[0][0] == "#":         # skip commented-out compounds
                continue
            compound_id = line[id_col]
            if compound_id in meas:
                print("Two or more entries for the same identifier: " + compound_id)
                print("Please fix. Exiting")
                exit(0)
            smiles = line[smi_col].replace("\\\\", "\\")
            if len(line) < len(props) + 2:
                print("Could not properly read line:")
                print(line)
                print("Exiting...")
                exit(0)
            try:
                mol = Chem.MolFromSmiles(smiles)
                res = remover.StripMol(mol)  # Remove Salts
                smiles = Chem.MolToSmiles(res, True)  # Canonicalize smiles
                mwt = Descriptors.MolWt(mol)
                if "." in smiles:
                    print("Found unknown salt in " + line[id_col] + ": " + smiles)
                    print("This compound will be ignored.")
                    continue
            except:
                print("Could not properly read SMILES " + smiles + "(see error SMILES file)")
                print("This compound will be ignored.")
                g.write(smiles + "\n")
                continue

            if smiles in smiles_registered:
                print("Two entries with the same structure: "+smiles_registered[smiles]+" and "+compound_id)
                print("Nonadd will use the first compound and discard the second.\n")
                continue
            else:
                smiles_registered[smiles] = compound_id

            meas[compound_id] = dict(smiles=smiles, Act=[], pAct=[], qualifiers=[], mwt=mwt, series=None)
            if series_column:
                meas[compound_id]["series"] = line[ser_col]
            for i, target in enumerate(props):
                if col_convert[i]:
                    u_conv = unit_conv["uM"]
                    if not len(units) == 0:
                        try:
                            u_conv = unit_conv[units[i]]
                        except:
                            print("Given unit " + units[i] + " has not been recognized.")
                            print("Please give one out of [M, mM, uM, nM, pM, noconv]")
                    if line[act_col[i]] in ["NA", "", "No Value"]:
                        meas[compound_id]["qualifiers"].append("")
                        meas[compound_id]["Act"].append("NA")
                        meas[compound_id]["pAct"].append("NA")
                    elif is_number(line[act_col[i]]):
                        if float(line[act_col[i]]) <= 0.0:
                            print("Cannot interpret measured activity of "+line[act_col[i]]+units[i]+" for compound "+compound_id)
                            print("Please fix. Exiting")
                            exit(0)
                        meas[compound_id]["qualifiers"].append("")
                        meas[compound_id]["Act"].append(float(line[act_col[i]]))
                        meas[compound_id]["pAct"].append((-1) * math.log10(float(line[act_col[i]]) * u_conv))
                    elif line[act_col[i]][0] in (">", "<", "*") and is_number(line[act_col[i]][1:]):
                        meas[compound_id]["qualifiers"].append(line[act_col[i]][0])
                        meas[compound_id]["Act"].append(float(line[act_col[i]][1:]))
                        meas[compound_id]["pAct"].append((-1) * math.log10(float(line[act_col[i]][1:]) * u_conv))
                    else:
                        print("Did not recognize number " + str(line[act_col[i]]))
                        print(" in line: " + " ".join(line))
                        print("Please fix. Exiting")
                        exit(0)
                elif log10[i]:
                    if line[act_col[i]] in ["NA", "", "No Value"]:
                        meas[compound_id]["qualifiers"].append("")
                        meas[compound_id]["Act"].append("NA")
                        meas[compound_id]["pAct"].append("NA")
                    elif is_number(line[act_col[i]]):
                        if float(line[act_col[i]]) <= 0.0:
                            print("Cannot interpret measured activity of "+line[act_col[i]]+units[i]+" for compound "+compound_id)
                            print("Please fix. Exiting")
                            exit(0)
                        meas[compound_id]["qualifiers"].append("")
                        meas[compound_id]["Act"].append(float(line[act_col[i]]))
                        meas[compound_id]["pAct"].append(math.log10(float(line[act_col[i]])))
                    elif line[act_col[i]][0] in (">", "<", "*") and is_number(line[act_col[i]][1:]):
                        meas[compound_id]["qualifiers"].append(line[act_col[i]][0])
                        meas[compound_id]["Act"].append(float(line[act_col[i]][1:]))
                        meas[compound_id]["pAct"].append((-1) * math.log10(float(line[act_col[i]][1:])))
                    else:
                        print("Did not recognize number " + str(line[act_col[i]]))
                        print(" in line: " + " ".join(line))
                        print("Please fix. Exiting")
                        exit(0)
                else:
                    if line[act_col[i]] in ["NA", "", "No Value"]:
                        meas[compound_id]["qualifiers"].append("")
                        meas[compound_id]["Act"].append("NA")
                        meas[compound_id]["pAct"].append("NA")
                    elif is_number(line[act_col[i]]):
                        meas[line[id_col]]["qualifiers"].append("")
                        if len(units) > 0:
                            if units[i] == "noconv":
                                meas[compound_id]["Act"].append("")
                            else:
                                meas[compound_id]["Act"].append(1E6 * 10 ** ((-1) * float(line[act_col[i]])))
                        else:
                            meas[compound_id]["Act"].append(1E6 * 10 ** ((-1) * float(line[act_col[i]])))
                        meas[compound_id]["pAct"].append(float(line[act_col[i]]))
                    elif line[act_col[i]][0] in (">", "<", "*") and is_number(line[act_col[i]][1:]):
                        meas[compound_id]["qualifiers"].append(line[act_col[i]][0])
                        if len(units) > 0:
                            if units[i] == "noconv":
                                meas[compound_id]["Act"].append("")
                            else:
                                meas[compound_id]["Act"].append(1E6 * 10 ** ((-1) * float(line[act_col[i]][1:])))
                        else:
                            meas[compound_id]["Act"].append(1E6 * 10 ** ((-1) * float(line[act_col[i]][1:])))
                        meas[compound_id]["pAct"].append(float(line[act_col[i]][1:]))
                    else:
                        print("Did not recognize number " + str(line[act_col[i]]))
                        print(" in line: " + " ".join(line))
                        print("Please fix. Exiting")
                        exit(0)

    if len(units) == 0:
        units = ["noconv" for i in props]

    return meas, props, units


def clean_image_folder(meas, props, infile):
    """
    Delete Cycle images where the properties have changed
    """
    print("Deleting Cycle images where property values changed since last run.")
    
    image_files = os.listdir("images")

    try:
        with open(infile[:infile.index(".")] + "_image_dat.pkl", "rb") as oldprops_file:
            oldprops = pickle.load(oldprops_file)
        remove_prop = {}
        for cid, a in meas.items():
            for idx, prp in enumerate(props):
                if (cid, prp) in oldprops:
                    if not "{:.3g}".format(a['pAct'][idx]) == oldprops[(cid, prp)]:
                        remove_prop[(cid, prp)] = True

        del_file = [False for i in image_files]
        for cid, prp in iter(remove_prop.keys()):
            for idx, fl in enumerate(image_files):
                if cid in fl and prp in fl:
                    del_file[idx] = True
    except:
        print("Could not find the file with properties from the previous run ("
              + infile[:infile.index(".")] + "_image_dat.pkl"+").")
        print("All pictures will be redone.")
        del_file = [True for i in image_files]

    print("Reusing " + str(len(del_file)-sum(del_file)) + " out of " + str(len(del_file))+" previous images.")
    
    for idx, delete in enumerate(del_file):
        if delete:
            os.remove("images/"+image_files[idx])

    return


def write_propdata_dict(meas, props, infile):
    """
    Store property values for future image update
    """

    dat = {}
    for cid, a in meas.items():
        for idx, prp in enumerate(props):
            try:
                dat[(cid, prp)] = "{:.3g}".format(a['pAct'][idx])
            except:
                continue
                
    with open(infile[:infile.index(".")]+"_image_dat.pkl", "wb") as output_file:
        pickle.dump(dat, output_file)
        
    return


def run_nonadd_calculation(run_control):
    """
    Main routine to run the Nonadditivity calculations.
    run_control needs to have the following attributes:
    - infile (required)
    - outfile (optional)
    - props (can be empty)
    - units (as many as props or empty)
    - shorts (as many as props or empty)
    - max_heavy (required)
    - update (Boolean)
    - no_chiral(Boolean)
    - write_images(Boolean)
    """
    meas, run_control.props, run_control.units = build_ligand_dictionary_from_infile(run_control.infile,
                                                                                     run_control.props,
                                                                                     run_control.units,
                                                                                     run_control.delimiter,
                                                                                     run_control.series_column)

    if run_control.update and not os.path.exists(run_control.infile[:run_control.infile.index(".")]+".fragments"):
        print("Was not able to locate results from previous fragmentation.")
        print("Will redo all fragmentation.")
        run_control.update = False

    write_smiles_id_file(meas, run_control.infile, run_control.max_heavy)
    calc_raw_MMPs(run_control.infile, run_control.update)

    mmps = read_raw_mmps(run_control.infile)
    neighs = build_neighbor_dictionary(mmps, no_chiral=run_control.no_chiral)
    circles = get_circles(neighs)

    if not run_control.outfile:
        if os.path.dirname(run_control.infile) == '':
            run_control.outfile = "Additivity_diffs_" + run_control.infile[:run_control.infile.index(".")]+".txt"
        else:
            run_control.outfile = os.path.dirname(run_control.infile) + '/' + "Additivity_diffs_"
            run_control.outfile = run_control.outfile + run_control.infile.split("/")[-1][:run_control.infile.split("/")[-1].index(".")]+".txt"
    else:
        if "." not in run_control.outfile:
            run_control.outfile = run_control.outfile +".txt"

    if run_control.shorts:
        if len(run_control.shorts) == len(run_control.props):
            run_control.props = run_control.shorts

    if run_control.update and run_control.write_images:
        if os.path.isdir("images"):
            clean_image_folder(meas, run_control.props, run_control.infile)

    write_propdata_dict(meas, run_control.props, run_control.infile)
    write_circles_to_output(circles=circles,
                            meas=meas,
                            neighs=neighs,
                            outfile=run_control.outfile,
                            props=run_control.props,
                            units=run_control.units,
                            images=run_control.write_images,
                            include_censored=run_control.include_censored,
                            update=run_control.update,
                            series_column=run_control.series_column)

    exit(0)
