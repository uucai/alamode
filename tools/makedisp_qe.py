import sys
import os
import collections
import numpy as np
import argparse
from pymatgen import Structure
from pymatgen.io.vasp import inputs
from pymatgen.core.periodic_table import get_el_sp
from interface.QE import QEParser
from GenDisplacement import AlamodeDisplace

ALAMODE_root = "~/Work/alamode"

parser = argparse.ArgumentParser()

parser.add_argument('--mag',
                    type=float, default=0.02,
                    help="Magnitude of displacement in units of \
                        Angstrom (default: 0.02)")

parser.add_argument('--prefix',
                    type=str, default="disp",
                    help="Prefix of the files to be created. (default: disp)")

parser.add_argument('-d', '--dim',
                    default=None, type=str,
                    help="Transformation matrix")

parser.add_argument('file_primitive', metavar='primitive.pw.in',
                    default=None,
                    help="Original primitive cell input file for pw.x")



def gen_kpoints_file(structure):

    kmesh = inputs.Kpoints.automatic_density_by_vol(structure, 450)
    kmesh_dict = kmesh.as_dict()

    return kmesh_dict['kpoints'], kmesh_dict['usershift']


def gen_species_dictionary(atomic_number_uniq):

    species_dict = {}
    counter = 1
    for num in atomic_number_uniq:
        species_dict[num] = counter
        counter += 1
    return species_dict


def gen_alm_input(filename, prefix, mode, structure, norder, str_cutoff,
                  ndata=0, dfset='DFSET'):

    Bohr = 0.52917721067

    if (mode != "suggest" and mode != "optimize"):
        raise RuntimeError("Invalid MODE: %s" % mode)

    atomic_numbers_uniq = list(
        collections.OrderedDict.fromkeys(structure.atomic_numbers))

    species_index = gen_species_dictionary(atomic_numbers_uniq)
    # Make input for ALM
    with open(filename, 'w') as f:
        f.write("&general\n")
        f.write(" PREFIX = %s\n" % prefix)
        f.write(" MODE = %s\n" % mode)
        f.write(" NAT = %i\n" % structure.num_sites)
        str_spec = ""
        for num in atomic_numbers_uniq:
            str_spec += str(get_el_sp(num)) + " "
        f.write(" NKD = %i; KD = %s\n" % (structure.ntypesp, str_spec))
        f.write(" TOLERANCE = 1.0e-3\n")
        f.write("/\n\n")
        f.write("&interaction\n")
        f.write(" NORDER = %i\n" % norder)
        f.write("/\n\n")
        f.write("&cutoff\n")
        f.write(" %s\n" % str_cutoff)
        f.write("/\n\n")
        f.write("&cell\n")
        f.write("%20.14f\n" % (1.0/Bohr))
        for i in range(3):
            for j in range(3):
                f.write("%20.13f" % structure.lattice.matrix[i][j])
            f.write("\n")
        f.write("/\n\n")
        f.write("&position\n")
        for i in range(len(structure.frac_coords)):
            f.write("%4i " % species_index[structure.atomic_numbers[i]])
            f.write("%20.14f " % structure.frac_coords[i][0])
            f.write("%20.14f " % structure.frac_coords[i][1])
            f.write("%20.14f\n" % structure.frac_coords[i][2])
        f.write("/\n\n")

        if mode == "optimize":
            f.write("&optimize\n")
            f.write(" DFSET = %s\n" % dfset)
            f.write("/\n\n")

def process_args(args):

    if args.dim is None:
        scaling_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    else:
        str_dim = args.dim.strip().split()
        if len(str_dim) == 1:
            nsize = int(str_dim[0])
            scaling_matrix = [[nsize, 0, 0], [0, nsize, 0], [0, 0, nsize]]
        elif len(str_dim) == 3:
            nsize1 = int(str_dim[0])
            nsize2 = int(str_dim[1])
            nsize3 = int(str_dim[2])
            scaling_matrix = [[nsize1, 0, 0], [0, nsize2, 0], [0, 0, nsize3]]
        elif len(str_dim) == 9:
            nsizes = [int(t) for t in str_dim]
            scaling_matrix = [[nsizes[0], nsizes[1], nsizes[2]],
                              [nsizes[3], nsizes[4], nsizes[5]],
                              [nsizes[6], nsizes[7], nsizes[8]]]
        else:
            raise RuntimeError("Invalid format of --dim.")

    return args.mag, args.prefix, scaling_matrix


def update_qeobj(qeparse_in, structure_in):

    kmesh, kshift = gen_kpoints_file(structure)

    # Update nat entry
    list_system_new = qeparse_in.list_system
    for i, elems in enumerate(list_system_new):
        if 'nat' in elems:
            list_system_new[i] = '    nat = %d\n' % structure_in.num_sites

    # Update Cell parameters
    list_CELL_PARAMETERS_new = []
    list_CELL_PARAMETERS_new.append(qeparse_in.list_cell_parameters[0])
    for i in range(3):
        str_tmp = ""
        for j in range(3):
            str_tmp += str("%20.15f" % structure_in.lattice.matrix[i][j])
        str_tmp += '\n'

        list_CELL_PARAMETERS_new.append(str_tmp)

    # Update Kpoint
    list_kmesh = []
    list_kmesh.extend([str(i) for i in kmesh[0]])
    list_kmesh.extend([str(i) for i in kshift])
    str_kmesh = ""
    for entry in list_kmesh:
        str_kmesh += "%3s" % entry

    str_kmesh += '\n'
    list_K_POINTS_new = []
    list_K_POINTS_new.append(qeparse_in.list_k_points[0])
    list_K_POINTS_new.append(str_kmesh)

    qeparse_in.lattice_vector = structure_in.lattice.matrix.transpose()
    qeparse_in.x_fractional = structure_in.frac_coords
    qeparse_in.kd_in_str = [str(t) for t in structure_in.species]
    qeparse_in.nat = structure_in.num_sites
    qeparse_in.list_system = list_system_new
    qeparse_in.list_cell_parameters = list_CELL_PARAMETERS_new
    qeparse_in.list_k_points = list_K_POINTS_new

    return qeparse_in


if __name__ == '__main__':

    args = parser.parse_args()
    disp_magnitude_angstrom, prefix, scaling_matrix = process_args(args)

    qeobj = QEParser()
    qeobj.load_initial_structure(args.file_primitive)

    structure = Structure(qeobj.lattice_vector.transpose(),
                          qeobj.kd_in_str,
                          qeobj.x_fractional)

    Structure.make_supercell(structure, scaling_matrix)

    print("Supercell generated. # Atoms: %i" % structure.num_sites)
    print("")

    suffix = 'scf.in'
    prefix0 = 'supercell'
    disp = np.zeros((structure.num_sites, 3))

    # update structural information of qeobj
    qeobj_mod = update_qeobj(qeobj, structure)

    # create the supercell structure
    qeobj_mod.generate_structures(prefix0, ['original'], [disp])
    # rename the file
    command = ("mv %s1.pw.in %s0.scf.in" % (prefix0, prefix0))
    os.system(command)

    # Generate displacement files
    gen_alm_input('ALM0.in', prefix0, 'suggest', structure, 1, "*-* None")
    command = ("%s/alm/alm ALM0.in > ALM0.log" % ALAMODE_root)
    os.system(command)

    dispobj = AlamodeDisplace("fd", qeobj_mod, verbosity=0)
    header_list, disp_list \
        = dispobj.generate(file_pattern=["%s.pattern_HARMONIC" % prefix0],
                           magnitude=disp_magnitude_angstrom)

    qeobj_mod.generate_structures(prefix,
                                  header_list,
                                  disp_list)

