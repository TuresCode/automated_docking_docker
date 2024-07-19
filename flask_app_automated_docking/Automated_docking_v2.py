# Automated docking script for docking with AutoDock Vina
# Author: Fabian Meyer
# Open Babel 3.0.0 -- Mar 11 2020 -- 11:52:36

# The following packages need to be installed:
# pip install meeko==0.1.dev2
# conda install -c conda-forge openbabel

import argparse
from rdkit import Chem  # conda create -c conda-forge -n my-rdkit-env rdkit
from rdkit.Chem.PandasTools import LoadSDF
from rdkit.Chem import Descriptors

import os
import sys
import pandas as pd
import subprocess

from biopandas.pdb import PandasPdb
from biopandas.mol2 import split_multimol2
from biopandas.mol2 import PandasMol2

# pdmol = PandasMol2()
import shutil
import pickle

import multiprocessing
import matplotlib.pyplot as plt
from natsort import natsorted

import glob
import time
import ast

from tqdm import tqdm
from biopandas.mol2 import PandasMol2

import numpy

# conda install -c conda-forge -c schrodinger pymol-bundle python needs to be 3.8.10 or lower

from pymol.cgo import *
from pymol import cmd


from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, QED
from rdkit.Chem.Draw import IPythonConsole
import numpy as np

# We suppres stdout from invalid smiles and validations
from rdkit import rdBase

rdBase.DisableLog("rdApp.*")
Chem.PandasTools.RenderImagesInAllDataFrames(images=True)
pd.options.mode.chained_assignment = None  # supress warnings


class config:
    def __init__(
        self,
        working_dir,
        center=[0, 0, 0],
        cofator="",
        plot=False,
        gpu_vina=False,
        result_path="Results",
        NADP_cofactor=False,
        metal_containing=False,
        align=False,
        output_formate="pdbqt",
        show_plots=False,
        fixed_torsion=False,
        size=[15, 15, 15],
        num_modes=9,
        exhaustiveness=64,
        energy_range=3,
    ):
        self.working_dir = working_dir
        self.center = center
        self.cofactor = cofator
        self.plot = plot
        self.gpu_vina = gpu_vina
        self.NADP_cofactor = NADP_cofactor
        self.metal_containing = metal_containing
        self.align = align
        self.output_formate = output_formate
        self.show_plots = show_plots
        self.size = size
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness
        self.energy_range = energy_range
        self.result_path = result_path
        self.fixed_torsion = fixed_torsion


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Automated docking script for docking with AutoDock Vina. Required parameter is the working directory containing the ligand(s) as mol2 file(s) or ligands3D.sdf and the receptor(s) as pdb file(s)"
    )
    parser.add_argument(
        "--working_dir",
        required=True,
        type=str,
        help="working directory containing the ligand(s) as mol2 file(s) and the receptor(s) as pdb file(s)",
    )
    parser.add_argument(
        "--center",
        nargs="+",
        default=[0, 0, 0],
        type=float,
        help="center coordinates for the docking box eg. --center -10.5 25 22",
    )
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--gpu_vina", action="store_true", help="Enable GPU Vina")
    parser.add_argument(
        "--result_path", type=str, default="Results", help="Path to result output"
    )
    parser.add_argument(
        "--NADP_cofactor", action="store_true", help="Enable NADP cofactor"
    )
    parser.add_argument(
        "--metal_containing", action="store_true", help="Enable metal-containing"
    )
    parser.add_argument("--align", action="store_true", help="Enable alignment")
    parser.add_argument(
        "--output_format", type=str, default="pdbqt", help="Output format"
    )
    parser.add_argument("--show_plots", action="store_true", help="Show plots")
    parser.add_argument(
        "--size", nargs="+", default=[15, 15, 15], type=float, help="Box size"
    )
    parser.add_argument("--num_modes", type=int, default=9, help="Number of modes")
    parser.add_argument("--exhaustiveness", type=int, default=64, help="Exhaustiveness")
    parser.add_argument("--energy_range", type=int, default=3, help="Energy range")
    parser.add_argument(
        "--fixed_torsion", action="store_true", help="Enable fixed torsion for ligand"
    )
    parser.add_argument(
        "--cofactor",
        type=str,
        default="",
        help="Cofactor to be docked not implemented yet",
    )

    args = parser.parse_args()

    return config(
        working_dir=args.working_dir,
        center=args.center,
        plot=args.plot,
        gpu_vina=args.gpu_vina,
        result_path=args.result_path,
        NADP_cofactor=args.NADP_cofactor,
        metal_containing=args.metal_containing,
        align=args.align,
        output_formate=args.output_format,
        show_plots=args.show_plots,
        size=args.size,
        num_modes=args.num_modes,
        exhaustiveness=args.exhaustiveness,
        energy_range=args.energy_range,
        fixed_torsion=args.fixed_torsion,
        cofator=args.cofactor,
    )


# Usage example:
parsed_config = parse_arguments()
for key, value in vars(parsed_config).items():
    print(key, value)

config = parsed_config


if config.metal_containing and config.NADP_cofactor:
    print("Please choose only one cofactor")
    sys.exit()
# if both False we need to give an error
if not config.metal_containing and not config.NADP_cofactor:
    if config.center == [0, 0, 0]:
        print("Please choose a cofactor or provide center coordinates")
        sys.exit()

center = [
    0,
    0,
    0,
]  # default center coordinates are overwritten if center coordinates are provided
if config.center != [0, 0, 0]:
    center = config.center

message = "Automated docking workflow. Required parameter include: molecules as smiles, sdf or mol2, pdb of enzymes and box centers or if nadph as cofactor C5N will be used as center"
print(message)

working_dir = (
    config.working_dir
)  # working directory containing the ligand(s) as mol2 file(s) and the receptor(s) as pdb file(s)

try:
    ligand_files = [
        file
        for file in os.listdir(working_dir)
        if file.endswith(".smiles") or file.endswith(".sdf")
    ][
        0
    ]  # file containing all ligands as smiles or sdf file
    print(".smiles or .sdf file detected namely: ", ligand_files)
except:
    ligand_files = ""


slopes_file = ""  # file containing all slopes
excel_sheet = ""  # sheet containing slopes

mol2_files = [
    file for file in os.listdir(working_dir) if file.endswith(".mol2")
]  # file containing all ligands as mol2 file make sure ligand_files is empty obtain the mol2 from Chemdraw 3D
if mol2_files != []:
    print("mol2 files detected: ", mol2_files)
# Set up docking parameters

print(config.gpu_vina, "config.gpu_vina")


def log(log_file_path, message):
    with open(log_file_path, "a+") as f:
        f.write(message + "\n")
    f.close()


def sdf_to_df(substances_db):
    """Function to convert 3d-sdf file to df_ligand"""
    os.system(f"obabel {substances_db} -O test.smiles")
    with open("test.smiles", "r") as f:
        lines = f.readlines()
        molecule_library = {}
        for line in lines:
            key = line.split(sep="\t")[1].strip("\n")
            item = line.split(sep="\t")[0]
            molecule_library[key] = item

    df_ligands = pd.DataFrame(columns=["ligand_names", "molecule_smiles"])
    for key, value in molecule_library.items():
        df_ligands = df_ligands.append(
            {"ligand_names": key, "molecule_smiles": value}, ignore_index=True
        )

    return df_ligands


# if list with smiles. smiles are converted into dictionary with smiles as key and molecule as value. In case of mol2 files, get all mol2 files in folder and convert them into smiles and create df_ligands
def mol2_to_df(mol2_files):
    """Function to convert mol2 files to df_ligand and create 3-dimensional sdf files"""
    # remove files with .smiles extension
    if os.path.exists("ligands.smiles"):
        print("removing old ligands.smiles")
        os.remove("ligands.smiles")
    if os.path.exists("ligands.sdf"):
        print("removing old ligands.sdf")
        os.remove("ligands.sdf")

    smiles_to_combine = ""
    for mol2_file in mol2_files:
        smiles_to_combine += mol2_file.strip(".mol2") + ".smiles "
        obabel_command = f"obabel {mol2_file} -O {mol2_file.strip('.mol2')}.smiles --title {mol2_file.strip('.mol2')}"
        print(obabel_command)
        os.system(obabel_command)

    obabel_command = f"obabel {smiles_to_combine} -O ligands.smiles"
    print(obabel_command)
    os.system(obabel_command)
    # clean up smiles files
    for file in smiles_to_combine.split(sep=" "):
        if file != "":
            os.remove(file)

    f = open("ligands.smiles", "r")
    lines = f.readlines()
    molecule_library = {}
    for line in lines:
        key = line.split(sep="\t")[1].strip("\n")
        item = line.split(sep="\t")[0]
        molecule_library[key] = item

    df_ligands = pd.DataFrame(columns=["ligand_names", "molecule_smiles"])
    for key, value in molecule_library.items():
        df_ligands = df_ligands.append(
            {"ligand_names": key, "molecule_smiles": value}, ignore_index=True
        )

        # generate molecules as sdf file
    mols = [Chem.MolFromSmiles(x) for x in df_ligands["molecule_smiles"]]
    with Chem.SDWriter("ligands.sdf") as w:
        for idx, m in enumerate(mols):
            w.write(m)

    if not glob.glob("ligands3D.sdf"):
        generate_3d_sdf = f"obabel ligands.sdf -O ligands3D.sdf --gen3D -aa"
        ps = subprocess.Popen(
            generate_3d_sdf,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        stdout, stderr = ps.communicate()
        print(stdout, stderr)
    return df_ligands


def make_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        m = Chem.MolToInchi(mol)
        mol_cured = Chem.MolFromInchi(m)

    except:
        print("Error with smile: ", smiles)

    return mol_cured


def find_box_center(filename, cent):
    # identify heteroatoms and take first chain_id, identify c4n only possible in case of NAD(maybe) or NADP needs to be adjusted according to pdb

    # enzyme_files = [file for file in os.listdir() if file.endswith('.pdb')]

    enz = PandasPdb().read_pdb(filename)
    try:
        cofactors = enz.df["HETATM"]
        gb = cofactors.groupby("chain_id")
        cofactor1 = [gb.get_group(x) for x in gb.groups][0]
        c5n = cofactor1[cofactor1.atom_name == "C5N"]
        if c5n.empty:
            c5n = cofactor1[cofactor1.atom_name == "C5"]

        ref_point1 = (
            c5n.x_coord.values[0],
            c5n.y_coord.values[0],
            c5n.z_coord.values[0],
        )
        return list(ref_point1)
    except:
        print("could not find cofactor in: ", filename)
        return cent


def prepare_docking_files():
    if os.path.exists(config.working_dir):
        os.chdir(config.working_dir)  # we are now in the docking folder
        print("Docking folder:", os.getcwd())
    else:
        print("no valid docking folder")

    if not os.path.exists("Experiment_log_file.txt"):
        log_file_path = os.path.join(config.working_dir, "Experiment_log_file.txt")
    else:
        log_file_path = os.path.join(config.working_dir, "Experiment_log_file.txt")
        print("Log file already there!")

    log(log_file_path, "***" * 20 + "New Experiment started" + "***" * 20)
    log(log_file_path, "Working directory: " + str(config.working_dir))
    log(log_file_path, message)
    log(log_file_path, "box size in angström: " + str(config.size))
    log(log_file_path, "num_modes: " + str(config.num_modes))
    log(log_file_path, "exhaustiveness: " + str(config.exhaustiveness))
    log(log_file_path, "energy_range: " + str(config.energy_range))
    log(log_file_path, "output_formate: " + config.output_formate)
    log(
        log_file_path,
        "center_coordinates only relevant if NADP_cofactor or metal_containing is False: "
        + str(config.center),
    )
    log(log_file_path, "metal containing: " + str(config.metal_containing))
    log(log_file_path, "NADP_cofactor: " + str(config.NADP_cofactor))

    log(log_file_path, "GPU: " + str(config.gpu_vina))


def receptor_list(formate):
    if formate == ".pdbqt":
        path = config.working_dir + "/pdbqt_folder"
    else:
        path = config.working_dir
    all_receptors = []
    for file in os.listdir(path):
        if file.endswith(formate):
            all_receptors.append(file)
    return all_receptors


def prepare_receptor(filename, new_filename):
    prepare_receptor = f"/flask_server/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor -r {filename} -o {new_filename}.pdbqt -A hydrogens -v -U nphs_lps_waters"
    print(prepare_receptor)

    ps = subprocess.Popen(
        [prepare_receptor], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    stdout, stderr = ps.communicate()
    # print(stdout)
    if stderr != None:
        print(stderr, "error")

    return


def receptor(filename):
    # print(filename, 'is converted to pdbqt')
    new_filename = "pdbqt_folder/" + filename[:-4]

    cmd.reinitialize()
    cmd.load(filename)
    area = cmd.get_area()

    center = (
        config.center
    )  # overwrite variable to use the same coordinates in all enzymes

    if config.metal_containing:
        print("Metal coordinates identified")
        cmd.select("metals")
        xyz = cmd.get_coords("sele")
        try:
            cen = xyz.tolist()
        except:
            print("no metal found")
            sys.exit()
        center = cen[0]

    if config.NADP_cofactor:
        center = find_box_center(filename, config.center)

    cmd.reinitialize()
    print("starting conversion of ", filename, " to pdbqt")
    prepare_receptor(filename, new_filename)

    return center


def prepare_all_receptors():
    all_receptors = receptor_list(".pdb")
    # enter here check for pdbqt files if len found eqaul to pdb don't do the below stuff and return bc from for loop

    if len(receptor_list(".pdbqt")) < len(all_receptors):
        print("pdb to pdbqt conversion of receptors:", len(all_receptors))

        if config.NADP_cofactor:
            print("Trying to identify atom C5N of cofactor NADP")

        bc = []
        for file in all_receptors:
            center = receptor(file)
            bc.append(center)

        # with multiprocessing.Pool(processes=4) as pool:
        #    bc = pool.starmap(receptor,[[x] for x in all_receptors])

        with open("box_centers", "wb") as fp:  # Pickling
            pickle.dump(bc, fp)

        print("conversion finished")

    else:
        print("receptor preparation skipped. Already done.")
        with open("box_centers", "rb") as fp:  # Unpickling
            bc = pickle.load(fp)

    print(bc)  # possible box coordinates
    return (all_receptors, bc)


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def analyse_affinities_all(receptor, rec_path, rec_sub_path, df_all):
    if len([file for file in os.listdir(rec_path) if file.endswith(".log")]) == 0:
        print(
            "no log files found. Docking failed. Make sure that there are all pdbqt files in the folder"
        )
        return
    error = 0
    for file in os.listdir(rec_path):
        if file.endswith(".log"):
            f = open(rec_sub_path + file)
            lines = f.readlines()
            parts = file.split(sep="_")

            try:
                if config.cofactor == "":
                    en = []
                    en1 = []
                    en2 = []
                    en3 = []
                    en4 = []

                    for line in lines:  # only include first docking pose
                        if line.startswith("Estimated Free Energy of Binding  "):
                            en.append(float(line.split()[-3]))

                        if line.startswith("(1)"):
                            en1.append(float(line.split()[-2]))
                        if line.startswith("(2)"):
                            en2.append(float(line.split()[-2]))
                        if line.startswith("(3)"):
                            en3.append(float(line.split()[-2]))
                        if line.startswith("(4)"):
                            en4.append(float(line.split()[-2]))

                    # sorting energies
                    data = {"en": en, "en1": en1, "en2": en2, "en3": en3, "en4": en4}
                    df = pd.DataFrame(data)
                    # Sort the DataFrame based on the values of en
                    df_sorted = df.sort_values(by="en")
                    # Reset the index of the sorted DataFrame
                    df_sorted.reset_index(drop=True, inplace=True)
                    en = df_sorted["en"].tolist()
                    en1 = df_sorted["en1"].tolist()
                    en2 = df_sorted["en2"].tolist()
                    en3 = df_sorted["en3"].tolist()
                    en4 = df_sorted["en4"].tolist()

                    df_energies = pd.DataFrame(
                        [en + en1 + en2 + en3 + en4],
                        columns=["Energies_" + str(x + 1) for x in range(len(en))]
                        + ["Energies1_" + str(x + 1) for x in range(len(en1))]
                        + ["Energies2_" + str(x + 1) for x in range(len(en2))]
                        + ["Energies3_" + str(x + 1) for x in range(len(en3))]
                        + ["Energies4_" + str(x + 1) for x in range(len(en4))],
                    ).astype(float)

                else:
                    print("needs to be done")

                if False == (df_energies.dtypes == np.float64).all():
                    error += 1
                    print("this log file looks strange", file)

            except:
                error += 1
                print("log error", file)
                en = 0
                continue

            df_energies["res_out"] = os.path.join(
                rec_path, file[:-8] + "_out." + config.output_formate
            )

            df_energies["ligand_names"] = parts[-1]
            # split of log file using underline but ligand defined with underline as well

            df_energies["Enzyme_group"] = parts[0]
            df_energies["Log_files"] = file
            df_energies["res_path"] = rec_sub_path
            df_energies["Enzyme_ligand"] = parts[0] + "_" + parts[-3]

            f.close()

            df_all = pd.concat([df_all, df_energies], sort=False)

    df_all = df_all.reset_index(drop=True)

    if config.plot:
        substances_db = str(config.mol2_files)

        fig = plt.figure(figsize=(20, 15))

        ax1 = fig.add_subplot()
        ax1.set_ylabel("binding affinity / (-kcal/mol)")
        ax1.set_title("Docking of " + receptor[:-6] + " with " + substances_db[:-6])

        heights = (-1) * (df_all.Energies_1)
        bars = df_all.Enzyme_group + "_" + df_all.ligand_names
        y_pos = range(len(bars))
        barlist = plt.bar(y_pos, heights)

        df_all["max"] = False

        ind = 0
        # ind = df['Energies'].argmin()
        barlist[ind].set_color("r")
        df_all["max"][ind] = True

        # Rotation of the bars names
        plt.xticks(y_pos, bars, rotation=90)

        plt.savefig(rec_sub_path + receptor[:-6] + "_results.png")
        plt.tight_layout()
        if config.show_plots == False:
            plt.close()

    df_all.to_excel(rec_sub_path + receptor[:-6] + "_results.xlsx")

    return df_all


def prepare_ligands(substances_db, result_path, df_ligands, receptors):
    if not glob.glob(config.result_path + "/" + "*.mol2"):
        if config.fixed_torsion:
            os.system(
                f"obabel {substances_db} -O {config.result_path+'/'}Ligand.mol2 -m"
            )
            lig_sdf = []
            lig_sdf_new = []
            for file in natsorted(glob.glob(os.path.join(result_path, "*.mol2"))):
                lig_sdf.append(file)
            print(len(lig_sdf))

            for idx, n in enumerate(lig_sdf):
                lig_sdf_new.append(
                    os.path.join(result_path, list(df_ligands.ligand_names)[idx])
                    + ".mol2"
                )
                os.rename(
                    n,
                    os.path.join(result_path, list(df_ligands.ligand_names)[idx])
                    + ".mol2",
                )
            os.chdir(config.result_path)
            for lig in lig_sdf_new:
                prepare_ligand = f"/flask_server/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand  -l {lig} -A hydrogens -v -U nphs_lps_waters -Z"
                # prepare_ligand=f'python3 /flask_server/mk_prepare_ligand_new.py -i {lig} --add_hydrogen --pH 7.4'
                # print(prepare_ligand)

                ps = subprocess.Popen(
                    [prepare_ligand],
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                stdout, stderr = ps.communicate()
                if len(stdout) != 0:
                    print(stdout)
                if stderr != None:
                    print(stderr, "error")
            os.chdir(working_dir)
        else:
            os.system(
                f"obabel {substances_db} -O {config.result_path+'/'}Ligand.sdf -m"
            )
            lig_sdf = []
            lig_sdf_new = []
            for file in natsorted(glob.glob(os.path.join(result_path, "*.sdf"))):
                lig_sdf.append(file)
            print(len(lig_sdf))

            for idx, n in enumerate(lig_sdf):
                lig_sdf_new.append(
                    os.path.join(result_path, list(df_ligands.ligand_names)[idx])
                    + ".sdf"
                )
                os.rename(
                    n,
                    os.path.join(result_path, list(df_ligands.ligand_names)[idx])
                    + ".sdf",
                )

            for lig in lig_sdf_new:
                prepare_ligand = f"python /flask_server/mk_prepare_ligand_new.py -i {lig} --add_hydrogen --pH 7.4"
                print(prepare_ligand)

                ps = subprocess.Popen(
                    [prepare_ligand],
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                stdout, stderr = ps.communicate()
                if len(stdout) != 0:
                    print(stdout)
                if stderr != None:
                    print(stderr, "error")

        all_lig_files = []
        full_paths = []
        for file in os.listdir(result_path):
            if file.endswith(".pdbqt") and (file not in receptors):
                full_paths.append(os.path.join(result_path, file))
                all_lig_files.append(file)
        all_lig_files = natsorted(all_lig_files)

        f = open(config.result_path + "/" + "Ligand.txt", "a+")
        f.writelines(["%s\n" % item for item in all_lig_files])
        f.close()
        print(len(all_lig_files))

        return all_lig_files

    else:
        print("Ligand files do already exist")


def prepare_cofactor():
    second_ligand = config.cofactor[:-5] + ".pdbqt"
    prepare_cofactor = f"python /flask_server/mk_prepare_ligand_new.py  -i {config.cofactor} -o {second_ligand}"  #  cd /flask_server/jupyter/Meeko/ and pip install .

    print(" CHECK AGAIN PLS")

    ps = subprocess.Popen(
        [prepare_cofactor], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    stdout, stderr = ps.communicate()
    # print(stdout)
    if stderr != None:
        print(stderr, "error")

    return second_ligand


# execute vina.exe with variables


def vina_docking(lig, file_name, receptor_v, center_v):
    center_x = center_v[0]
    center_y = center_v[1]
    center_z = center_v[2]

    (size_x, size_y, size_z) = config.size
    if config.cofactor != "":
        second_ligand = config.cofactor[:-5] + ".pdbqt"
    else:
        second_ligand = ""

    if config.gpu_vina:
        file_name = working_dir + "/" + file_name
        lig = working_dir + "/" + lig
        receptor_v = working_dir + "/" + receptor_v
        vina_docking = f"/flask_server/QuickVina2-GPU/QuickVina2-GPU --thread 8000 --receptor {receptor_v} --ligand {lig} {second_ligand} --seed 42 --center_x {center_x} --center_y {center_y} --center_z {center_z} --out {file_name}_out.{config.output_formate} --size_x {size_x} --size_y {size_y} --size_z {size_z} --num_modes {config.num_modes} --energy_range {config.energy_range}"

    else:
        # check if lig is present
        if not os.path.exists(lig):
            print("ligand not found: ", lig)
            return
        vina_docking = f"/flask_server/vina_1.2.3_linux_x86_64 --verbosity 2 --receptor {receptor_v} --ligand {lig} {second_ligand} --cpu 64 --seed 42 --center_x {center_x} --center_y {center_y} --center_z {center_z} --out {file_name}_out.{config.output_formate} --size_x {size_x} --size_y {size_y} --size_z {size_z} --num_modes {config.num_modes} --exhaustiveness {config.exhaustiveness} --energy_range {config.energy_range}"
        print(vina_docking)

    log_file = f"{file_name}_log.log"
    # print(vina_docking)
    ps = subprocess.Popen(
        [vina_docking], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    stdout, stderr = ps.communicate()
    os.chdir(working_dir)
    try:
        decoded_text = stdout.decode("utf-8", errors="replace")
    except UnicodeDecodeError as e:
        print("Error decoding text:", e)
    log(log_file, decoded_text)
    if stderr != None:
        try:
            decoded_text = stderr.decode("utf-8", errors="replace")
        except UnicodeDecodeError as e:
            print("Error decoding text:", e)
        log(log_file, decoded_text)


def virtual_screen(receptor, center):
    sub_path = config.result_path + "/"

    rec_path = os.path.join(sub_path, receptor[:-6])
    rec_sub_path = sub_path + receptor[:-6] + "/"

    if not os.path.exists(rec_path):
        os.makedirs(rec_path)

    with open(sub_path + "Ligand.txt", "r") as f:
        lines = f.readlines()
        print("number of ligands", len(lines))
        f.close()

    fn = []
    lig = []
    for i in lines:
        i = i.rstrip("\n")
        file_name = i.split(".")[0]  # give file name without formate ending
        file_name = receptor[:-6] + "_" + file_name
        file_name = (
            rec_sub_path + file_name
        )  # log file into specific receptor result folder
        i = (
            sub_path + i
        )  # take ligand from sub_path rather than result file for receptor
        lig.append(i)
        fn.append(file_name)

    receptor_v = ["pdbqt_folder/" + receptor] * len(
        lig
    )  # preparing receptor as list for zipping
    center_v = [center] * len(lig)  # preparing center as list for zipping
    if config.gpu_vina:
        os.chdir("/flask_server/QuickVina2-GPU")

    if len(lig) > len(config.receptors):
        print("more ligands than receptors")
        processes = int(50 / config.exhaustiveness)
        if processes == 0:
            processes = 1
        with multiprocessing.Pool(
            processes=processes
        ) as pool:  # limit multiprocessing to 50 cores depending on exhaustiveness
            pool.starmap(
                vina_docking, tqdm(zip(lig, fn, receptor_v, center_v), total=len(lig))
            )  # tqdm and total = len(lig) for progressbar

    else:
        print("more receptors than ligands")

        for i in range(len(lig)):
            vina_docking(lig[i], fn[i], receptor_v[i], center_v[i])

    if config.gpu_vina:
        os.chdir(working_dir)

    df = pd.DataFrame()
    if config.gpu_vina:
        return df
    else:
        df = analyse_affinities_all(receptor, rec_path, rec_sub_path, df)
        return df


def delete_all_pdbqt_files():
    for file in os.listdir():
        if file.endswith(".pdbqt"):
            os.remove(file)
    print("all pdbqt files deleted")


def plot_binding_affinity(df_xlsx, substances_db, destination):
    fig = plt.figure(figsize=(35, 15))

    ax1 = fig.add_subplot()
    ax1.set_ylabel("binding affinity / (-kcal/mol)")
    ax1.set_title("Docking of all receptors with " + substances_db[:-4])

    heights = (-1) * (df_xlsx.Energies_1)
    bars = df_xlsx.Enzyme_ligand
    y_pos = range(len(bars))
    barlist = plt.bar(y_pos, heights)

    print(
        "for "
        + substances_db[:-4]
        + " the best affinity is "
        + str(df_xlsx.Energies_1.min())
        + " kcal/mol"
    )

    # color best affinities in red
    best = df_xlsx[df_xlsx.Energies_1 == df_xlsx.Energies_1.min()]
    best_list = list(best.index)

    for ind in best_list:
        barlist[ind].set_color("r")

    # Rotation of the bars names
    plt.xticks(y_pos, bars, rotation=90)
    plt.tight_layout()

    plt.savefig(os.path.join(destination, "All_results.png"))
    if config.show_plots == False:
        plt.close()
    return best


def summarize_results(res1, substances_db):
    df_xlsx = pd.DataFrame()
    for i in range(len(res1)):
        # res1[i]['Enzyme_group']=i
        df_xlsx = df_xlsx.append(res1[i])

    df_xlsx["out"] = df_xlsx["res_path"] + df_xlsx["Log_files"]
    try:
        df_xlsx["enzyme_reaction"] = df_xlsx.apply(
            lambda row: row.Enzyme_group + "_" + row.ligand_names, axis=1
        )
    except:
        print("no enzyme groups or ligand names")
    df_xlsx = df_xlsx.reset_index()

    # Destination path create summary folder which harbors all the results summarized
    path = os.getcwd()
    destination = os.path.join(config.result_path + "/", "Results_summary")

    # move pdb files which end with aligned.pdb or pymol.pdb to results_summary folder
    if not os.path.exists(destination):
        os.makedirs(destination)

    for file in os.listdir(path):
        if file.endswith("aligned.pdb") or file.endswith("pymol.pdb"):
            shutil.move(file, config.result_path + "/")
        if file.endswith(".txt"):
            shutil.move(file, config.result_path + "/")

    best = plot_binding_affinity(df_xlsx, substances_db, destination)
    df_xlsx.to_excel(os.path.join(destination, "Overview_binding_affinity.xlsx"))

    # copy best files to Results_summary folder

    for i in range(len(best)):
        log_file = os.path.join(path, best.iloc[i]["out"])
        log_file_dest = os.path.join(destination, best.iloc[i]["Log_files"])

        out_file = os.path.join(
            path, best.iloc[i]["out"][:-8] + "_out.pdbqt"
        )  # maybe add into result file instead
        out_file_dest = os.path.join(
            destination, best.iloc[i]["Log_files"][:-8] + "_out.pdbqt"
        )  # maybe add into result file instead
        try:
            shutil.copyfile(log_file, log_file_dest)  # copy from source to destination
            shutil.copyfile(out_file, out_file_dest)  # copy from source to destination
        except:
            print("copy error")

        # convert all result files into mol2

    res_outs = df_xlsx["res_out"]
    res_mols = [[]] * len(res_outs)
    ligand_mols = glob.glob("Results/*.mol2", recursive=True)
    mol_files = glob.glob("Results/**/*.mol2", recursive=True)

    for i in range(len(res_outs)):
        if len(mol_files) < (len(res_outs) + len(ligand_mols)):
            pdbqt_to_mol = f"obabel {res_outs[i]} -O {res_outs[i][:-6]}.mol2"

            ps = subprocess.Popen(
                pdbqt_to_mol,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            stdout, stderr = ps.communicate()
            if stderr != None:
                print(stderr)

        res_mols[i] = res_outs[i][:-10] + ".mol2"

    df_xlsx["res_mols"] = res_mols
    # df_xlsx['Enzyme_group'] = df_xlsx.apply(lambda row: row.Ligand.split(sep='_')[0],axis = 1)
    print("adding no distancese")
    # df_xlsx = add_distances(df_xlsx)

    return df_xlsx


def gpu_vina_summarize():
    """summarize the results of the docking using gpu_vina. Only the best scoring pose is taken into account"""
    log_files = glob.glob(
        f"{config.result_path}/**/*.log", recursive=False
    )  # list of all log files in the folder
    df_log = pd.DataFrame(
        columns=["variant", "state", "affinity", "vina_rmsd_lb", "vina_rmsd_ub"]
    )  # create a dataframe to store the log file data
    for log_file in log_files:
        with open(log_file, "r") as f:
            lines = f.readlines()
            for count, line in enumerate(lines):
                if "-----+------------+----------+----------" in line:
                    while lines[count + 1] != "Writing output ... done.\n":
                        state = lines[count + 1].split()[0]
                        affinity = lines[count + 1].split()[1]
                        vina_rmsd_lb = lines[count + 1].split()[2]
                        vina_rmsd_ub = lines[count + 1].split()[3]
                        df_log = df_log.append(
                            {
                                "variant": log_file.split("/")[-2],
                                "ligand": log_file.split("/")[-1].split("_")[-2],
                                "state": state,
                                "affinity": affinity,
                                "vina_rmsd_lb": vina_rmsd_lb,
                                "vina_rmsd_ub": vina_rmsd_ub,
                            },
                            ignore_index=True,
                        )
                        count += 1
    return df_log


def post_processing():
    path = os.getcwd()

    def receptor_list(formate):
        all_receptors = []
        for file in os.listdir(path):
            if file.endswith(formate):
                all_receptors.append(file)
        return all_receptors

    def post_analysis(receptor):
        sub_path = config.result_path + "/"
        rec_path = os.path.join(sub_path, receptor[:-6])
        rec_sub_path = sub_path + receptor[:-6] + "/"
        df = pd.DataFrame()
        df = analyse_affinities_all(receptor, rec_path, rec_sub_path, df)
        return df

    all_receptors = receptor_list(".pdbqt")

    res1 = [[]] * len(all_receptors)
    for i in range(len(all_receptors)):
        res1[i] = post_analysis(all_receptors[i])
    print("data loaded")
    return res1


def smiles_to_dict(file):
    molecule_library = {}
    f = open(file)
    lines = f.readlines()
    for line in lines:
        try:
            key = line.split(sep="\t")[1].strip(".log\n")
            item = line.split(sep="\t")[0]
            molecule_library[key] = item
        except:
            key = line.split(sep=" ")[1].strip("\n")
            item = line.split(sep=" ")[0]
            molecule_library[key] = item

    df_ligands = pd.DataFrame(columns=["ligand_names", "molecule_smiles"])
    for key, value in molecule_library.items():
        df_ligands = df_ligands.append(
            {"ligand_names": key, "molecule_smiles": value}, ignore_index=True
        )

        # generate molecules as sdf file
    mols = [Chem.MolFromSmiles(x) for x in df_ligands["molecule_smiles"]]
    with Chem.SDWriter("ligands.sdf") as w:
        for idx, m in enumerate(mols):
            w.write(m)

    if not glob.glob("ligands3D.sdf"):
        generate_3d_sdf = f"obabel ligands.sdf -O ligands3D.sdf --gen3D -aa"
        ps = subprocess.Popen(
            generate_3d_sdf,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        stdout, stderr = ps.communicate()
        print(stdout, stderr)
    return df_ligands


def align_pdbs():
    if not os.path.exists("original_pdbs"):
        os.makedirs("original_pdbs")
    all_pdb = [
        pdb
        for pdb in os.listdir()
        if pdb.endswith(".pdb") and not pdb.endswith("aligned.pdb")
    ]
    t = all_pdb[1:]  # get all pdbs without first one
    cmd.reinitialize()
    for pdb in all_pdb:
        cmd.load(pdb)

    for pdb in t:
        cmd.align(pdb[:-4], all_pdb[0][:-4])
    for pdb in all_pdb:
        os.system(f"mv {pdb} original_pdbs/{pdb}")
        cmd.save(pdb[:-4] + "_aligned.pdb", selection=pdb[:-4])
    cmd.reinitialize()

    print(len(all_pdb), " pdbs aligned and original moved to orignal_pdb folder")


def pymol_selection_pdbs():
    """select all pdbs in pymol and save them as pdb files. This is done as a sanitization step and to be coherent with AlphaDock"""
    if not os.path.exists("original_pdbs"):
        os.makedirs("original_pdbs")
    all_pdb = [
        pdb
        for pdb in os.listdir()
        if pdb.endswith(".pdb")
        and not pdb.endswith("aligned.pdb")
        and not pdb.endswith("pymol.pdb")
    ]
    cmd.reinitialize()
    for pdb in all_pdb:
        cmd.load(pdb)
        os.system(f"mv {pdb} original_pdbs/{pdb}")
        cmd.save(pdb[:-4] + "_pymol.pdb", selection=pdb[:-4])
    cmd.reinitialize()
    print(
        len(all_pdb), " pdbs selected by pymol and original moved to orignal_pdb folder"
    )


# main function for docking start all
start = time.time()
substances_db = ""

os.chdir(
    config.working_dir
)  # move all homology models to the current working directory containing all receptors as pdb files including ligand library as sdf
(center_x, center_y, center_z) = config.center

pdbqt_folder = os.path.join(config.working_dir, "pdbqt_folder")
if not os.path.exists(pdbqt_folder):
    os.makedirs(pdbqt_folder)

pymol_selection_pdbs()
print("pymol selection done")
if config.align:
    # if multiple pdb files exist and they are homologs to each other it might make sense to align them before docking to set the box coordinates at the same postion in all cases
    # if true this will be done using pymol
    align_pdbs()
else:
    print("no alignment performed")

if ligand_files.endswith(".smiles"):
    df_ligands = smiles_to_dict(ligand_files)
    substances_db = "ligands3D.sdf"
if ligand_files.endswith(".sdf"):
    print(
        "sdf file conversion function experimental. Make sure file is named ligands3D.sdf contains all ligands and the ligand names are given"
    )
    substances_db = "ligands3D.sdf"
    print(substances_db)
    df_ligands = sdf_to_df(substances_db)

if ligand_files == "":
    print("no ligand library provided as smiles or sdf file")

if mol2_files != []:
    print("mol2 files provided, converting to sdf")
    df_ligands = mol2_to_df(mol2_files)
    substances_db = "ligands3D.sdf"


if substances_db == "":
    print("no mol2/sdf/smiles files as ligand provided")
    sys.exit()

if slopes_file != "":
    print("slopes file provided, code missing")

prepare_docking_files()


result_path = os.path.join(config.working_dir, config.result_path)


all_receptors, bc = prepare_all_receptors()


receptors = []
for i in all_receptors:
    re = i[:-4] + ".pdbqt"
    receptors.append(re)
config.receptors = receptors
print("receptors", receptors)
print("length receptors", len(receptors))
print("ligand database", substances_db)

res1 = len(receptors) * [[]]

print(result_path)
if not os.path.exists(result_path):
    os.makedirs(result_path)

    if config.cofactor != "":
        if config.cofactor.endswith(".mol2"):
            second_ligand = prepare_cofactor()
        else:
            print("cofactor formate wrong. needs to be mol2")
    else:
        second_ligand = ""

    # new way of preparing ligands using meeko approach
    ligands = prepare_ligands(substances_db, result_path, df_ligands, receptors)

    if len(ligands) > len(receptors):
        print("more ligands than receptors")
        for i in tqdm(range(len(receptors)), total=len(receptors)):
            res1[i] = virtual_screen(receptors[i], bc[i])

    else:
        print("more receptors than ligands")

        if config.gpu_vina:
            for i in tqdm(range(len(receptors)), total=len(receptors)):
                res1[i] = virtual_screen(receptors[i], bc[i])

        else:
            with multiprocessing.Pool(
                processes=4
            ) as pool:  # limit multiprocessing to 50 cores
                for result in pool.starmap(
                    virtual_screen, tqdm(zip(receptors, bc), total=len(receptors))
                ):
                    res1.append(result)

    print("Multidocking finished !!!!")

    if config.gpu_vina:
        df = gpu_vina_summarize()
        path = os.getcwd()
        destination = os.path.join(config.result_path + "/", "Results_summary")
        if not os.path.exists(destination):
            os.makedirs(destination)
        df.to_excel(os.path.join(destination, "Overview_binding_affinity.xlsx"))
        print("GPU Vina used data differently summarized. Check Data")

    else:
        df_xlsx = summarize_results(res1, substances_db)

        with open("df_xlsx.pkl", "wb") as fp:  # Pickling
            pickle.dump(df_xlsx, fp)
    done_text = f"It took {round(time.time()-start,2)} seconds for docking of {len(receptors)} receptors and {len(ligands)} ligands."
    print(done_text)
    log(os.path.join(config.working_dir, "Experiment_log_file.txt"), done_text)

    # zip Results folder
    shutil.make_archive("Results", "zip", "Results")

else:
    if os.path.exists("df_xlsx.pkl"):
        with open("df_xlsx.pkl", "rb") as fp:  # Unpickling
            df_xlsx = pickle.load(fp)
    else:
        res1 = post_processing()
        df_xlsx = summarize_results(res1, substances_db)
        with open("df_xlsx.pkl", "wb") as fp:  # Pickling
            pickle.dump(df_xlsx, fp)

    print("It took", round(time.time() - start, 1), "seconds. To re-load results.")

    print("\n Docking already done. Do you want to redo it? Delete Results folder.")
    # zip Results folder
    shutil.make_archive("Results", "zip", "Results")
