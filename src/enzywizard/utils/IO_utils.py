from __future__ import annotations

from Bio.PDB import MMCIFParser, PDBParser, MMCIFIO, PDBIO
from Bio.PDB.Structure import Structure
from pathlib import Path
from Bio.PDB.DSSP import DSSP
from ..utils.logging_utils import Logger
import json
import tempfile
from ..utils.common_utils import convert_to_json_serializable, InlineJSONEncoder, wrap_leaf_lists_as_rawjson, get_clean_filename
from openmm.app import PDBFile,PDBxFile, Modeller
from ..utils.structure_utils import get_single_chain,get_residues_by_chain,get_sequence
from ..utils.conservation_utils import load_msa_sto,load_msa_aligned_fasta,load_msa_a3m, write_sto,write_aligned_fasta,write_a3m
from typing import List, Dict,Any, Optional, Tuple
import subprocess
from rdkit import Chem
from ..utils.substrate_utils import is_valid_mol_3d
from Bio.PDB import StructureBuilder
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
import copy
import numpy as np

def file_exists(path: str | Path) -> bool:
    p = Path(path)
    return p.exists() and p.is_file()

def get_stem(input_path: str | Path) -> str:
    return Path(input_path).stem

MAXFILENAME=100

def check_filename_length(name: str, logger: Logger) -> bool:
    if len(name) > MAXFILENAME:
        logger.print(f"[ERROR] Filename too long (>{MAXFILENAME}): {name}")
        return False
    return True

def load_protein_structure(path: str | Path, protein_name:str, logger: Logger) -> Structure | None:
    p = Path(path)

    try:
        if p.suffix.lower() in {".cif", ".mmcif"}:
            parser = MMCIFParser(QUIET=True)
        elif p.suffix.lower() == ".pdb":
            parser = PDBParser(QUIET=True)
        else:
            logger.print(f"[ERROR] Unsupported format: {p}")
            return None

        return parser.get_structure(protein_name, str(p))

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading structure for {str(p)}: {e}")
        return None

def load_openmm_structure(path: str | Path, logger: Logger) -> PDBFile | PDBxFile | None:
    p = Path(path)

    try:
        if p.suffix.lower() in {".cif", ".mmcif"}:
            return PDBxFile(str(p))
        elif p.suffix.lower() == ".pdb":
            return PDBFile(str(p))
        else:
            logger.print(f"[ERROR] Unsupported format: {p}")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading OpenMM structure for {str(p)}: {e}")
        return None

def structure_to_pdbfile(struct: Structure, logger: Logger, protein_name: str = "structure") -> PDBFile | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            io = PDBIO()
            io.set_structure(struct)
            io.save(str(pdb_path))

            pdb = PDBFile(str(pdb_path))
            return pdb

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Structure to PDBFile: {e}")
        return None

def pdbfile_to_structure(pdb: PDBFile, logger: Logger, protein_name: str = "structure") -> Structure | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            with open(pdb_path, "w") as f:
                PDBFile.writeFile(pdb.topology, pdb.positions, f)

            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(protein_name, str(pdb_path))

            return struct

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert PDBFile to Structure: {e}")
        return None

def modeller_to_structure(modeller: Modeller, logger: Logger, protein_name: str = "structure") -> Structure | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            with open(pdb_path, "w") as f:
                PDBFile.writeFile(modeller.topology, modeller.positions, f)

            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(protein_name, str(pdb_path))

            return struct

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Modeller to Structure: {e}")
        return None

def load_dssp(struct: Structure, logger: Logger) -> DSSP | None:
    try:
        model = next(struct.get_models())

        with tempfile.TemporaryDirectory() as td:
            cif_path = Path(td) / "tmp.cif"

            io = MMCIFIO()
            io.set_structure(model)
            io.save(str(cif_path))

            dssp = DSSP(model,str(cif_path),dssp="mkdssp")
            return dssp

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading dssp: {e}")
        return None

def write_cif(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = MMCIFIO()
    io.set_structure(struct)
    io.save(str(output_path))

def write_pdb(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(struct)
    io.save(str(output_path))

def write_json_from_dict(dict_data: dict, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_data=convert_to_json_serializable(dict_data)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(dict_data, f, indent=2, ensure_ascii=False)

def write_json_from_dict_inline_leaf_lists(dict_data: dict, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_data = convert_to_json_serializable(dict_data)
    dict_data = wrap_leaf_lists_as_rawjson(dict_data)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(
            dict_data,
            f,
            cls=InlineJSONEncoder,
            indent=2,
            ensure_ascii=False
        )

def write_fasta(struct: Structure, output_path: str | Path, logger: Logger) -> bool:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chain = get_single_chain(struct, logger)
    if chain is None:
        return False

    residues = get_residues_by_chain(chain, logger)
    if residues is None:
        return False

    seq = get_sequence(residues, logger)
    if seq is None:
        return False

    header = str(struct.id).strip() if getattr(struct, "id", None) else None
    if header is None:
        logger.print("[ERROR] Structure header is None")
        return False

    try:
        with output_path.open("w", encoding="utf-8") as f:
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to write FASTA to {output_path}: {e}")
        return False

def load_msa(path: str | Path, logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return load_msa_sto(p, logger)

        elif suffix in {".fa", ".fasta", ".afa"}:
            return load_msa_aligned_fasta(p, logger)

        elif suffix == ".a3m":
            return load_msa_a3m(p, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {str(p)}")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Exception in load_msa from {str(p)}: {e}")
        return None

def load_fasta(path: str | Path, logger: Logger) -> Dict[str, str] | None:
    p = Path(path)

    try:
        header: str | None = None
        seq_parts: list[str] = []

        with p.open("r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.rstrip("\n").strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        logger.print(f"[ERROR] Multiple sequences found in FASTA file: {str(p)}")
                        return None

                    header = line[1:].strip()
                else:
                    if header is None:
                        logger.print(f"[ERROR] Invalid FASTA format in {str(p)}: sequence line appears before header.")
                        return None
                    seq_parts.append(line)

        if header is None:
            logger.print(f"[ERROR] No header found in FASTA file: {str(p)}")
            return None

        sequence = "".join(seq_parts)

        if sequence.strip() == "":
            logger.print(f"[ERROR] Empty sequence found in FASTA file: {str(p)}")
            return None

        return {"header": header, "sequence": sequence}

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading FASTA from {str(p)}: {e}")
        return None

def write_msa(msa_list: List[Dict[str, str]], output_path: str | Path, logger: Logger) -> bool:
    p = Path(output_path)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return write_sto(msa_list, p, logger)

        elif suffix in {".fa", ".fasta", ".afa"}:
            return write_aligned_fasta(msa_list, p, logger)

        elif suffix == ".a3m":
            return write_a3m(msa_list, p, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {str(p)}")
            return False

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_msa: {e}")
        return False

def write_hmm(sto_path: str | Path, output_path: str | Path, logger: Logger) -> bool:
    sto_file = Path(sto_path)
    hmm_file = Path(output_path)

    try:
        if not sto_file.exists():
            logger.print(f"[ERROR] Input Stockholm file not found: {str(sto_file)}")
            return False

        hmm_file.parent.mkdir(parents=True, exist_ok=True)

        p = subprocess.run(
            [
                "hmmbuild",
                "--hand",
                str(hmm_file),
                str(sto_file),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if p.returncode != 0:
            logger.print(f"[ERROR] hmmbuild failed for {str(sto_file)}: {p.stderr.strip()}")
            return False

        if (not hmm_file.exists()) or hmm_file.stat().st_size == 0:
            logger.print(f"[ERROR] Failed to generate HMM file: {str(hmm_file)}")
            return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_hmm: {e}")
        return False


def write_sdf(mol_3d: Chem.Mol, sdf_path: str | Path, logger: Logger,) -> bool:
    if not is_valid_mol_3d(mol_3d, logger):
        return False

    try:
        sdf_path = Path(sdf_path)
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        writer = Chem.SDWriter(str(sdf_path))
        conf_id = mol_3d.GetConformer().GetId()
        writer.write(mol_3d, confId=conf_id)
        writer.close()

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to save SDF file.")
            return False

        return True
    except Exception:
        logger.print("[ERROR] Failed to save Mol(3D) to SDF file.")
        return False


def save_substrate_structures(substrate_feature_list: List[Dict[str, Any]],output_dir: str | Path,logger: Logger) -> bool:

    if not isinstance(substrate_feature_list, list):
        logger.print("[ERROR] substrate_feature_list must be a list.")
        return False

    try:
        output_dir = Path(output_dir)

        tasks: List[Tuple[Chem.Mol, Path]] = []

        for item in substrate_feature_list:
            structures = item.get("structures", [])

            for s in structures:
                mol = s.get("structure_mol")
                name = s.get("structure_name")

                if not mol or not name:
                    logger.print("[ERROR] Invalid structure entry.")
                    return False

                clean_name = get_clean_filename(name)
                path = output_dir / f"{clean_name}.sdf"

                tasks.append((mol, path))

        for mol, path in tasks:
            if not write_sdf(mol, path, logger):
                return False

        return True

    except Exception:
        logger.print("[ERROR] Failed to save substrate structures.")
        return False

def write_protein_pdbqt(struct: Structure,pdbqt_path: str | Path,logger: Logger) -> bool:
    try:
        pdbqt_path = Path(pdbqt_path)
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_pdb = Path(tmp_dir) / "protein.pdb"

            # 写 PDB
            write_pdb(struct, tmp_pdb)

            if not tmp_pdb.exists() or tmp_pdb.stat().st_size <= 0:
                logger.print("[ERROR] Failed to write temporary PDB file.")
                return False

            # 调用 Meeko CLI
            p = subprocess.run(
                [
                    "mk_prepare_receptor.py",
                    "--read_pdb", str(tmp_pdb),
                    "--write_pdbqt", str(pdbqt_path),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

            if p.returncode != 0:
                logger.print("[ERROR] mk_prepare_receptor.py failed.")
                return False

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to generate PDBQT file.")
            return False

        return True

    except Exception:
        logger.print("[ERROR] Failed to convert Structure to PDBQT.")
        return False

def write_substrate_pdbqt_from_sdf(sdf_path: str | Path,pdbqt_path: str | Path,logger: Logger) -> bool:
    try:
        sdf_path = Path(sdf_path)
        pdbqt_path = Path(pdbqt_path)
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input SDF file.")
            return False

        # 调用 Meeko CLI
        p = subprocess.run(
            [
                "mk_prepare_ligand.py",
                "-i", str(sdf_path),
                "-o", str(pdbqt_path),
                "--add_index_map",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if p.returncode != 0:
            logger.print("[ERROR] mk_prepare_ligand.py failed.")
            return False

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to generate PDBQT file.")
            return False

        return True

    except Exception:
        logger.print("[ERROR] Failed to convert SDF to PDBQT.")
        return False

def load_sdf_mol_3d(sdf_path: str | Path, logger: Logger) -> Chem.Mol | None:
    try:
        sdf_path = Path(sdf_path)

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input SDF file.")
            return None

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        if supplier is None or len(supplier) == 0:
            logger.print("[ERROR] Failed to load SDF file.")
            return None

        mol = supplier[0]
        if mol is None:
            logger.print("[ERROR] Failed to parse Mol from SDF file.")
            return None

        if mol.GetNumConformers() <= 0:
            logger.print("[ERROR] Input SDF does not contain 3D coordinates.")
            return None

        return mol

    except Exception:
        logger.print("[ERROR] Failed to read Mol(3D) from SDF file.")
        return None

def write_docked_sdf_from_atom_info(
    original_mol_3d: Chem.Mol,
    docked_atom_info_list: List[Dict[str, Any]],
    sdf_path: str | Path,
    logger: Logger,
) -> bool:
    if original_mol_3d is None or original_mol_3d.GetNumConformers() <= 0:
        logger.print("[ERROR] Invalid original Mol(3D).")
        return False

    if not isinstance(docked_atom_info_list, list) or len(docked_atom_info_list) == 0:
        logger.print("[ERROR] Invalid docked_atom_info_list.")
        return False

    try:
        sdf_path = Path(sdf_path)
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        atom_num = original_mol_3d.GetNumAtoms()

        used_original_atom_index_set = set()
        kept_original_atom_index_list: List[int] = []

        for item in docked_atom_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))

            if original_atom_index <= 0 or original_atom_index > atom_num:
                logger.print("[ERROR] Invalid original atom index in docked_atom_info_list.")
                return False

            if original_atom_index in used_original_atom_index_set:
                logger.print("[ERROR] Duplicate original atom index in docked_atom_info_list.")
                return False

            used_original_atom_index_set.add(original_atom_index)
            kept_original_atom_index_list.append(original_atom_index)

        kept_original_atom_index_list.sort()

        old_to_new_index_dict: Dict[int, int] = {}
        for new_index, old_index in enumerate(kept_original_atom_index_list):
            old_to_new_index_dict[old_index] = new_index

        rw_mol = Chem.RWMol()
        new_conf = Chem.Conformer(len(kept_original_atom_index_list))

        for old_index in kept_original_atom_index_list:
            old_atom = original_mol_3d.GetAtomWithIdx(old_index - 1)
            new_atom = Chem.Atom(old_atom.GetAtomicNum())
            new_atom.SetFormalCharge(old_atom.GetFormalCharge())
            new_atom.SetIsAromatic(old_atom.GetIsAromatic())
            new_atom.SetChiralTag(old_atom.GetChiralTag())
            new_atom.SetNumExplicitHs(old_atom.GetNumExplicitHs())
            new_atom.SetNoImplicit(old_atom.GetNoImplicit())
            new_atom.SetNumRadicalElectrons(old_atom.GetNumRadicalElectrons())

            rw_mol.AddAtom(new_atom)

        kept_old_index_set = set(kept_original_atom_index_list)

        for bond in original_mol_3d.GetBonds():
            begin_old_index = int(bond.GetBeginAtomIdx()) + 1
            end_old_index = int(bond.GetEndAtomIdx()) + 1

            if begin_old_index in kept_old_index_set and end_old_index in kept_old_index_set:
                begin_new_index = old_to_new_index_dict[begin_old_index]
                end_new_index = old_to_new_index_dict[end_old_index]
                rw_mol.AddBond(begin_new_index, end_new_index, bond.GetBondType())

        for item in docked_atom_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))
            x = float(item.get("x", 0.0))
            y = float(item.get("y", 0.0))
            z = float(item.get("z", 0.0))

            new_atom_index = old_to_new_index_dict[original_atom_index]
            new_conf.SetAtomPosition(new_atom_index, (x, y, z))

            original_atom_name = str(item.get("original_atom_name", "")).upper()
            new_atom = rw_mol.GetAtomWithIdx(new_atom_index)

            if original_atom_name and new_atom.GetSymbol().upper() != original_atom_name:
                logger.print("[ERROR] Atom name mismatch when writing docked SDF.")
                return False

        mol = rw_mol.GetMol()
        mol.RemoveAllConformers()
        mol.AddConformer(new_conf, assignId=True)

        writer = Chem.SDWriter(str(sdf_path))
        conf_id = new_conf.GetId()
        writer.write(mol, confId=conf_id)
        writer.close()

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to save docked SDF file.")
            return False

        return True

    except Exception as e:
        import traceback
        logger.print(f"[ERROR] Failed to write docked atom information to SDF file: {e}")
        logger.print(traceback.format_exc())
        return False

def write_docked_complex_cif(
    docking_result: Dict[str, Any],
    struct: Structure,
    protein_name: str,
    output_dir: str | Path,
    logger: Logger,
) -> str | None:
    if not isinstance(docking_result, dict):
        logger.print("[ERROR] Invalid docking_result.")
        return None

    if struct is None:
        logger.print("[ERROR] struct is None.")
        return None

    if not isinstance(protein_name, str) or not protein_name.strip():
        logger.print("[ERROR] Invalid protein_name.")
        return None

    if not isinstance(output_dir, (str, Path)):
        logger.print("[ERROR] output_dir must be a str or Path.")
        return None

    try:
        substrate_names = docking_result.get("substrate_names", "")
        docked_list = docking_result.get("docked_substrate_info_list", [])

        if not isinstance(substrate_names, str) or not substrate_names.strip():
            logger.print("[ERROR] Invalid docking_result['substrate_names'].")
            return None

        if not isinstance(docked_list, list) or len(docked_list) == 0:
            logger.print("[ERROR] Invalid docking_result['docked_substrate_info_list'].")
            return None

        protein_chain = get_single_chain(struct, logger)
        if protein_chain is None:
            return None

        residue_info_list = get_residues_by_chain(protein_chain, logger)
        if residue_info_list is None or len(residue_info_list) == 0:
            logger.print("[ERROR] Failed to get valid protein residues.")
            return None

        max_protein_resseq = max(res_id[1] for res_id, _, _ in residue_info_list)
        ligand_resseq_start = max_protein_resseq + 1

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        complex_name = f"docked_{protein_name}_{substrate_names}"
        cif_path = output_dir / f"{complex_name}.cif"

        builder = StructureBuilder.StructureBuilder()
        builder.init_structure("complex")
        builder.init_model(0)
        builder.init_chain("A")

        new_struct = builder.get_structure()
        new_model: Model = new_struct[0]
        new_protein_chain: Chain = new_model["A"]

        for residue in protein_chain.get_residues():
            new_protein_chain.add(copy.deepcopy(residue))

        ligand_chain = Chain("L")
        new_model.add(ligand_chain)

        for i, ligand_info in enumerate(docked_list):
            atom_info_list = ligand_info.get("atom_info_list", [])

            if not isinstance(atom_info_list, list) or len(atom_info_list) == 0:
                logger.print(f"[ERROR] Invalid atom_info_list for ligand index {i}.")
                return None

            ligand_res_id = ("H", ligand_resseq_start + i, " ")
            ligand_residue = Residue(ligand_res_id, "LIG", " ")

            for atom_j, atom_item in enumerate(atom_info_list, start=1):
                if not isinstance(atom_item, dict):
                    logger.print(f"[ERROR] Invalid atom item for ligand index {i}.")
                    return None

                original_atom_name = str(atom_item.get("original_atom_name", "")).strip().upper()
                if not original_atom_name:
                    logger.print(f"[ERROR] Missing original_atom_name for ligand index {i}.")
                    return None

                try:
                    x = float(atom_item["x"])
                    y = float(atom_item["y"])
                    z = float(atom_item["z"])
                except Exception:
                    logger.print(
                        f"[ERROR] Invalid atom coordinates for ligand index {i}, atom {original_atom_name}."
                    )
                    return None

                original_atom_index = atom_item.get("original_atom_index", atom_j)
                try:
                    serial_number = int(original_atom_index)
                except Exception:
                    serial_number = atom_j

                # residue 内 atom name 必须唯一
                unique_atom_name = f"A{atom_j}"
                if len(unique_atom_name) > 4:
                    unique_atom_name = unique_atom_name[:4]

                # element 尽量从 original_atom_name 推断
                element = original_atom_name[:2].strip().capitalize()
                if len(element) == 0:
                    element = "X"
                elif len(element) == 2 and element[1].isdigit():
                    element = element[0]
                elif len(element) == 2 and not element[1].islower():
                    element = element[0]

                atom = Atom(
                    name=unique_atom_name,
                    coord=np.array([x, y, z], dtype=float),
                    bfactor=1.0,
                    occupancy=1.0,
                    altloc=" ",
                    fullname=unique_atom_name.rjust(4),
                    serial_number=serial_number,
                    element=element,
                )
                ligand_residue.add(atom)

            ligand_chain.add(ligand_residue)

        write_cif(new_struct, cif_path)

        if not cif_path.exists() or cif_path.stat().st_size == 0:
            logger.print("[ERROR] Failed to write complex CIF.")
            return None

        return str(cif_path)

    except Exception:
        logger.print(f"[ERROR] Failed to build docked complex CIF:")
        return None