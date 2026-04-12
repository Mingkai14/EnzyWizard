from typing import TypedDict, List, Dict, Literal, Required, NotRequired

CLEAN_TYPE = "enzywizard_clean"
MUTCLEAN_TYPE = "enzywizard_mutclean"
AAPROPS_TYPE = "enzywizard_aaprops"
HYDROCLUSTER_TYPE = "enzywizard_hydrocluster"
ENERGY_TYPE = "enzywizard_energy"
FLEXIBILITY_TYPE = "enzywizard_flexibility"
DISORDER_TYPE = "enzywizard_disorder"
CONSERVATION_TYPE = "enzywizard_conservation"
EMBEDDING_TYPE = "enzywizard_embedding"
POCKET_TYPE = "enzywizard_pocket"
SUBSTRATE_TYPE = "enzywizard_substrate"
DOCK_TYPE = "enzywizard_dock"
INTERACTION_TYPE = "enzywizard_interaction"
INTEGRATE_TYPE = "enzywizard_integrate"


'''
clean_report
'''


class ResidueInfo(TypedDict):
    aa_id: int
    aa_name: str
    hydrogen_atom_count: int


class ResidueMapping(TypedDict):
    old_residue: ResidueInfo
    new_residue: ResidueInfo


class CleanStatistics(TypedDict):
    changed_resname: int
    removed_nonstd: int
    removed_missing_bb: int
    removed_missing_heavy_atoms: int
    removed_bad_occ: int
    removed_inscodes: int
    kept_residues: int


class EnzyWizardCleanOutput(TypedDict):
    output_type: Literal["enzywizard_clean"]
    amino_acid_mapping_old_to_new: List[ResidueMapping]
    clean_statistics: CleanStatistics

'''
'''

'''
mutclean_report
'''

class EnzyWizardMutCleanOutput(TypedDict):
    output_type: Literal["enzywizard_mutclean"]

    amino_acid_substitution: str
    cleaned_amino_acid_substitution: str

    wt_amino_acid_mapping_old_to_new: List[ResidueMapping]
    wt_clean_statistics: CleanStatistics

    mut_amino_acid_mapping_old_to_new: List[ResidueMapping]
    mut_clean_statistics: CleanStatistics

'''
'''

'''
aaprops_report
'''

class AAPropsEntry(TypedDict):
    aa_id: int
    aa_name: str

    aa_name_one_hot: List[int]

    aa_class: str
    aa_class_one_hot: List[int]

    aa_ss: str
    aa_ss_one_hot: List[int]

    aa_rsa: float
    aa_phi: float
    aa_psi: float

    aa_net_charge: float
    aa_pka: float
    aa_volume: float
    aa_hydrophobicity: float
    aa_molecular_weight: float
    aa_pi: float

    aa_coord: List[float]  # [x, y, z]


class AAPropsStatistics(TypedDict):
    aa_name_statistics: Dict[str, int]
    aa_class_statistics: Dict[str, int]
    aa_ss_statistics: Dict[str, int]


class EnzyWizardAAPropsOutput(TypedDict):
    output_type: Literal["enzywizard_aaprops"]
    aa_props: List[AAPropsEntry]
    aa_props_statistics: AAPropsStatistics

'''
'''

'''
hydrocluster_report
'''
class ResidueIDName(TypedDict):
    aa_id: int
    aa_name: str


class HydrophobicClusterEntry(TypedDict):
    area: float
    residues: List[ResidueIDName]


class EnzyWizardHydroClusterOutput(TypedDict):
    output_type: Literal["enzywizard_hydrocluster"]
    hydrophobic_cluster: List[HydrophobicClusterEntry]

'''
'''

'''
energy_report
'''

class EnergyTerms(TypedDict):
    total_potential_energy: float
    harmonic_bond_force: float
    harmonic_angle_force: float
    custom_bond_force: float
    custom_torsion_force: float
    custom_nonbonded_force: float
    nonbonded_force: float
    periodic_torsion_force: float
    cmap_torsion_force: float


class EnzyWizardEnergyOutput(TypedDict):
    output_type: Literal["enzywizard_energy"]
    energy_terms: EnergyTerms

'''
'''

'''
flexibility_report
'''

class ResidueRMSFEntry(TypedDict):
    aa_id: int
    aa_name: str
    rmsf: float


class EnzyWizardFlexibilityOutput(TypedDict):
    output_type: Literal["enzywizard_flexibility"]
    protein_rmsf: List[ResidueRMSFEntry]

'''
'''

'''
disorder_report
'''



class DisorderRegionEntry(TypedDict):
    length: int
    residues: List[ResidueIDName]


class EnzyWizardDisorderOutput(TypedDict):
    output_type: Literal["enzywizard_disorder"]
    disorder_regions: List[DisorderRegionEntry]

'''
'''

'''
conservation_report
'''

class ConservationEntry(TypedDict):
    aa_id: int
    aa_name: str

    hmm_emission_log_score: float
    emission_probability: float
    conservation_score: float


class EnzyWizardConservationOutput(TypedDict):
    output_type: Literal["enzywizard_conservation"]
    conservation_scores: List[ConservationEntry]

'''
'''

'''
embedding_report
'''


class EmbeddingEntry(TypedDict):
    aa_id: int
    aa_name: str
    embedding: List[float]


class EnzyWizardEmbeddingOutput(TypedDict):
    output_type: Literal["enzywizard_embedding"]
    embeddings: List[EmbeddingEntry]

'''
'''

'''
pocket_report
'''


class PocketRegionEntry(TypedDict):
    volume: float
    n_spheres: int
    residues: List[ResidueIDName]
    pocket_center_coord: List[float]
    pocket_box_boundaries: List[float]


class EnzyWizardPocketOutput(TypedDict):
    output_type: Literal["enzywizard_pocket"]
    pocket_regions: List[PocketRegionEntry]

'''
'''


'''
substrate_report
'''

class SubstrateStructureEntry(TypedDict):
    structure_name: str
    structure_energy: float


class SubstrateEntry(TypedDict):
    substrate_name: str
    smiles: str
    fingerprint: List[Literal[0, 1]]

    num_atoms: int
    mol_weight: float
    logp: float

    structures: List[SubstrateStructureEntry]


class EnzyWizardSubstrateOutput(TypedDict):
    output_type: Literal["enzywizard_substrate"]
    substrates: List[SubstrateEntry]

'''
'''


'''
dock_report
'''


class DockedSubstrateEntry(TypedDict):
    substrate_name: str
    conformation_name: str
    docked_center_coord: List[float]


class DockedResult(TypedDict):
    complex_name: str
    docking_score: float
    substrate_names: str
    docking_box_center: List[float]
    docking_box_size: List[float]
    docked_substrates: List[DockedSubstrateEntry]


class EnzyWizardDockOutput(TypedDict):
    output_type: Literal["enzywizard_dock"]
    docked_result: DockedResult

'''
'''

'''
interaction_report
'''


class AminoAcidNode(TypedDict):
    aa_index: int
    aa_name: str
    node_type: Literal["amino_acid"]


class SubstrateNode(TypedDict):
    substrate_index: int
    substrate_name: str
    node_type: Literal["substrate"]


InteractionNode = AminoAcidNode | SubstrateNode


class InteractionEntry(TypedDict):
    interaction: Literal["HBOND", "IONIC", "VDW", "PIPISTACK", "PICATION", "SSBOND"]
    node1: InteractionNode
    node2: InteractionNode


class InteractionTypeCount(TypedDict):
    HBOND: int
    IONIC: int
    VDW: int
    PIPISTACK: int
    PICATION: int
    SSBOND: int


class InteractionStatisticsBlock(TypedDict):
    count: InteractionTypeCount
    unique_pair_count: InteractionTypeCount


class InteractionStatistics(TypedDict):
    overall: InteractionStatisticsBlock
    intra_protein: InteractionStatisticsBlock
    protein_substrate: InteractionStatisticsBlock


class EnzyWizardInteractionOutput(TypedDict):
    output_type: Literal["enzywizard_interaction"]
    interactions: List[InteractionEntry]
    interactions_statistics: InteractionStatistics

'''
'''

'''
integrate_report
'''


class IntegratedOverallStatistics(TypedDict):
    aa_name_count: NotRequired[List[int]]
    aa_class_count: NotRequired[List[int]]
    aa_ss_count: NotRequired[List[int]]

    total_potential_energy: NotRequired[float]
    harmonic_bond_force: NotRequired[float]
    harmonic_angle_force: NotRequired[float]
    custom_bond_force: NotRequired[float]
    custom_torsion_force: NotRequired[float]
    custom_nonbonded_force: NotRequired[float]
    nonbonded_force: NotRequired[float]
    periodic_torsion_force: NotRequired[float]
    cmap_torsion_force: NotRequired[float]

    docking_score: NotRequired[float]

    hbond_count: NotRequired[int]
    ionic_count: NotRequired[int]
    vdw_count: NotRequired[int]
    pipistack_count: NotRequired[int]
    pication_count: NotRequired[int]
    ssbond_count: NotRequired[int]


class IntegratedEdge(TypedDict):
    interaction: Literal["HBOND", "IONIC", "VDW", "PIPISTACK", "PICATION", "SSBOND"]
    interaction_one_hot: List[int]
    interaction_count: int


class IntegratedAminoAcidNode(TypedDict):
    node_id: Required[int]
    node_type: Required[Literal["amino_acid"]]
    node_type_one_hot: Required[List[int]]

    aa_index: Required[int]
    aa_name: Required[str]
    aa_coord: NotRequired[List[float]]

    aa_class: NotRequired[str]
    aa_ss: NotRequired[str]
    aa_rsa: NotRequired[float]
    aa_phi: NotRequired[float]
    aa_psi: NotRequired[float]
    aa_net_charge: NotRequired[float]
    aa_pka: NotRequired[float]
    aa_volume: NotRequired[float]
    aa_hydrophobicity: NotRequired[float]
    aa_molecular_weight: NotRequired[float]
    aa_pi: NotRequired[float]

    rmsf: NotRequired[float]
    conservation_score: NotRequired[float]

    aa_name_one_hot: NotRequired[List[int]]
    aa_class_one_hot: NotRequired[List[int]]
    aa_ss_one_hot: NotRequired[List[int]]
    embedding: NotRequired[List[float]]

    is_in_hydrophobic_cluster: NotRequired[bool]
    is_in_disorder_region: NotRequired[bool]
    is_in_pocket: NotRequired[bool]


class IntegratedSubstrateNode(TypedDict):
    node_id: Required[int]
    node_type: Required[Literal["substrate"]]
    node_type_one_hot: Required[List[int]]

    substrate_index: Required[int]
    substrate_name: Required[str]

    smiles: NotRequired[str]
    num_atoms: NotRequired[int]
    mol_weight: NotRequired[float]
    logp: NotRequired[float]
    docked_center_coord: NotRequired[List[float]]
    fingerprint: NotRequired[List[Literal[0, 1]]]


IntegratedNode = IntegratedAminoAcidNode | IntegratedSubstrateNode


class IntegratedSingleNodeEntry(TypedDict):
    node_1: Required[IntegratedNode]


class IntegratedEdgeNodeEntry(TypedDict):
    edge: Required[IntegratedEdge]
    node_1: Required[IntegratedNode]
    node_2: Required[IntegratedNode]


IntegratedGraphEntry = IntegratedSingleNodeEntry | IntegratedEdgeNodeEntry


class EnzyWizardIntegrateOutput(TypedDict):
    output_type: Required[Literal["enzywizard_integrate"]]
    overall_statistics: Required[IntegratedOverallStatistics]
    integrated_graph: Required[List[IntegratedGraphEntry]]

'''
'''