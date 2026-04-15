from typing import TypedDict, List, Dict, Literal, Required, NotRequired

CLEAN_TYPE = "enzywizard_clean"
MUTCLEAN_TYPE = "enzywizard_mut_clean"
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
MUT_INTEGRATE_TYPE = "enzywizard_mut_integrate"


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
    removed_unexpected_heavy_atoms: int
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
    output_type: Literal["enzywizard_mut_clean"]

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


class HydrophobicClusterStatistics(TypedDict):
    cluster_num: int
    max_cluster_area: float
    total_cluster_area: float


class EnzyWizardHydroClusterOutput(TypedDict):
    output_type: Literal["enzywizard_hydrocluster"]
    hydrophobic_cluster_statistics: HydrophobicClusterStatistics
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


class DisorderRegionStatistics(TypedDict):
    region_num: int
    max_region_length: int
    total_region_length: int


class EnzyWizardDisorderOutput(TypedDict):
    output_type: Literal["enzywizard_disorder"]
    disorder_region_statistics: DisorderRegionStatistics
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


class PocketRegionStatistics(TypedDict):
    pocket_num: int
    max_pocket_volume: float
    total_pocket_volume: float


class EnzyWizardPocketOutput(TypedDict):
    output_type: Literal["enzywizard_pocket"]
    pocket_region_statistics: PocketRegionStatistics
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

    hydrophobic_cluster_count: NotRequired[int]
    max_hydrophobic_cluster_area: NotRequired[float]
    total_hydrophobic_cluster_area: NotRequired[float]

    disorder_region_count: NotRequired[int]
    max_disorder_region_length: NotRequired[int]
    total_disorder_region_length: NotRequired[int]

    pocket_region_count: NotRequired[int]
    max_pocket_region_volume: NotRequired[float]
    total_pocket_region_volume: NotRequired[float]

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
    aa_name_one_hot: NotRequired[List[int]]
    aa_coord: NotRequired[List[float]]

    aa_class: NotRequired[str]
    aa_class_one_hot: NotRequired[List[int]]
    aa_ss: NotRequired[str]
    aa_ss_one_hot: NotRequired[List[int]]
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

'''
mut_integrate_report
'''
class MutIntegratedOverallStatistics(TypedDict):
    wt_aa_name_count: NotRequired[List[int]]
    mut_aa_name_count: NotRequired[List[int]]
    diff_aa_name_count: NotRequired[List[int]]

    wt_aa_class_count: NotRequired[List[int]]
    mut_aa_class_count: NotRequired[List[int]]
    diff_aa_class_count: NotRequired[List[int]]

    wt_aa_ss_count: NotRequired[List[int]]
    mut_aa_ss_count: NotRequired[List[int]]
    diff_aa_ss_count: NotRequired[List[int]]

    wt_hydrophobic_cluster_count: NotRequired[int]
    mut_hydrophobic_cluster_count: NotRequired[int]
    diff_hydrophobic_cluster_count: NotRequired[int]

    wt_max_hydrophobic_cluster_area: NotRequired[float]
    mut_max_hydrophobic_cluster_area: NotRequired[float]
    diff_max_hydrophobic_cluster_area: NotRequired[float]

    wt_total_hydrophobic_cluster_area: NotRequired[float]
    mut_total_hydrophobic_cluster_area: NotRequired[float]
    diff_total_hydrophobic_cluster_area: NotRequired[float]

    wt_disorder_region_count: NotRequired[int]
    mut_disorder_region_count: NotRequired[int]
    diff_disorder_region_count: NotRequired[int]

    wt_max_disorder_region_length: NotRequired[int]
    mut_max_disorder_region_length: NotRequired[int]
    diff_max_disorder_region_length: NotRequired[int]

    wt_total_disorder_region_length: NotRequired[int]
    mut_total_disorder_region_length: NotRequired[int]
    diff_total_disorder_region_length: NotRequired[int]

    wt_pocket_region_count: NotRequired[int]
    mut_pocket_region_count: NotRequired[int]
    diff_pocket_region_count: NotRequired[int]

    wt_max_pocket_region_volume: NotRequired[float]
    mut_max_pocket_region_volume: NotRequired[float]
    diff_max_pocket_region_volume: NotRequired[float]

    wt_total_pocket_region_volume: NotRequired[float]
    mut_total_pocket_region_volume: NotRequired[float]
    diff_total_pocket_region_volume: NotRequired[float]

    wt_total_potential_energy: NotRequired[float]
    mut_total_potential_energy: NotRequired[float]
    diff_total_potential_energy: NotRequired[float]

    wt_harmonic_bond_force: NotRequired[float]
    mut_harmonic_bond_force: NotRequired[float]
    diff_harmonic_bond_force: NotRequired[float]

    wt_harmonic_angle_force: NotRequired[float]
    mut_harmonic_angle_force: NotRequired[float]
    diff_harmonic_angle_force: NotRequired[float]

    wt_custom_bond_force: NotRequired[float]
    mut_custom_bond_force: NotRequired[float]
    diff_custom_bond_force: NotRequired[float]

    wt_custom_torsion_force: NotRequired[float]
    mut_custom_torsion_force: NotRequired[float]
    diff_custom_torsion_force: NotRequired[float]

    wt_custom_nonbonded_force: NotRequired[float]
    mut_custom_nonbonded_force: NotRequired[float]
    diff_custom_nonbonded_force: NotRequired[float]

    wt_nonbonded_force: NotRequired[float]
    mut_nonbonded_force: NotRequired[float]
    diff_nonbonded_force: NotRequired[float]

    wt_periodic_torsion_force: NotRequired[float]
    mut_periodic_torsion_force: NotRequired[float]
    diff_periodic_torsion_force: NotRequired[float]

    wt_cmap_torsion_force: NotRequired[float]
    mut_cmap_torsion_force: NotRequired[float]
    diff_cmap_torsion_force: NotRequired[float]

    wt_docking_score: NotRequired[float]
    mut_docking_score: NotRequired[float]
    diff_docking_score: NotRequired[float]

    wt_hbond_count: NotRequired[int]
    mut_hbond_count: NotRequired[int]
    diff_hbond_count: NotRequired[int]

    wt_ionic_count: NotRequired[int]
    mut_ionic_count: NotRequired[int]
    diff_ionic_count: NotRequired[int]

    wt_vdw_count: NotRequired[int]
    mut_vdw_count: NotRequired[int]
    diff_vdw_count: NotRequired[int]

    wt_pipistack_count: NotRequired[int]
    mut_pipistack_count: NotRequired[int]
    diff_pipistack_count: NotRequired[int]

    wt_pication_count: NotRequired[int]
    mut_pication_count: NotRequired[int]
    diff_pication_count: NotRequired[int]

    wt_ssbond_count: NotRequired[int]
    mut_ssbond_count: NotRequired[int]
    diff_ssbond_count: NotRequired[int]


class MutationSiteFeatures(TypedDict):
    wt_aa_name: NotRequired[str]
    mut_aa_name: NotRequired[str]

    wt_aa_class: NotRequired[str]
    mut_aa_class: NotRequired[str]

    wt_aa_ss: NotRequired[str]
    mut_aa_ss: NotRequired[str]

    wt_aa_rsa: NotRequired[float]
    mut_aa_rsa: NotRequired[float]
    diff_aa_rsa: NotRequired[float]

    wt_aa_phi: NotRequired[float]
    mut_aa_phi: NotRequired[float]
    diff_aa_phi: NotRequired[float]

    wt_aa_psi: NotRequired[float]
    mut_aa_psi: NotRequired[float]
    diff_aa_psi: NotRequired[float]

    wt_aa_net_charge: NotRequired[float]
    mut_aa_net_charge: NotRequired[float]
    diff_aa_net_charge: NotRequired[float]

    wt_aa_pka: NotRequired[float]
    mut_aa_pka: NotRequired[float]
    diff_aa_pka: NotRequired[float]

    wt_aa_volume: NotRequired[float]
    mut_aa_volume: NotRequired[float]
    diff_aa_volume: NotRequired[float]

    wt_aa_hydrophobicity: NotRequired[float]
    mut_aa_hydrophobicity: NotRequired[float]
    diff_aa_hydrophobicity: NotRequired[float]

    wt_aa_molecular_weight: NotRequired[float]
    mut_aa_molecular_weight: NotRequired[float]
    diff_aa_molecular_weight: NotRequired[float]

    wt_aa_pi: NotRequired[float]
    mut_aa_pi: NotRequired[float]
    diff_aa_pi: NotRequired[float]

    wt_rmsf: NotRequired[float]
    mut_rmsf: NotRequired[float]
    diff_rmsf: NotRequired[float]

    wt_conservation_score: NotRequired[float]
    mut_conservation_score: NotRequired[float]
    diff_conservation_score: NotRequired[float]

    wt_aa_name_one_hot: NotRequired[List[int] | List[float]]
    mut_aa_name_one_hot: NotRequired[List[int] | List[float]]
    diff_aa_name_one_hot: NotRequired[List[int] | List[float]]

    wt_aa_class_one_hot: NotRequired[List[int] | List[float]]
    mut_aa_class_one_hot: NotRequired[List[int] | List[float]]
    diff_aa_class_one_hot: NotRequired[List[int] | List[float]]

    wt_aa_ss_one_hot: NotRequired[List[int] | List[float]]
    mut_aa_ss_one_hot: NotRequired[List[int] | List[float]]
    diff_aa_ss_one_hot: NotRequired[List[int] | List[float]]


class EnzyWizardMutIntegrateOutput(TypedDict):
    output_type: Required[Literal["enzywizard_mut_integrate"]]
    cleaned_amino_acid_substitution: Required[str]

    overall_statistics: Required[MutIntegratedOverallStatistics]
    mutation_site_features: Required[MutationSiteFeatures]

    wt_integrated_graph: Required[List[IntegratedGraphEntry]]
    mut_integrated_graph: Required[List[IntegratedGraphEntry]]

'''
'''