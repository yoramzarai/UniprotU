"""
Few UniprotU usage example.

[This file is not included in the build .whl file]
"""
import uniprotu as u  # (__init__.py contain "from uniprotu.uniprot_utils import *")

def main() -> None:
    """Testing few uniprotu utils."""
    print("Ensembl ID to UniProt ID:")
    ensmbl_id: str = "ENST00000559488"
    uniprot_id = u.ensembl_id2uniprot_id(ensmbl_id)
    print(f"{ensmbl_id} => {uniprot_id}.")

    print("\nDomains:")
    # Set the list of features to retrieve. Set to [] to retrieve all features.
    uniprot_features: list[str] = ["Topological domain", "Transmembrane", "Domain", "Repeat", "Region"]
    df_domains = u.retrieve_protein_data_features_subset(uniprot_id, uniprot_features)
    print(df_domains.to_string())

    print("\nAA sequence:")
    print(u.AA_seq(uniprot_id))

    print("\nCross-reference:")
    df_cross = u.get_CrossReferences_databases_info(uniprot_id)
    print(df_cross.to_string(max_rows=10))

    print("\nGet all protein data:")
    protein_data = u.lookup_protein_data(uniprot_id)
    print(f"protein_data contains {len(protein_data):,} entries...")

    print("\nGet protein features:")
    df_features = u.retrieve_protein_data_features_subset(uniprot_id, [])
    print(df_features.to_string(max_rows=10))

    print("\nUniProt ID to Ensembl ID:")
    uid: str = 'P42336'
    eid = u.uniprot_id2ensembl_id(uid)
    print(f"{uid} => {eid}")

if __name__ == '__main__':
    main()
