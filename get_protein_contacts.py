from database import db_session
from models import UnitPairsDistances, UnitInfo, ChainInfo
from sqlalchemy import tuple_

test_units = ['5J7L|1|AA|A|964', '5J7L|1|AA|A|968', '5J7L|1|AA|A|969', '5J7L|1|AA|C|970', '5J7L|1|AA|C|972', '5J7L|1|AA|G|963', '5J7L|1|AA|G|971', '5J7L|1|AA|U|965']

test_chain = '5J7L|1|AA|'
test_chain = '5UYM|1|01|'

def get_protein_contacts(units):
    """A function that returns the protein residues that are located
       within 8.5 Angstroms of the loop nts

    Args:
        units (list of strings): A list containing the unitid of loop nts

    Returns:
        list of tuples: A list of 3-elements tuple in the form (rna_unit_id, protein_unit_id, distance)
    """

    contacts_list = []
    with db_session() as session:
        query = session.query(UnitPairsDistances) \
                       .join(UnitInfo, UnitPairsDistances.unit_id_2 == UnitInfo.unit_id) \
                       .filter(UnitInfo.unit_type_id == 'aa') \
                       .filter(UnitPairsDistances.distance <= 8.5) \
					   .filter(UnitPairsDistances.unit_id_1.in_(units))

        return [(row.unit_id_1, row.unit_id_2, row.distance) for row in query]

def get_bound_protein_chains(rna_chain):
    """A function that returns all the protein chains bound to a particular RNA chain

    Args:
        rna_chain (string): RNA chain

    Returns:
        set of tuples: A set of 2-elements tuple in the form of (pdbid, chain)
    """

    bound_protein_chains = set()
    search_keyword = rna_chain + "%"
    with db_session() as session:
        query = session.query(UnitPairsDistances) \
                       .join(UnitInfo, UnitPairsDistances.unit_id_2 == UnitInfo.unit_id) \
                       .filter(UnitPairsDistances.unit_id_1.like(search_keyword)) \
                       .filter(UnitInfo.unit_type_id == 'aa') \
                       .filter(UnitPairsDistances.distance <= 8.5) 
    
        for row in query:
            protein_info_list = row.unit_id_2.split("|")
            pdb, _, chain, _ = protein_info_list[0], protein_info_list[1], protein_info_list[2], protein_info_list[3:]
            protein_info = (pdb, chain)
            bound_protein_chains.add(protein_info)

        return bound_protein_chains

def get_protein_names(protein_info):
    """A function that returns the names of protein chains

    Args:
        protein_info (set of tuples): A set of 2-elements tuple in the form of (pdbid, chain)

    Returns:
        dictionary: The key is the chain_id while the value is the compound name
    """

    chain_names = {}
    with db_session() as session:
        query = session.query(ChainInfo) \
                       .filter(tuple_(ChainInfo.pdb_id, ChainInfo.chain_name) \
					   .in_(protein_info))

        for row in query:
            chain_names[row.chain_name] = row.compound

        return chain_names

#loop_protein_contacts = get_protein_contacts(test_units)

protein_chains = get_bound_protein_chains(test_chain)
protein_chains = list(protein_chains)

protein_names = get_protein_names(protein_chains)

for k, v in protein_names.iteritems():
    print k + ": " + v
