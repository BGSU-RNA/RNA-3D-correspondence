from database import db_session
from models import NrChains, NrClasses, NrReleases, IfeInfo, PDBInfo, ChainInfo 

REJECT_LIST = ['5LZE|1|a+5LZE|1|y+5LZE|1|v+5LZE|1|x']

def get_source_organism(chain_id):

	pdb_id, model, chain = chain_id.split('|')

	with db_session() as session:
		query = session.query(ChainInfo).filter_by(pdb_id=pdb_id) \
										.filter_by(chain_name=chain)
		return query[0].source


def get_class_id(ife, resolution):

	with db_session() as session:
		query = session.query(NrChains).join(NrClasses, NrReleases) \
	        					       .filter(NrChains.ife_id == ife) \
	        					       .filter(NrClasses.resolution == resolution) \
	                                   .order_by(NrReleases.date.desc()) \
	                                   .limit(1)
		return query[0].nr_class_id


def get_members(class_id, exp_method):

	if exp_method != 'all':
	
		with db_session() as session:
			query = session.query(NrChains) \
						   .join(IfeInfo, NrChains.ife_id == IfeInfo.ife_id) \
						   .join(PDBInfo, IfeInfo.pdb_id == PDBInfo.pdb_id)	\
			               .filter(NrChains.nr_class_id == class_id) \
			               .filter(PDBInfo.experimental_technique == exp_method)
			return [row.ife_id for row in query]

	else:

		with db_session() as session:
			query = session.query(NrChains).filter(NrChains.nr_class_id == class_id)
			return [row.ife_id for row in query]


def get_ec_info(class_id):
	
	with db_session() as session:
		query = session.query(NrClasses).filter_by(nr_class_id=class_id)
		return query[0].name, query[0].nr_release_id
	


def remove_ife(members, query_ife):

	REJECT_LIST.append(query_ife)

	for ife in REJECT_LIST:
		members.remove(ife)
	
	return members


# Potential bug. Consider records with only 1 chain
def get_chains(members, index_val):

	chains = []

	for elem in members:
		if '+' in elem:
			chains.append(elem.split('+')[index_val])
		else:
			chains.append(elem)

	return chains


def process_members(members, id_type, index_val):

	if id_type == 'ife':
		# remove entries with only single chain
		members = [elem for elem in members if '+' in elem]
		members = get_chains(members, index_val)

	members_processed = []
	for ife in members:
		pdb = ife.split('|')[0]
		chain = ife.split('|')[-1]
		members_processed.append((pdb, chain))

	return members_processed


def get_complete_id(id):

	#search = "{}%".format(id)
	search = str(id) + "%"

	with db_session() as session:
		query = session.query(IfeInfo).filter(IfeInfo.ife_id.like(search))
		return query[0].ife_id


def get_exp_method(pdbid):

	with db_session() as session:
		query = session.query(PDBInfo).filter(PDBInfo.pdb_id == pdbid)
		return query[0].experimental_technique
	

def check_chain_or_ife(id):

	if '+' in id:
		return 'ife'
	else:
		return 'chain'


def get_id_index(query_id, complete_id, id_type):

	if id_type == 'ife':
		all_chains = complete_id.split('+')
		return all_chains.index(query_id)
	else:
		return None


def get_id_type(query_id):

	if query_id != 'unitid':
		complete_id = get_complete_id(query_id)
		id_type = check_chain_or_ife(complete_id)
		query_id_index = get_id_index(query_id, complete_id, id_type)
	return complete_id, id_type, query_id_index


def get_ec_members(resolution, exp_method, query_id):
	
	#query_ife = '|'.join(units[0].split('|')[:3])

	error_msg = "No members found in equivalence class"

	try:
		error_msg = "Query id type not found"

		print("get_ec_members query_id:",query_id)

		id, id_type, query_id_index = get_id_type(query_id)

		print("get_ec_members id, id_type, query_id_index:",id,id_type,query_id_index)

		error_msg = "Equivalence class id not found, check resolution threshold"
		class_id = get_class_id(id, resolution)
		error_msg = "Equivalence class name/release not found"
		ec_name, nr_release = get_ec_info(class_id)
		error_msg = "Equivalence class members not found"
		members = get_members(class_id, exp_method)
		if id in members: members.remove(id)

		members = process_members(members, id_type, query_id_index)

		if len(members) == 0:
			error_msg = "No members found in equivalence class"
		else:
			error_msg = ""

		return members, ec_name, nr_release, error_msg

	except:
		return [], "", "", error_msg


def check_valid_membership(members, query_data, exp_method):
	
	pdbid = query_data[0]
	selection_method = get_exp_method(pdbid)

	if not members:
		empty_members = True
	else:
		empty_members = False

	if selection_method == exp_method:
		method_equality = True
	else:
		method_equality = False

	# Check this logic
	if exp_method == "all":
		method_equality = True

	return empty_members, method_equality


def get_pdb_resolution(pdb_list):
	
	resolution_dict = {}
	with db_session() as session:
		query = session.query(PDBInfo).filter(PDBInfo.pdb_id.in_(pdb_list))

		for row in query:
			resolution_dict[row.pdb_id] = '{0:.2f}'.format(row.resolution)

		return resolution_dict




