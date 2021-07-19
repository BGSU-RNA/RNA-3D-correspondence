from database import db_session
from models import UnitCorrespondence, UnitInfo, LoopInfo
from sqlalchemy import case, tuple_
from itertools import groupby
import copy

def query(units, members):
	""" Database query for getting correspondences
		:param list units: The list of query unit-ids.
		:param list of tuples members: (pdb_id, chain) for each instance
		
		returns: The list of correspondences        
	"""

	ordering = case (
		#{unit: index for index, unit in enumerate(units)},
		dict((unit, index) for index, unit in enumerate(units)),
		value=UnitCorrespondence.unit_id_1
	)


	with db_session() as session:
		query = session.query(UnitCorrespondence) \
					   .filter(UnitCorrespondence.unit_id_1.in_(units)) \
					   .order_by(ordering) \
					   .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
					   .in_(members))
		
		return [row.unit_id_2 for row in query]


def get_loop_correspondence(loop_names):

	with db_session() as session:
		query = session.query(LoopInfo) \
					   .filter(tuple_(LoopInfo.pdb_id, LoopInfo.loop_name) \
					   .in_(loop_names))


		return [(row.loop_id, row.loop_name) for row in query]	


def group_correspondence(correspondence):
	""" Group correspondence according to ifes
		:param list correspondence: The list of correspondences
		
		returns: The dict containing ife as key & list of correspondences as value
	"""

	key_ife = lambda x: '|'.join(x.split('|')[:3])
	#sort according to ife
	correspondence.sort(key=key_ife)
	#return {grp: list(items) for grp, items in groupby(correspondence, key_ife)}
	return dict((grp, list(items)) for grp, items in groupby(correspondence, key_ife))


def check_modifications(correspondence_grouped):
	""" Check for the presence of modified nts and remove their corresponding positions
		:param dict correspondence_grouped: The dict containing ife as key & list of correspondences as value
		
		returns
		:dict correspondence_std_grouped: The dict containing ife as key & list of correspondences as value
										  after removing modified nts position	        

	"""

	accepted_sequences = ['A', 'C', 'G', 'U']
	modified_indices = []
	correspondence_std_grouped = copy.deepcopy(correspondence_grouped)

	for ife, correspondence in correspondence_grouped.iteritems():
		for unit in correspondence:
			seq = unit.split('|')[3]
			if seq not in accepted_sequences:
				modified_indices.append(correspondence.index(unit))

	nr_modified_indices = list(set(modified_indices))

	for ife, correspondence in correspondence_std_grouped.iteritems():
		for idx in sorted(nr_modified_indices, reverse=True):
			del correspondence[idx]

	return correspondence_std_grouped


def pairwise_structure_correspondence(ife1, ife2):

	pdb_id1, _, chain1 = ife1.split('|')
	pdb_id2, _, chain2 = ife2.split('|')
	
	
	with db_session() as session:
		query = session.query(UnitCorrespondence) \
					   .join(UnitInfo, UnitCorrespondence.unit_id_1 == UnitInfo.unit_id) \
					   .filter(UnitCorrespondence.pdb_id_1 == pdb_id1) \
					   .filter(UnitCorrespondence.chain_name_1 == chain1) \
					   .filter(UnitCorrespondence.pdb_id_2 == pdb_id2) \
					   .filter(UnitCorrespondence.chain_name_2 == chain2) \
					   .order_by(UnitInfo.chain_index.asc())

		return [(row.unit_id_1, row.unit_id_2) for row in query]

def get_correspondence(units, members, method_equality=True):
	""" Main method that deals with getting correspondences
		:param list units: The list of query unit-ids.
		:param list of tuples members: (pdb_id, chain) for each instance
		
		returns
		:list correspondence: The list of correspondence ungrouped
		:dict correspondence_grouped: The dict containing ife as key & list of correspondences as value
		:dict correspondence_std: The dict containing ife as key & list of correspondences as value
								  after removing modified nts position	        

	"""

	correspondence = query(units, members)	
	if method_equality is True: correspondence.extend(units)
	correspondence_grouped = group_correspondence(correspondence)
	correspondence_std_grouped = check_modifications(correspondence_grouped)

	return correspondence, correspondence_grouped, correspondence_std_grouped

