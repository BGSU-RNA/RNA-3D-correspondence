from database import db_session
from models import UnitPairsInteractions
from collections import OrderedDict


def get_pairwise_interactions(correspondence_data):

	with db_session() as session:

		pairwise_data = {}
		pairwise_residue_pairs_reference = set()

		for ife, correspondence in correspondence_data.iteritems():
			corr_len = len(correspondence)
			pairwise_interactions = {}
			for idx1 in range(0, corr_len):
				for idx2 in range(idx1 + 1, corr_len):
					pos1 = idx1 + 1
					pos2 = idx2 + 1
					key = "nt" + str(pos1) + "-nt" + str(pos2)

					query = session.query(UnitPairsInteractions) \
					               .filter(UnitPairsInteractions.unit_id_1 == correspondence[idx1]) \
					               .filter(UnitPairsInteractions.unit_id_2 == correspondence[idx2])


					for row in query:
						if (row.f_lwbp is not None) or (row.f_stacks is not None) or (row.f_brbs is not None) or (row.f_bphs is not None):
							interaction = list((row.f_lwbp, row.f_stacks, row.f_brbs, row.f_bphs))
							interaction = ','.join(filter(None, interaction))
							pairwise_interactions[key] = interaction
							pairwise_residue_pairs_reference.add(key)


			pairwise_data[ife] = pairwise_interactions

		'''
		pairwise_residue_pairs = list(pairwise_residue_pairs_reference)
		pairwise_residue_pairs.sort()

		chain_list = pairwise_data.keys()

		pairwise_interactions_collection = {}

		for chain in chain_list:
			pairwise_interactions_ordered = OrderedDict()
			for unique_res_pair in pairwise_residue_pairs:
				pairwise_interactions_ordered[unique_res_pair] = ''

			pairwise_interactions_collection[chain] = pairwise_interactions_ordered

		
		for chain in chain_list:
			for res_pair in pairwise_residue_pairs:
				try:
					pairwise_interactions_collection[chain][res_pair] = pairwise_data.get(chain, {}).get(res_pair)
				except:
					pass

		display_str = ''

		for k, v in pairwise_interactions_collection.iteritems():
			display_str = display_str + ','.join(str(x) for x in v.values()) + "</br>"
		'''
		
		return pairwise_data, pairwise_residue_pairs_reference


def get_pairwise_interactions_single(correspondence):

	with db_session() as session:

		pairwise_residue_pairs_reference = set()
		corr_len = len(correspondence)
		pairwise_interactions = {}
		for idx1 in range(0, corr_len):
			for idx2 in range(idx1 + 1, corr_len):
				pos1 = idx1 + 1
				pos2 = idx2 + 1
				#key = "Nt" + str(pos1) + "-Nt" + str(pos2)
				key = (int(pos1), int(pos2))

				query = session.query(UnitPairsInteractions) \
					            .filter(UnitPairsInteractions.unit_id_1 == correspondence[idx1]) \
					            .filter(UnitPairsInteractions.unit_id_2 == correspondence[idx2])


				for row in query:
					if (row.f_lwbp is not None) or (row.f_stacks is not None) or (row.f_brbs is not None) or (row.f_bphs is not None):
						interaction = list((row.f_lwbp, row.f_stacks, row.f_brbs, row.f_bphs))
						interaction = ','.join(filter(None, interaction))
						pairwise_interactions[key] = interaction
						pairwise_residue_pairs_reference.add(key)


		return pairwise_interactions, pairwise_residue_pairs_reference