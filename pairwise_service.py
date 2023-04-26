from database import db_session
from models import UnitPairsInteractions
from collections import OrderedDict
from sqlalchemy import tuple_


def get_pairwise_interactions(correspondences):

	with db_session() as session:
		pairwise_data = {}
		nt_pairs_positions = set()

		for ife, correspondence in correspondences.iteritems():
			corr_len = len(correspondence)
			pairwise_interactions = {}
			for idx1 in range(0, corr_len):
				for idx2 in range(idx1 + 1, corr_len):
					pos1 = idx1 + 1
					pos2 = idx2 + 1
					key = str(pos1) + "--" + str(pos2)

					query = session.query(UnitPairsInteractions) \
					               .filter(UnitPairsInteractions.unit_id_1 == correspondence[idx1]) \
					               .filter(UnitPairsInteractions.unit_id_2 == correspondence[idx2])

					for row in query:
						if (row.f_lwbp is not None) or (row.f_stacks is not None) or (row.f_brbs is not None) or (row.f_bphs is not None):
							interaction = list((row.f_lwbp, row.f_stacks, row.f_brbs, row.f_bphs))
							interaction = ','.join(filter(None, interaction))
							pairwise_interactions[key] = interaction
							nt_pairs_positions.add(key)

			pairwise_data[ife] = pairwise_interactions
		
		return (pairwise_data, nt_pairs_positions)

def get_pairwise_interactions_new(nt_pairs_dict, nt_position_dict):
	with db_session() as session:
		pairwise_data = {}
		nt_pairs_positions = set()

		for ife, nt_pairs in nt_pairs_dict.iteritems():
			pairwise_interactions = {}
			pdb = ife.split("|")[0]
			query = session.query(UnitPairsInteractions) \
						   .filter(tuple_(UnitPairsInteractions.unit_id_1, UnitPairsInteractions.unit_id_2).in_(nt_pairs)) \
						   .filter(UnitPairsInteractions.pdb_id == pdb)

			for row in query:
				if (row.f_lwbp is not None) or (row.f_stacks is not None):
					pos1 = nt_position_dict.get(row.unit_id_1)
					pos2 = nt_position_dict.get(row.unit_id_2)
					key = str(pos1) + "--" + str(pos2)
					interaction = list((row.f_lwbp, row.f_stacks))
					interaction = ','.join(filter(None, interaction))
					pairwise_interactions[key] = interaction
					nt_pairs_positions.add(key)
					
			pairwise_data[ife] = pairwise_interactions

	return (pairwise_data, nt_pairs_positions)


def get_pairwise_interactions_single(correspondence):

	with db_session() as session:

		nt_pairs_positions = set()
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
						nt_pairs_positions.add(key)


		return pairwise_interactions, nt_pairs_positions