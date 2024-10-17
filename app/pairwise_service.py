from database import db_session
from models import UnitPairsInteractions2024
from collections import OrderedDict
from sqlalchemy import tuple_


def get_pairwise_interactions(correspondences):

	with db_session() as session:
		pairwise_data = {}
		nt_pairs_positions = set()

	    # revise for Python 2/3 compatibility
		for ife, correspondence in correspondences.items():
			corr_len = len(correspondence)
			pairwise_interactions = {}
			for idx1 in range(0, corr_len):
				for idx2 in range(idx1 + 1, corr_len):
					pos1 = idx1 + 1
					pos2 = idx2 + 1
					key = str(pos1) + "--" + str(pos2)

					query = session.query(UnitPairsInteractions2024) \
					               .filter(UnitPairsInteractions2024.unit_id_1 == correspondence[idx1]) \
					               .filter(UnitPairsInteractions2024.unit_id_2 == correspondence[idx2])

					for row in query:
						if (row.f_lwbp_detail is not None) or (row.f_stacks is not None) or (row.f_brbs is not None) or (row.f_bphs is not None):
							interaction = list((row.f_lwbp_detail, row.f_stacks, row.f_brbs, row.f_bphs))
							interaction = ','.join(filter(None, interaction))
							pairwise_interactions[key] = interaction
							nt_pairs_positions.add(key)

			pairwise_data[ife] = pairwise_interactions

		return (pairwise_data, nt_pairs_positions)

def get_pairwise_interactions_new(nt_pairs_dict, nt_position_dict):
	with db_session() as session:
		pairwise_data = {}
		nt_pairs_positions = set()

	    # revise for Python 2/3 compatibility
		for ife, nt_pairs in nt_pairs_dict.items():
			pairwise_interactions = {}
			pdb = ife.split("|")[0]
			query = session.query(UnitPairsInteractions2024) \
						   .filter(tuple_(UnitPairsInteractions2024.unit_id_1, UnitPairsInteractions2024.unit_id_2).in_(nt_pairs)) \
						   .filter(UnitPairsInteractions2024.pdb_id == pdb)

			for row in query:
				if (row.f_lwbp_detail is not None) or (row.f_stacks is not None):
					pos1 = nt_position_dict.get(row.unit_id_1)
					pos2 = nt_position_dict.get(row.unit_id_2)
					key = str(pos1) + "--" + str(pos2)
					interaction = list((row.f_lwbp_detail, row.f_stacks))
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

				query = session.query(UnitPairsInteractions2024) \
					            .filter(UnitPairsInteractions2024.unit_id_1 == correspondence[idx1]) \
					            .filter(UnitPairsInteractions2024.unit_id_2 == correspondence[idx2])


				for row in query:
					if (row.f_lwbp_detail is not None) or (row.f_stacks is not None) or (row.f_brbs is not None) or (row.f_bphs is not None):
						interaction = list((row.f_lwbp_detail, row.f_stacks, row.f_brbs, row.f_bphs))
						interaction = ','.join(filter(None, interaction))
						pairwise_interactions[key] = interaction
						nt_pairs_positions.add(key)


		return pairwise_interactions, nt_pairs_positions