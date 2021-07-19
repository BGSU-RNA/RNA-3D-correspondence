from database import db_session
from models import UnitRotation
import numpy as np


def query(corr_ungrouped):
	""" Database query for getting rotation data
        :param list corr_ungrouped: The list of correspondences
        
        returns: The dict containing unit_id as key & rotation data (in np.array) as value       
    """

	data = {}

	with db_session() as session:
		query = session.query(UnitRotation) \
		               .filter(UnitRotation.unit_id.in_(corr_ungrouped))

		for row in query:
			data[row.unit_id] = np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2],
	                                            [row.cell_1_0, row.cell_1_1, row.cell_1_2],
	                                            [row.cell_2_0, row.cell_2_1, row.cell_2_2]])

		return data


def group_rotation(rotation, corr):
	""" Group rotation data according to ife and correspondence
        :param dict rotation: The dict containing unit_id as key & rotation data (in np.array) as value
        :param dict corr: The dict containing ife as key & list of correspondences as value
		                  after removing modified nts position

        returns: The dict containing ife as key & list of np_arrays rotation data as value
    """
        
	rotation_grouped = {}
	for ife, correspondence in corr.iteritems():
		rotation_corr = []
		for unit in correspondence:
			rotation_data = rotation.get(unit, None)
			rotation_corr.append(rotation_data)
		
		rotation_grouped[ife] = rotation_corr	

	return rotation_grouped


def get_rotation(corr_ungrouped, corr_std):
	""" Main method that deals with getting rotation data
        :param list corr_ungrouped: The list of correspondences
        :param dict corr_std: The dict containing ife as key & list of correspondences as value
		                      after removing modified nts position
        
        returns: The dictionary containing ife as key & list of rotation data (in np.array) as value        
    """
	
	rotation_ungrouped = query(corr_ungrouped)
	rotation_grouped = group_rotation(rotation_ungrouped, corr_std)

	return rotation_grouped