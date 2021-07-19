from database import db_session
from models import UnitCenter
import numpy as np



def query(corr_ungrouped):
	""" Database query for getting center data
        :param list corr_ungrouped: The list of correspondences
        
        returns: The dict containing unit_id as key & center data (in np.array) as value       
    """
	
	data = {}
	with db_session() as session:
		query = session.query(UnitCenter) \
		               .filter(UnitCenter.unit_id.in_(corr_ungrouped), \
		                  	      UnitCenter.name == 'base')

		for row in query:
			data[row.unit_id] = np.array([row.x, row.y, row.z])

		return data


def group_center(center, corr):
	""" Group center data according to ife and correspondence
        :param dict center: The dict containing unit_id as key & center data (in np.array) as value
        :param dict corr: The dict containing ife as key & list of correspondences as value
		                  after removing modified nts position

        returns: The dict containing ife as key & list of center data (in np.array) as value
    """
	
	center_grouped = {}
	for ife, correspondence in corr.iteritems():
		center_corr = []
		for unit in correspondence:
			center_data = center.get(unit, None)
			center_corr.append(center_data)
		
		center_grouped[ife] = center_corr	

	return center_grouped


def get_center(corr_ungrouped, corr_std):
	""" Main method that deals with getting center data
        :param list corr_ungrouped: The list of correspondences
        :param dict corr_std: The dict containing ife as key & list of correspondences as value
		                      after removing modified nts position
        
        returns: The dict containing ife as key & list of center data (in np.array) as value       
    """

	center_ungrouped = query(corr_ungrouped)
	center_grouped = group_center(center_ungrouped, corr_std)

	return center_grouped