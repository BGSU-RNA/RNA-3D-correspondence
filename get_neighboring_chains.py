from contextlib import contextmanager
from models import UnitCenter, ChainInfo
from database import db_session
from sqlalchemy import tuple_
import math
import timeit

def filter_neighboring_residues(centers_coord, potential_neighboring_units, distance, nt_ids):
    # use the x,y,z coordinates in centers_coord to check if x,y,z coordinates in query results
    # potential_neighboring_units are within distance, avoiding those in nt_ids

    output_nt_ids = []
    #output_distance_list = []
    distance_squared = distance * distance

    for unit_arr in potential_neighboring_units:
        # if unit id of this potential unit is in the query, don't check distances
        if unit_arr['unit_id'] not in nt_ids:
            # if the unit id of this potential unit is already in the output, don't check distances
            # That misses the possibility of finding an even closer match with a different center.
            if unit_arr['unit_id'] not in output_nt_ids:

                # keep track of minimum distance of potential unit unit_arr to query x,y,z locations
                d2min = 10*distance_squared

                # loop over query x,y,z locations
                # for i in range(len(centers_coord[0])):
                for i, _ in enumerate(centers_coord[0]):
                    # calculate squared distance, keep track of minimum squared distance
                    d2 = math.pow((unit_arr['x'] - centers_coord[0][i]), 2)  + math.pow((unit_arr['y'] - centers_coord[1][i]), 2) + math.pow((unit_arr['z'] - centers_coord[2][i]), 2)
                    if d2 < d2min:
                        d2min = d2

                if d2min < distance_squared:
                    output_nt_ids.append(str(unit_arr['unit_id']))
                    # currently not returning output_distance_list, so don't compute it
                    #output_distance_list.append(math.sqrt(d2min))

    # currently not returning output_distance_list

    return output_nt_ids


def get_xyz_coordinates_between_limits(pdb_id, model_num, coord_limits):
    # Find all units in the given PDB file whose center is within the given limits.
    # In the future, if we want to, we could
    # only find residues where the base center and the amino acid functional group center
    # is within the limits
    #center_type = ['base', 'aa_fg']

    with db_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.x >= coord_limits[0]) \
                       .filter(UnitCenter.x <= coord_limits[1]) \
                       .filter(UnitCenter.y >= coord_limits[2]) \
                       .filter(UnitCenter.y <= coord_limits[3]) \
                       .filter(UnitCenter.z >= coord_limits[4]) \
                       .filter(UnitCenter.z <= coord_limits[5]) 

        unit_coord_arr = []
        for row in query:
            mn = row.unit_id.split("|")
            if mn[1] == model_num: # neighbors from same model
                unit_coord = {
                    "unit_id": row.unit_id,
                    "x": float(row.x),
                    "y": float(row.y),
                    "z": float(row.z),
                    "name": row.name
                }
                unit_coord_arr.append(unit_coord)
        return unit_coord_arr


def get_xyz_coordinates(unit_ids, pdb_id):
    # retrieve the x,y,z coordinates of all centers of all units in unit_ids
    # for nucleotides, that will include the base center, glycosidic atom, sugar center, phosphate center
    # for amino acids, that will include the functional group center and backbone center

    with db_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.unit_id.in_(unit_ids))
        
        if not query:
            return False

        given_x = []
        given_y = []
        given_z = []
        for row in query:
            given_x.append(float(row.x))
            given_y.append(float(row.y))
            given_z.append(float(row.z))
        
        centers_coord = [given_x, given_y, given_z]
        return centers_coord

def get_possible_chains_list(neighboring_residues, query_chain):
    """
    This function takes in a list of neighboring residues and the query chain
    and returns a list of unique tuples containing the PDB ID and chain id 
    of all the neighboring residues except from the input chain.

    Input:
    neighboring_residues: list of neighboring unitids
    chain: string, the query chain to be excluded from the output list

    Output:
    A list of tuples containing the PDBID and chain id of the neighboring 
    chains. If the input list is empty, the function returns None.
    """
    if neighboring_residues:
        # get a list of tuples where the first element is the pdb id while the second element is the chain id
        return list(set([ (unit.split("|")[0], unit.split("|")[2]) for unit in neighboring_residues if unit.split("|")[2] != query_chain]))
    else:
        return None

def get_chain_name(chains_list):
    """
    This function takes in a list of chains and returns a dictionary containing the chain names and their 
    corresponding compounds.

    Input:
    chains_list: list of tuples, where each tuple contains a PDB ID and a chain name

    Output:
    chain_name_dict: a dictionary where the key is the chain name and the value is the corresponding compound 
    name. If the input list is empty, the function returns None.
    """
    if chains_list:
        with db_session() as session:
            chain_data_list = []
            for chain in chains_list:
                pdb_id = chain[0]
                chain_name = chain [1]
                query = session.query(ChainInfo) \
                               .filter(ChainInfo.pdb_id == pdb_id) \
                               .filter(ChainInfo.chain_name == chain_name) \
                
                for row in query:
                    chain_data = (str(row.chain_name), str(row.compound))
                    chain_data_list.append(chain_data)
               
            return chain_data_list
    else:
        return None

def get_chains(unit_ids, distance=10):
    # Starting with unit_ids, find the x,y,z coordinates of their centers,
    # expand by distance to a rectangular box around them,
    # then find other units within that box,
    # then filter down to ones that are within distance of one of the centers in unit_ids

    # get the pdb id and model number from the first unit id
    unit_ids = unit_ids.split(",")
    fields = unit_ids[0].split("|")
    pdb_id = fields[0]
    model_num = fields[1]
    chain = fields[2]

    # Get all centers of unit_ids, including base, sugar, phosphate, aa_fg
    centers_xyz_coord = get_xyz_coordinates(unit_ids, pdb_id)

    # # Find the maxima and minima and expand by distance
    x_min = min(centers_xyz_coord[0]) - distance
    x_max = max(centers_xyz_coord[0]) + distance
    y_min = min(centers_xyz_coord[1]) - distance
    y_max = max(centers_xyz_coord[1]) + distance
    z_min = min(centers_xyz_coord[2]) - distance
    z_max = max(centers_xyz_coord[2]) + distance

    # store the limits in an array
    coord_limits = [x_min, x_max, y_min, y_max, z_min, z_max]

    # query to find all units whose x, y, z coordinates are between the limits
    potential_neighboring_units = get_xyz_coordinates_between_limits(pdb_id, model_num, coord_limits)

    # # the following line is unlikely to happen
    if len(potential_neighboring_units) == 0:
        return []

    # skipping units already in unit_ids, check distances of potential to centers_xyz_coord
    # calculate distances to units in unit_ids, record the smallest
    neighboring_residues = filter_neighboring_residues(centers_xyz_coord, potential_neighboring_units, distance, unit_ids)
    neighboring_chains = get_possible_chains_list(neighboring_residues, chain)
    neighboring_chains_with_names = get_chain_name(neighboring_chains)

    # return neighboring_residues
    return neighboring_chains_with_names


def test_run(units_list):
    # start_time = timeit.default_timer()
    # with terminating(Pool(processes=4)) as pool:
    #     result = pool.map(get_neighboring_chains, units_list)
    #     return result
    result = [get_chains(x) for x in units_list]
    # elapsed = timeit.default_timer() - start_time
    return result


# units = "5J7L|1|AA|G|1491,5J7L|1|AA|A|1492,5J7L|1|AA|A|1493,5J7L|1|AA|G|1494,5J7L|1|AA|U|1495,5J7L|1|AA|C|1496"
# test_units = ["5J7L|1|AA|G|1491,5J7L|1|AA|A|1492,5J7L|1|AA|A|1493,5J7L|1|AA|G|1494,5J7L|1|AA|U|1495,5J7L|1|AA|C|1496", "5J7L|1|AA|A|1492,5J7L|1|AA|A|1493"]
# test_units = [units] * 200
# # print(len(test_units))
# result, result_time = test_run(test_units)
# print result_time
