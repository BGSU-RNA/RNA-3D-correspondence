from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
# from sqlalchemy import tuple_
from sqlalchemy.orm import aliased
import logging
import math

from models import UnitCenter, ChainInfo, ChainPropertyValue
from database import db_session

# List of substrings to check for
SUBSTRINGS_TO_CHECK = ['mrna', 'messenger rna', 'messenger-rna']

@contextmanager
def create_session():
    """
    Context manager to create a new session for each thread.
    This ensures that the database session is thread-safe.
    """
    with db_session() as session:
        yield session

def filter_neighboring_residues(centers_coord, potential_neighboring_units, distance, nt_ids):
    # use the x,y,z coordinates in centers_coord to check if x,y,z coordinates in query results
    # potential_neighboring_units are within distance, avoiding those in nt_ids

    output_nt_ids = []
    #output_distance_list = []
    distance_squared = distance * distance

    for unit in potential_neighboring_units:
        # if unit id of this potential unit is in the query, don't check distances
        if unit['unit_id'] not in nt_ids:
            # if the unit id of this potential unit is already in the output, don't check distances
            # That misses the possibility of finding an even closer match with a different center.
            if unit['unit_id'] not in output_nt_ids:

                # keep track of minimum distance of potential unit unit to query x,y,z locations
                d2min = 10*distance_squared

                # loop over query x,y,z locations
                # for i in range(len(centers_coord[0])):
                for i, _ in enumerate(centers_coord[0]):
                    # calculate squared distance
                    d2 = math.pow((unit['x'] - centers_coord[0][i]), 2)  + math.pow((unit['y'] - centers_coord[1][i]), 2) + math.pow((unit['z'] - centers_coord[2][i]), 2)
                    # keep track of minimum squared distance
                    if d2 < d2min:
                        d2min = d2

                if d2min < distance_squared:
                    output_nt_ids.append(str(unit['unit_id']))
                    # currently not returning output_distance_list, so don't compute it
                    #output_distance_list.append(math.sqrt(d2min))

    # currently not returning output_distance_list

    return output_nt_ids


def filter_neighboring_residues_unique_chain(centers_coord, potential_neighboring_units, distance):
    # use the x,y,z coordinates in centers_coord to check if x,y,z coordinates in query results
    # potential_neighboring_units are within distance
    # once a neighbor from a chain is found, do not check any more units from that chain

    output_nt_ids = []
    chains_found = set()
    distance_squared = distance * distance

    # loop over potential neighboring units
    for unit in potential_neighboring_units:
        fields = unit['unit_id'].split("|")
        chain = fields[2]
        # once a neighbor in a chain is found, do not check any more units from that chain
        if not chain in chains_found:
            x = unit['x']
            y = unit['y']
            z = unit['z']
            # x = float(unit['x'])
            # y = float(unit['y'])
            # z = float(unit['z'])
            # loop over query x,y,z locations
            for i, _ in enumerate(centers_coord[0]):
                # calculate squared distance between the units
                d2 = math.pow((x - centers_coord[0][i]), 2)  + math.pow((y - centers_coord[1][i]), 2) + math.pow((z - centers_coord[2][i]), 2)
                # compare to target squared distance
                if d2 < distance_squared:
                    # record the unit id
                    output_nt_ids.append(str(unit['unit_id']))
                    # record that at least one unit from this chain is close enough
                    chains_found.add(chain)
                    # no need to check distances of this unit to any other query x,y,z locations
                    break

    return output_nt_ids


def filter_neighboring_residues_unique_chain_tuple(centers_coord, potential_neighboring_units, distance):
    # use the x,y,z coordinates in centers_coord to check if x,y,z coordinates in query results
    # potential_neighboring_units are within distance
    # once a neighbor from a chain is found, do not check any more units from that chain

    output_nt_ids = []
    chains_found = set()
    distance_squared = distance * distance

    # loop over potential neighboring units
    for unit in potential_neighboring_units:
        fields = unit[0].split("|")
        chain = fields[2]
        # once a neighbor in a chain is found, do not check any more units from that chain
        if not chain in chains_found:
            x = unit[1]
            y = unit[2]
            z = unit[3]
            # x = float(unit['x'])
            # y = float(unit['y'])
            # z = float(unit['z'])
            # loop over query x,y,z locations
            for i, _ in enumerate(centers_coord[0]):
                # calculate squared distance between the units
                d2 = math.pow((x - centers_coord[0][i]), 2)  + math.pow((y - centers_coord[1][i]), 2) + math.pow((z - centers_coord[2][i]), 2)
                # compare to target squared distance
                if d2 < distance_squared:
                    # record the unit id
                    output_nt_ids.append(str(unit[0]))
                    # record that at least one unit from this chain is close enough
                    chains_found.add(chain)
                    # no need to check distances of this unit to any other query x,y,z locations
                    break

    return output_nt_ids


def get_xyz_coordinates_between_limits(pdb_id, model_num, coord_limits):
    # Find all units in the given PDB file whose center is within the given limits.
    # In the future, if we want to, we could
    # only find residues where the base center and the amino acid functional group center
    # is within the limits
    #center_type = ['base', 'aa_fg']

    with create_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.x >= coord_limits[0]) \
                       .filter(UnitCenter.x <= coord_limits[1]) \
                       .filter(UnitCenter.y >= coord_limits[2]) \
                       .filter(UnitCenter.y <= coord_limits[3]) \
                       .filter(UnitCenter.z >= coord_limits[4]) \
                       .filter(UnitCenter.z <= coord_limits[5])

        units = []
        for row in query:
            mn = row.unit_id.split("|")
            if mn[1] == model_num: # neighbors from same model
                unit_coord = {
                    "unit_id": row.unit_id,
                    # "x": float(row.x),
                    # "y": float(row.y),
                    # "z": float(row.z)
                    "x": row.x,
                    "y": row.y,
                    "z": row.z
                }
                units.append(unit_coord)
        return units


def get_xyz_coordinates_between_limits_different_chain(pdb_id, model_num, coord_limits, query_chains):
    # Find all units in the given PDB file whose center is within the given limits.
    # Only keep those whose chain is different from the query
    # In the future, if we want to, we could
    # only find residues where the base center and the amino acid functional group center
    # is within the limits
    #center_type = ['base', 'aa_fg']

    avoid_chain = "%s|%s|%s" % (pdb_id, model_num, sorted(query_chains)[0])

    with create_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.x >= coord_limits[0]) \
                       .filter(UnitCenter.x <= coord_limits[1]) \
                       .filter(UnitCenter.y >= coord_limits[2]) \
                       .filter(UnitCenter.y <= coord_limits[3]) \
                       .filter(UnitCenter.z >= coord_limits[4]) \
                       .filter(UnitCenter.z <= coord_limits[5]) \
                       .filter(~UnitCenter.unit_id.startswith(avoid_chain)) \
                       .order_by(UnitCenter.name.asc())

        units = []
        for row in query:
            fields = row.unit_id.split("|")
            if fields[2] not in query_chains:    # different chain than the query
                if fields[1] == model_num:       # same model as the query
                    unit_coord = {
                        "unit_id": row.unit_id,
                        "x": row.x,              # already a float
                        "y": row.y,
                        "z": row.z
                        # "x": float(row.x),
                        # "y": float(row.y),
                        # "z": float(row.z)
                    }
                    units.append(unit_coord)
        return units


def get_xyz_coordinates_between_limits_different_chain_tuple(pdb_id, model_num, coord_limits, query_chains):
    # Find all units in the given PDB file whose center is within the given limits.
    # Only keep those whose chain is different from the query
    # In the future, if we want to, we could
    # only find residues where the base center and the amino acid functional group center
    # is within the limits
    #center_type = ['base', 'aa_fg']

    avoid_chain = "%s|%s|%s" % (pdb_id, model_num, sorted(query_chains)[0])

    with create_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.x >= coord_limits[0]) \
                       .filter(UnitCenter.x <= coord_limits[1]) \
                       .filter(UnitCenter.y >= coord_limits[2]) \
                       .filter(UnitCenter.y <= coord_limits[3]) \
                       .filter(UnitCenter.z >= coord_limits[4]) \
                       .filter(UnitCenter.z <= coord_limits[5]) \
                       .filter(~UnitCenter.unit_id.startswith(avoid_chain)) \
                       .order_by(UnitCenter.name.asc())

        units = []
        for row in query:
            fields = row.unit_id.split("|")
            if fields[2] not in query_chains:    # different chain than the query
                if fields[1] == model_num:       # same model as the query
                    units.append((row.unit_id, row.x, row.y, row.z))

        # arithmetic apparently takes too long and does not help
        # if len(units) > 10:
        #     xm = (coord_limits[0] + coord_limits[1]) / 2.0
        #     ym = (coord_limits[2] + coord_limits[3]) / 2.0
        #     zm = (coord_limits[4] + coord_limits[5]) / 2.0
        #     units = sorted(units, key=lambda x: abs(x[1]-xm)+abs(x[2]-ym)+abs(x[3]-zm))   # sort by distance from x center

        return units


def get_xyz_coordinates(unit_ids, pdb_id):
    # retrieve the x,y,z coordinates of all centers of all units in unit_ids
    # for nucleotides, that will include the base center, glycosidic atom, sugar center, phosphate center
    # for amino acids, that will include the functional group center and backbone center

    with create_session() as session:
        query = session.query(UnitCenter) \
                       .filter(UnitCenter.pdb_id == pdb_id) \
                       .filter(UnitCenter.unit_id.in_(unit_ids))

        if not query:
            return False

        given_x = []
        given_y = []
        given_z = []
        for row in query:
            # given_x.append(float(row.x))
            # given_y.append(float(row.y))
            # given_z.append(float(row.z))
            given_x.append(row.x)
            given_y.append(row.y)
            given_z.append(row.z)

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
        return []

# def get_chain_name(chains_list):
#     """
#     This function takes in a list of chains and returns a dictionary containing the chain names and their
#     corresponding compounds.

#     Input:
#     chains_list: list of tuples, where each tuple contains a PDB ID and a chain name

#     Output:
#     chain_name_dict: a dictionary where the key is the chain name and the value is the corresponding compound
#     name. If the input list is empty, the function returns None.
#     """
#     if chains_list:
#         with db_session() as session:
#             chain_data_list = []
#             for chain in chains_list:
#                 pdb_id = chain[0]
#                 chain_name = chain [1]
#                 query = session.query(ChainInfo) \
#                                .filter(ChainInfo.pdb_id == pdb_id) \
#                                .filter(ChainInfo.chain_name == chain_name) \

#                 for row in query:
#                     chain_data = (str(row.chain_name), str(row.compound))
#                     chain_data_list.append(chain_data)



#             return chain_data_list
#     else:
#         return None

def get_chain_name_revised(pdb_chain_list):
    # Create aliases for the models to differentiate them in the query
    cpv = aliased(ChainPropertyValue)
    ci = aliased(ChainInfo)

    if pdb_chain_list:
        with create_session() as session:
            chain_data_list = []

            for pdb_id, chain in pdb_chain_list:
                # Query the ChainPropertyValue table for 'unp_name' or 'standardized_name'
                q_cpv_unp_name = session.query(cpv.value.label('NAME'), cpv.chain.label('CHAIN')).\
                    filter(cpv.pdb_id == pdb_id, cpv.chain == chain, cpv.property.in_(['unp_name', 'standardized_name'])).\
                    first()

                if q_cpv_unp_name:
                    chain_id = str(q_cpv_unp_name.CHAIN)
                    chain_name = str(q_cpv_unp_name.NAME)
                    if ";" in chain_name:
                        chain_name = chain_name.split(";")[0]
                    chain_data_list.append((chain_id, chain_name))
                else:
                    # If 'unp_name' or 'standardized_name' is not found, get the name from the ChainInfo table
                    q_ci = session.query(ci.compound.label('NAME'), ci.chain_name.label('CHAIN')).\
                        filter(ci.pdb_id == pdb_id, ci.chain_name == chain).\
                        first()

                    if q_ci:
                        chain_id = str(q_ci.CHAIN)
                        chain_name = str(q_ci.NAME)
                        chain_data_list.append((chain_id, chain_name))

            return chain_data_list if chain_data_list else []
    else:
        return []


def contains_mrna_substring(lst):
    for tup in lst:
        # Convert the second element of the tuple to lowercase
        second_element = tup[1].lower()

        # Check if any substring from the substrings_to_check list is present in the second_element
        if any(substring in second_element for substring in SUBSTRINGS_TO_CHECK):
            return True

    # If no tuple satisfies the condition, return False
    return False

def replace_with_mRNA(lst):
    for index, tup in enumerate(lst):
        # Convert the second element of the tuple to lowercase
        second_element = tup[1].lower()

        # Check if any substring from the substrings_to_check list is present in the second_element
        if any(substring in second_element for substring in SUBSTRINGS_TO_CHECK):
            # Replace the second element with "mRNA"
            lst[index] = (tup[0], "mRNA")

    return lst


def make_xyz_cubes(x, y, z, distance):
    """
    Build cubes with side length equal to distance.
    Cubes are named by a rounded value of x,y,z.
    All 26 neighboring cubes are generated.
    Each cube maps to the list of points in the cube.
    Each neighboring cube maps back to the cubes it neighbors.
    """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    key_to_points = {}
    key_to_neighbors = {}

    # build a set of cubes and record which bases are in which cube
    for t in range(len(x)):
        i = math.floor(x[t]/distance)
        j = math.floor(y[t]/distance)
        k = math.floor(z[t]/distance)
        key = "%d,%d,%d" % (i,j,k)
        if key in key_to_points:
            key_to_points[key].append((x[t],y[t],z[t]))
        else:
            key_to_points[key] = [((x[t],y[t],z[t]))]
            # map this new cube to all of its neighbors, starting with itself
            key_to_neighbors[key] = []
            for a in [0,-1,1]:
                for b in [0,-1,1]:
                    for c in [0,-1,1]:
                        neighbor_key = "%d,%d,%d" % (i+a,j+b,k+c)
                        key_to_neighbors[key].append(neighbor_key)

    # map cubes without points back to their neighbors
    # loop over cubes with points
    for key in key_to_points.keys():
        # loop over neighbors of that cube
        for neighbor_key in key_to_neighbors[key]:
            # if the neighbor has no points,
            if not neighbor_key in key_to_points:
                # if the neighbor cube has already been seen, map back to this neighbor
                if neighbor_key in key_to_neighbors:
                    key_to_neighbors[neighbor_key].append(key)
                else:
                    # otherwise, this is the first key to map back to
                    key_to_neighbors[neighbor_key] = [key]

    return key_to_points, key_to_neighbors


def filter_neighboring_residues_with_cubes(key_to_points, key_to_neighbors, potential_neighboring_units, distance):
    # for each potential neighboring unit, find the cube it is in,
    # then loop over that cube and all neighbors to check distances to query units
    distance_squared = distance * distance
    neighboring_residues = []
    chains_found = set()
    for unit in potential_neighboring_units:
        fields = unit['unit_id'].split("|")
        chain = fields[2]
        # once a neighbor in a chain is found, do not check any more units from that chain
        if not chain in chains_found:
            x = unit['x']
            y = unit['y']
            z = unit['z']
            # x = float(unit['x'])
            # y = float(unit['y'])
            # z = float(unit['z'])
            i = math.floor(x/distance)
            j = math.floor(y/distance)
            k = math.floor(z/distance)
            key = "%d,%d,%d" % (i,j,k)
            match_found = False
            # loop over this cube and up to 26 other neighboring cubes
            for neighbor_key in key_to_neighbors.get(key,[]):
                # loop over query nucleotides from the cube
                for point in key_to_points.get(neighbor_key,[]):
                    d2 = math.pow((x - point[0]), 2)  + math.pow((y - point[1]), 2) + math.pow((z - point[2]), 2)
                    if d2 < distance_squared:
                        # unit is near a query nucleotide, note its unit id and chain
                        neighboring_residues.append(unit['unit_id'])
                        chains_found.add(chain)
                        match_found = True
                        # stopping checking points in this cube
                        break
                if match_found:
                    # stopping checking neighboring cubes
                    break

    return neighboring_residues


def filter_neighboring_residues_with_cubes_tuple(key_to_points, key_to_neighbors, potential_neighboring_units, distance):
    # for each potential neighboring unit, find the cube it is in,
    # then loop over that cube and all neighbors to check distances to query units
    distance_squared = distance * distance
    neighboring_residues = []
    chains_found = set()
    for unit in potential_neighboring_units:
        id = unit[0]
        fields = id.split("|")
        chain = fields[2]
        # once a neighbor in a chain is found, do not check any more units from that chain
        if not chain in chains_found:
            x = unit[1]
            y = unit[2]
            z = unit[3]
            # x = float(unit['x'])
            # y = float(unit['y'])
            # z = float(unit['z'])
            i = math.floor(x/distance)
            j = math.floor(y/distance)
            k = math.floor(z/distance)
            key = "%d,%d,%d" % (i,j,k)
            match_found = False
            # loop over this cube and up to 26 other neighboring cubes
            for neighbor_key in key_to_neighbors.get(key,[]):
                # loop over query nucleotides from the cube
                for point in key_to_points.get(neighbor_key,[]):
                # for a,b,c in key_to_points.get(neighbor_key,[]):
                    d2 = math.pow((x - point[0]), 2) + math.pow((y - point[1]), 2) + math.pow((z - point[2]), 2)
                    # d2 = (x - point[0])**2 + (y - point[1])**2 + (z - point[2])**2
                    # d2 = (x-a)*(x-a) + (y-b)*(y-b) + (z-c)*(z-c)
                    # d2 = (x-a)**2 + (y-b)**2 + (z-c)**2
                    if d2 < distance_squared:
                        # unit is near a query nucleotide, note its unit id and chain
                        neighboring_residues.append(id)
                        chains_found.add(chain)
                        match_found = True
                        # stopping checking points in this cube
                        break
                if match_found:
                    # stopping checking neighboring cubes
                    break

    return neighboring_residues


def get_chains(unit_ids, distance=10.0):
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

    # identify the chain(s) in the query nucleotides
    # some day, we will have input units from more than one chain. This code is ready.
    chains = set()
    for unit_id in unit_ids:
        fields = unit_id.split("|")
        chains.add(fields[2])

    # Get all centers of query unit_ids, including base, sugar, phosphate
    # Returns one list of x coordinates, one list of y coordinates, one list of z coordinates
    centers_xyz_coord = get_xyz_coordinates(unit_ids, pdb_id)

    # Check if any of the coordinate lists is empty and return an empty list or handle it accordingly
    # Very unlikely to happen
    if any(len(coords) == 0 for coords in centers_xyz_coord):
        print("No coordinates found for the given unit IDs {}".format(unit_ids))
        return []

    # Find the maxima and minima of query coordinates and expand by distance
    x_min = min(centers_xyz_coord[0]) - distance
    x_max = max(centers_xyz_coord[0]) + distance
    y_min = min(centers_xyz_coord[1]) - distance
    y_max = max(centers_xyz_coord[1]) + distance
    z_min = min(centers_xyz_coord[2]) - distance
    z_max = max(centers_xyz_coord[2]) + distance

    # store the limits in an array
    coord_limits = [x_min, x_max, y_min, y_max, z_min, z_max]

    # query to find all units whose x, y, z coordinates are between the limits
    # that is all units in a rectangular box, which is quite large
    # units from the same chain(s) as the query are excluded
    # potential_neighboring_units = get_xyz_coordinates_between_limits(pdb_id, model_num, coord_limits)
    # potential_neighboring_units = get_xyz_coordinates_between_limits_different_chain(pdb_id, model_num, coord_limits, chains)

    potential_neighboring_units = get_xyz_coordinates_between_limits_different_chain_tuple(pdb_id, model_num, coord_limits, chains)

    if len(potential_neighboring_units) == 0:
        return []

    # make it easier to test two methods of finding the actual nearby chains
    method = [2]
    neighboring_chains = []

    if len(method) > 1:
        logging.info(" ")

    if 1 in method:
        # check distances of potential to centers_xyz_coord
        # calculate distances to units in unit_ids
        # neighboring_residues = filter_neighboring_residues_unique_chain(centers_xyz_coord, potential_neighboring_units, distance)
        neighboring_residues = filter_neighboring_residues_unique_chain_tuple(centers_xyz_coord, potential_neighboring_units, distance)

        # simplify the list of chains
        neighboring_chains = sorted(get_possible_chains_list(neighboring_residues, chain))
        if len(method) > 1:
            logging.info("Neighboring chains method 1: %s", neighboring_chains)

    if 2 in method:
        # make xyz cubes for the query units
        key_to_points, key_to_neighbors = make_xyz_cubes(centers_xyz_coord[0], centers_xyz_coord[1], centers_xyz_coord[2], distance)

        # find the neighboring units that are within distance of the query units
        # neighboring_residues = filter_neighboring_residues_with_cubes(key_to_points, key_to_neighbors, potential_neighboring_units, distance)
        neighboring_residues = filter_neighboring_residues_with_cubes_tuple(key_to_points, key_to_neighbors, potential_neighboring_units, distance)
        neighboring_chains = sorted(get_possible_chains_list(neighboring_residues, chain))
        if len(method) > 1:
            logging.info("Neighboring chains method 2: %s", neighboring_chains)


    # query the database to get the chain names
    neighboring_chains_with_names = get_chain_name_revised(neighboring_chains)

    if contains_mrna_substring(neighboring_chains_with_names):
        formatted_names = replace_with_mRNA(neighboring_chains_with_names)
    else:
        formatted_names = neighboring_chains_with_names

    return sorted(formatted_names, key=lambda x: x[1])

def get_neighbors(units_list, max_workers=4):
    """
    This function runs `get_chains` in parallel using ThreadPoolExecutor.
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        result = list(executor.map(get_chains, units_list))
    return result