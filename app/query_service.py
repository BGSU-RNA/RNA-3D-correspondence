import itertools
import logging
import re

logging.info("Importing query_service 1")

# from sqlalchemy import or_

from models import UnitInfo, LoopPositions
logging.info("Importing query_service 2")
from database import db_session
logging.info("Importing query_service 3")
import utility as ui

logging.info("Importing query_service 4")

LOOP_REGEX_PATTERN = r'^(IL|HL|J3)_[0-9A-Z]{4}_\d{3}$'

def get_units(incomplete_units):

    with db_session() as session:

        complete_units = []
        for unit in incomplete_units:
            query = session.query(UnitInfo).filter(UnitInfo.unit_id.like(unit))
            for row in query:
                complete_units.append(row.unit_id)

        return complete_units

def get_units_new(incomplete_unit):

    with db_session() as session:
        query = session.query(UnitInfo).filter(UnitInfo.unit_id.like(incomplete_unit))

        output_list = [result for result in query]

        if len(output_list) == 0:
            return None
        else:
            return query[0].unit_id

def get_single_range_units(range_positions, pdb_id, chain):

    # Previously was using sequence position, known as chain_index.
    # .filter(UnitInfo.chain_index.between(range_positions[0], range_positions[1])) \


    with db_session() as session:

        complete_units = []
        query = session.query(UnitInfo).filter_by(pdb_id=pdb_id) \
                                       .filter_by(chain=chain) \
                                       .filter(UnitInfo.number.between(range_positions[0], range_positions[1])) \
                                       .order_by(UnitInfo.chain_index)

        for row in query:
            complete_units.append(row.unit_id)

        return complete_units


def get_multiple_range_units(range_positions_list, pdb_id, chain):

    with db_session() as session:

        complete_units = []
        for single_range in range_positions_list:
            single_range_units = []
            query = session.query(UnitInfo).filter_by(pdb_id=pdb_id) \
                                       .filter_by(chain=chain) \
                                       .filter(UnitInfo.chain_index.between(single_range[0], single_range[1])) \
                                       .order_by(UnitInfo.chain_index)

            for row in query:
                single_range_units.append(row.unit_id)

            complete_units.append(single_range_units)

        #merged_units = list(itertools.chain(*complete_units))
        # return the flattened list of lists
        return list(itertools.chain(*complete_units))

def get_loop_units(loop_id):

    with db_session() as session:

        complete_units = []
        query = session.query(LoopPositions).filter_by(loop_id=loop_id).order_by(LoopPositions.position_2023)

        for row in query:
            complete_units.append(row.unit_id)

        return complete_units

def get_query_units(query_type, query_list, query_ife):

    if query_type == 'units_str':
        incomplete_units = [unit[0] for unit in query_list]
        ife = query_ife + '|%|'
        incomplete_units = [ife + unit for unit in incomplete_units]
        complete_units = get_units(incomplete_units)

    elif query_type == 'single_range':
        pass

    elif query_type == 'multiple_ranges':
        pass

    elif query_type == 'loop_id':
        loop_id = query_list[0][0]
        # unsorted_units = get_loop_units(loop_id)
        # complete_units = ui.get_sorted_units(unsorted_units)
        complete_units = get_loop_units(loop_id)

    return complete_units


def get_query_units_new(input_type, selection, chain_id):

    if input_type == 'loop_id':
        loop_id = selection
        complete_units = get_loop_units(loop_id)
        # complete_units = ui.get_sorted_units(unsorted_units)

    elif input_type == 'unit_id':
        # complete_units = selection.split(",")
        complete_units = re.split(',|\t', selection)

    elif input_type == 'res_num':
        incomplete_units = selection.split(",")
        chain = chain_id + '|%|'
        incomplete_units = [chain + unit for unit in incomplete_units]
        complete_units = get_units(incomplete_units)

    elif input_type == 'single_range':
        pdb_id, model, chain = chain_id.split('|')
        range_positions = selection.split(':')
        complete_units = get_single_range_units(range_positions, pdb_id, chain)

    elif input_type == 'multiple_ranges':
        pdb_id, model, chain = chain_id.split('|')
        ranges_selection = selection.split(",")
        range_positions_list = [single_range.split(':') for single_range in ranges_selection]
        complete_units = get_multiple_range_units(range_positions_list, pdb_id, chain)

    return complete_units

def is_unit_id(item):

    fields = item.split("|")

    if 4 <= len(fields) <= 9:
        if fields[3] == "HOH":
            return False
        else:
            return True
    else:
        return False

def is_loop_id(item):
    if re.match(LOOP_REGEX_PATTERN, item):
        return True
    else:
        return False

def is_residue_num(item):
    if item.isdigit():
        num = int(item)
        if 1 <= num <= 99999:
            return True
        else:
            return False
    else:
        return False

def get_query_units_modified(param_dict):

    complete_units = []
    error_message = ""
    range_separator = ":"

    selection = param_dict['selection']

    # replace space colon with colon
    selection = re.sub(r'\s+:', ':', selection)
    # replace colon space with colon
    selection = re.sub(r':\s+', ':', selection)

    # replace whitespace, tab with comma
    selection = re.sub(r'\s+', ',', selection)
    # replace multiple commas with single comma
    selection = re.sub(r',+', ',', selection)

    logging.info('Selection: ' + selection)

    items = selection.split(",")
    for item in items:

        logging.info('Item: ' + item)

        if range_separator in item:
            if (param_dict['pdb'] is None) or (param_dict['chain'] is None):
                error_message = "Need pdb and chain information when using %s operator" % range_separator
                return None, error_message
            else:
                range_positions = item.split(range_separator)
                if not len(range_positions) == 2:
                    error_message = "Need exactly two range positions"
                    return None, error_message
                for position in range_positions:
                    if not position.isdigit():
                        error_message = "Range positions must be integers"
                        return None, error_message

                complete_units.extend(get_single_range_units(range_positions, param_dict['pdb'], param_dict['chain']))

        elif is_loop_id(item):
            loop_units = get_loop_units(item)
            if not loop_units or len(loop_units) == 0:
                error_message = "Invalid loop id %s" % item
                return None, error_message
            complete_units.extend(loop_units)
        elif is_unit_id(item):
            # check to see if this is a known unit id; if not, give an error message
            complete_unit = get_units_new(item)
            if not complete_unit:
                error_message = "Invalid unit id %s" % item
                return None, error_message
            if is_unit_id(complete_unit):
                complete_units.append(item)
        elif is_residue_num(item):
            # use pdb id, chain, and residue number to get the unit id
            # assumption here is that the model num is always 1

            if (param_dict['pdb'] is None) or (param_dict['chain'] is None):
                error_message = "Need pdb and chain information when using residue number"
                return None, error_message

            incomplete_unit = param_dict['pdb'] + "|1|" + param_dict['chain'] + "|%|" + str(item)
            complete_unit = get_units_new(incomplete_unit)

            if not complete_unit:
                if param_dict['pdb'] and param_dict['chain']:
                    error_message = "Invalid residue number %s in %s chain %s" % (item,param_dict['pdb'],param_dict['chain'])
                else:
                    error_message = "Invalid residue number %s" % item
                return None, error_message

            if is_unit_id(complete_unit):
                complete_units.append(complete_unit)
            else:
                if param_dict['pdb'] and param_dict['chain']:
                    error_message = "Invalid residue number %s in %s chain %s" % (item,param_dict['pdb'],param_dict['chain'])
                else:
                    error_message = "Invalid residue number %s" % item
                return None, error_message

    return complete_units, error_message




