"""
This file will contain the logic for determining the type of input
Generally, input can be divided into three types: string of unit_ids,
ranges and loop_id
"""

import logging
import collections as coll

logging.info("Importing process_input")

from werkzeug.exceptions import BadRequest

PARTS_LIMIT = 40


def input_type(data):
    parts = data.split(',')
    if len(parts) == 0 or len(parts) > PARTS_LIMIT:
        raise BadRequest("Must give 1 to %s parts in a collection" %
                         PARTS_LIMIT)

    processed = []
    for part in parts:
        if not part:
            raise BadRequest("Cannot give empty part of a list")
        units = part.split(':')
        if len(units) == 1:
            if not units[0]:
                raise BadRequest("Cannot give empty unit")
            processed.append(tuple([units[0], units[0]]))
        elif len(units) == 2:
            if not units[0] or not units[1]:
                raise BadRequest("Both units in range must exist")
            processed.append(tuple(units))
        else:
            raise BadRequest("Range should must have 1 or 2 endpoints")

    return processed


def check_query(query):
    """
    This check the type of query. It can be of three types:
    i)   loop_id
    ii)  single range (for a HL)
    iii) multiples ranges (for a IL)
    iv)  string of unit ids
    """

    if len(query) == 1:
        if query[0][0] == query[0][1]:
            query_type = 'loop_id'
        else:
            query_type = 'single_range'
    else:
        if query[0][0] != query[0][1]:
            query_type = 'multiple_ranges'
        else:
            query_type = 'units_str'

    return query_type


def check_input_type(loop_id, unit_id, res_num):

    if loop_id is not None and unit_id is None and res_num is None:
        return "loop_id"
    elif unit_id is not None and loop_id is None and res_num is None:
        return "unit_id"
    elif res_num is not None and loop_id is None and unit_id is None:
        return "res_num"
