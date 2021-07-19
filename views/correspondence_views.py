import flask
from flask import render_template, request

blueprint = flask.Blueprint('correspondence_list', __name__, template_folder='templates')


@blueprint.route('/corr-server/list/<query>')

	return query
