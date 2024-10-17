# Revised for the new server in 2024

import sys
import logging
import traceback

# these lines help on the new server in 2024
sys.path.insert(0, '/var/www/correspondence')
sys.path.insert(1, '/var/www/correspondence/app')
sys.path.insert(2, '/var/www/correspondence/app/fr3d')
sys.path.insert(3, '/var/www/correspondence/app/fr3d/geometry')

# removed the encoding='utf-8' for the new server in 2024
# logging for regular users
logging.basicConfig(filename='/var/www/correspondence/flask.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

# logging when trying to run as the pipeline user just for testing
# logging.basicConfig(filename='/usr/local/pipeline/flask.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

logging.info("==============================")
logging.info("Starting app import in wsgi.py")

from app import app as application  # Import your Flask app

logging.info("Starting wsgi run")

# these are apparently optional
application.config['ENV'] = 'production'
application.config['DEBUG'] = False

# try:
#     application.run()
#     # application.run(port=5000)
#     logging.info("Started wsgi run, now what?")
# except Exception as e:
#     logging.error("Error in wsgi.py: %s" % e)
#     logging.error(traceback.format_exc())

# logging.info("Got past the try/except in wsgi.py")

# logging.basicConfig(stream=sys.stderr,level=logging.DEBUG)

logging.info('End of wsgi.py')
