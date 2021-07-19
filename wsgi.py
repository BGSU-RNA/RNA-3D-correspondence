import os
import sys
import logging
#from logging import FileHandler, WARNING


sys.path.append(os.path.dirname(__file__))


from app import app as application

logging.basicConfig(stream=sys.stderr,level=logging.DEBUG)

