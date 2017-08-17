import os
CONF_ROOT = os.path.dirname(os.path.abspath(__file__))
os.environ['CONF_ROOT'] = CONF_ROOT

# from random import randint
# randomizer = (randint(0, 10))
# if randomizer == 8:
# 	import requests
# 	import json
# 	# Error types
# 	class authorisationException(Exception):
# 		pass
# 	# Check auth status
# 	headers = {'Content-Type': 'application/json'}
# 	USER = os.environ.get('USERNAME')
# 	PASS = os.environ.get('PASSWORD')
# 	try:
# 		authorisation = requests.get('https://rest.variantvalidator.org/api/myaccount', auth=(USER, PASS))
# 	except:
# 		raise authorisationException('Unable to reach Auth server')
# 	else:
# 		auth = (authorisation.text)
# 		try:
# 			auth_data = json.loads(auth)
# 		except:
# 			raise authorisationException('Authorisation required')
# 		if auth_data['status'] == 'LIVE':
# 			pass
# 		else:
# 			raise authorisationException('Account expired on %s') %(auth_data['expiry'])
					