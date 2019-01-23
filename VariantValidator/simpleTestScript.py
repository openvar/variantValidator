import json
import os

from VariantValidator import Validator

#variant = 'NM_000088.3:c.589G>T'
variant = 'NC_000012.11:g.122064776delG'
select_transcripts = 'all'
selected_assembly = 'GRCh37'

validator=Validator()
out=Validator().validate(variant, selected_assembly, select_transcripts)

print json.dumps(out, sort_keys=True, indent=4, separators=(',', ': '))
