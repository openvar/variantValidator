from Bio import Entrez,SeqIO
import httplib2 as http
import json
from urlparse import urlparse #Python 2
import functools
#from urllib.parse import urlparse #Python 3

def handleCursor(func):
    #Decorator function for handling opening and closing cursors.
    @functools.wraps(func)
    def wrapper(self,*args,**kwargs):
        try:
            self.cursor = self.conn.cursor(buffered=True)
            out=func(*args,**kwargs)
            self.cursor.close()
            self.cursor=None
            return out
        except:
            try:
                self.cursor.close()
                self.cursor=None
            except:
                self.cursor=None
            raise
    return wrapper

def entrez_efetch(val, db, id, rettype, retmode):
    Entrez.email = val.entrezID
    handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
    record = SeqIO.read(handle, "gb")
    handle.close()
    return record

def hgnc_rest(path):
    data = {
        'record': '',
        'error': 'false'
    }
    # HGNC server
    headers = {
        'Accept': 'application/json',
    }
    uri = 'http://rest.genenames.org'
    target = urlparse(uri + path)
    method = 'GET'
    body = ''
    h = http.Http()
    # collect the response
    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)
    if response['status'] == '200':
        # assume that content is a json reply
        # parse content with the json module
        data['record'] = json.loads(content)
    else:
        data['error'] = "Unable to contact the HGNC database: Please try again later"
    return data

# method for final validation and stringifying parsed hgvs variants prior to printing/passing to html
def valstr(hgvs_variant):
    """
    Function to ensure the required number of reference bases are displayed in descriptions
    """
    cp_hgvs_variant = copy.deepcopy(hgvs_variant)
    if cp_hgvs_variant.posedit.edit.type == 'identity':
        if len(cp_hgvs_variant.posedit.edit.ref) > 1:
            cp_hgvs_variant = output_formatter.remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    else:
        cp_hgvs_variant = output_formatter.remove_reference(cp_hgvs_variant)
        cp_hgvs_variant = str(cp_hgvs_variant)
    return cp_hgvs_variant
