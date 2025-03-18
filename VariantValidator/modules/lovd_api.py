import requests

def lovd_syntax_check(variant_description):
    base_url = "https://api.lovd.nl/v2/checkHGVS"
    url = f"{base_url}/{variant_description}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        json_data = response.json()
        json_data["url"] = url
        return json_data
    except requests.RequestException as e:
        return {"lovd_api_error": str(e)}

# <LICENSE>
# Copyright (C) 2016-2025 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>