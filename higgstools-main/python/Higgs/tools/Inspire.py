import requests
import re
from collections import OrderedDict
from typing import Union


def apiUrl(id: Union[str, int]):
    """
    Convert the id into a valid inspirehep API url.

    Parameters:
    -----------
      id: str or int
        The url of an inspire entry, a numeric inspire id, an
        inspire cite key, or an arxiv number.
    """
    from urllib.parse import urlparse, urljoin

    id = str(id)
    search = "https://inspirehep.net/api/literature/?q="
    if ":" in id:  # the id is a citekey
        search = search + "texkeys:" + id
    if "." in id or id.startswith("hep-"):  # the id is an arxiv id
        search = search + "e " + id
    else:
        path = urlparse(id).path
        if not path.startswith("/"):
            path = "/literature/" + path
        if not path.startswith("/api"):
            path = "/api" + path
        return urljoin("https://inspirehep.net/", path)
    try:
        return requests.get(search).json()["hits"]["hits"][0]["links"]["json"]
    except IndexError:
        raise RuntimeError(
            "Could not find inspire entry based on arXiv id. Give the numeric inspire id instead."
        )


def experiment(j: dict, collider: str):
    """
    Obtain the experiment from the inspire entry.
    """
    collabs = []
    if "collaborations" in j.keys():
        collabs = [x["value"] for x in j["collaborations"]]
        if len(collabs) > 1:
            if collider.startswith("LHC"):
                return "LHCComb"
            else:
                return collider + "Comb"
    if len(collabs) == 0:
        print("Could not detect experiment")
        return ""
    return collabs[0]


def luminosity(j: dict):
    """
    Parse the luminosity from the abstract.
    """
    if "abstracts" in j:
        abstracts = " ".join([x["value"] for x in j["abstracts"]])
        lumiRE = re.compile(
            r"([0-9\.]+)[\s~$]*(?:(?:fb|\\mathrm{fb})[$^{\s]*[-−]1|inverse femtobarn)"
        )
        lumis = [float(x) for x in set(lumiRE.findall(abstracts))]
        if len(lumis) > 1:
            print(
                "Multiple different luminosity values:",
                lumis,
                "using",
                lumis[0],
                "unless you set it manually",
            )
        if len(lumis) > 0:
            return lumis[0]
    return -1.0


def reference(j: dict, experiment: str):
    """
    Get the arxiv id from inspire. The experiment is used to select the
    appropriate report number, if no arxiv number is available.
    """
    try:
        return j["arxiv_eprints"][0]["value"]
    except:
        print("No arxiv eprint found.")
    try:
        return [x["value"] for x in j["report_numbers"] if experiment in x["value"]][0]
    except:
        if "report_numbers" in j.keys():
            print(
                "none of the report numbers",
                [x["value"] for x in j["report_numbers"]],
                "matches the detected experiment",
                experiment,
            )
        else:
            print("no report number found")
    return ""


def citeKey(j: dict):
    return j["texkeys"][0]


def collider(j: dict):
    """
    Parse the collider (and energy for the LHC) from the abstract.
    """
    if "abstracts" in j:
        abstracts = " ".join([x["value"] for x in j["abstracts"]])
        if re.search(r"\bLEP\b", abstracts):
            return "LEP"
        else:
            LHCE = re.compile(r"([0-9]+)[\s~$]*(?:\\mathrm{)?TeV")
            energies = set(LHCE.findall(abstracts))
            if "13" in energies:
                return "LHC13"
            elif "8" in energies:
                return "LHC8"
            elif "7" in energies:
                return "LHC7"
    return ""


def getMetadata(inspireLink: Union[str, int]):
    """
    Use the inspire data to obtain limit metadata. This will return dictionary
    that contains all the top-level entries required for a HiggsBounds limit. If
    possible, their values are filled using the information available on
    InspireHEP.

    inspireLink: str or int
      The url of an inspire entry, a numeric inspire id, an
      inspire cite key, or an arxiv number.

    Entries:
    --------
      "limitClass": must be appropriately filled by the user
      "reference": The arxiv id, if it exists, a report number that matches the experiment, otherwise.
      "id": A default value is generated from the arxiv id. If a
      single paper contains multiple limits the ids of those limits have to be
      made distinct (e.g. by adding another digit to enumerate them).
      "citeKey": The inspireHep cite key.
      "source": Should be manually filled with the concrete data source (table or figure number).
      "collider": The collider id. This information is parsed from the abstract and may be unreliable.
      "experiment": The experimental collaboration.
      "luminosity": The luminosity used for this limit. This information is parsed from the abstract and may be unreliable.
      "process": Just an empty dict, needs to be filled appropriately with the process definition.
      "analysis": Just an empty dict, needs to be filled appropriately with the date.
    """
    response = requests.get(apiUrl(inspireLink))
    j = response.json()["metadata"]

    data = OrderedDict()
    coll = collider(j)
    exp = experiment(j, coll)
    ref = reference(j, exp)
    data["limitClass"] = ""
    if ref:
        try:
            data["id"] = int(ref.replace(".", ""))
        except:
            pass
        try:
            data["id"] = int("".join(filter(lambda x: x.isdigit(), ref)))
        except:
            data["id"] = 0
    else:
        data["id"] = j["control_number"]
    data["reference"] = ref
    data["source"] = ""
    data["citeKey"] = citeKey(j)
    data["collider"] = coll
    data["experiment"] = exp
    data["luminosity"] = luminosity(j)
    return data
