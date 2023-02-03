#!/usr/bin/env python3

# author: Ole Schuett

from typing import Any, Dict, Set, List, Optional, NewType, cast
from urllib.request import urlopen
from datetime import datetime
from glob import glob
import configparser
import dataclasses
import numpy as np
import gzip
import sys
import re

# ======================================================================================
Report = NewType("Report", Dict[str, str])

# ======================================================================================
@dataclasses.dataclass
class TestDef:
    test_type: int
    flags: List[str]
    tolerance: float
    ref_value: str


# ======================================================================================
def main() -> None:
    if len(sys.argv) != 3:
        print("Usage generate_regtest_survey.py <dashboard.conf> <output-dir>")
        sys.exit(1)

    config_fn, outdir = sys.argv[1:]
    assert outdir.endswith("/")

    # parse ../../tests/*/*/TEST_FILES
    test_defs = parse_test_files()

    # parse ../../tests/TEST_TYPES
    test_types = parse_test_types()

    # find eligible testers by parsing latest reports
    tester_names: List[str] = list()
    tester_values: Dict[str, Report] = dict()
    inp_names: Set[str] = set()
    config = configparser.ConfigParser()
    config.read(config_fn)

    def get_sortkey(s: str) -> int:
        return config.getint(s, "sortkey")

    for tname in sorted(config.sections(), key=get_sortkey):
        fn = outdir + "archive/%s/list_recent.txt" % tname
        list_recent = open(fn, encoding="utf8").readlines()
        if not list_recent:
            continue  # list_recent is empty
        latest_report_url = list_recent[0].strip()
        report = parse_report(latest_report_url)
        if report:
            inp_names.update(report.keys())
            tester_values[tname] = report
            tester_names.append(tname)

    # remove outdated inp-names
    inp_names = inp_names.intersection(test_defs.keys())

    # html-header
    output = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
    output += "<html><head>\n"
    output += '<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n'
    output += '<script type="text/javascript" src="https://code.jquery.com/jquery-2.1.4.min.js"></script>\n'
    output += '<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jquery.tablesorter/2.23.2/js/jquery.tablesorter.min.js"></script>\n'
    output += '<script type="text/javascript">\n'
    output += "$(document).ready(function(){\n"
    output += '    $("table").tablesorter();\n'
    output += '    $("table").bind("sortStart",function(){ $(".waitmsg").show(); });\n'
    output += '    $("table").bind("sortEnd",  function(){ $(".waitmsg").hide(); });\n'
    output += "  }\n"
    output += ");\n"
    output += "</script>\n"
    output += '<style type="text/css">\n'
    output += ".nowrap { white-space: nowrap; }\n"
    output += "tr:hover { background-color: #ffff99; }\n"
    output += ".waitmsg {\n"
    output += "  color: red;\n"
    output += "  font: bold 20pt sans-serif;\n"
    output += "  display: none;\n"
    output += "}\n"
    output += "#factbox {\n"
    output += "  display: inline-block;\n"
    output += "  border-radius: 1em;\n"
    output += "  box-shadow: .2em .2em .7em 0 #777;\n"
    output += "  background: #f7f7f0;\n"
    output += "  padding: 1em;\n"
    output += "  margin: 20px;\n"
    output += "}\n"
    output += "#factbox h2 { margin: 0 0 0.5em 0; }\n"
    output += "</style>\n"
    output += "<title>CP2K Regtest Survey</title>\n"
    output += "</head><body>\n"
    output += "<center><h1>CP2K REGTEST SURVEY</h1></center>\n"

    # fun facts
    output += '<div id="factbox"><table>\n'
    ntests = len(test_defs)
    output += "<tr><td>Total number of test-cases</td>"
    output += '<td align="right">%d</td><td align="right">100.0%%</tr>\n' % ntests
    n = sum(len(t.flags) != 0 for t in test_defs.values() if t.flags)
    pc = n / (0.01 * ntests)
    output += "<tr><td>Tests which require flags</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (n, pc)
    n = sum(t.test_type != 0 for t in test_defs.values())
    pc = n / (0.01 * ntests)
    output += "<tr><td>Numeric tests, ie. type &ne; 0</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (n, pc)
    n = sum(t.ref_value != "" for t in test_defs.values())
    pc = n / (0.01 * ntests)
    output += "<tr><td>Numeric tests with fixed reference</td>"
    output += '<td align="right">%d</td><td align="right">%.1f%%</td></tr>\n' % (n, pc)
    for i in range(14, 5, -1):
        tol = float("1.e-%d" % i)
        n = sum(t.test_type != 0 and t.tolerance <= tol for t in test_defs.values())
        pc = n / (0.01 * ntests)
        output += "<tr><td>Numeric tests with tolerance &le; 10<sup>-%d</sup></td>" % i
        output += '<td align="right">%d</td><td align="right">%.1f%%</td>' % (n, pc)
        output += "</tr>\n"
    output += "</table></div>\n"

    # table-body
    tester_nskipped = dict([(n, 0) for n in tester_names])
    tester_nfailed = dict([(n, 0) for n in tester_names])
    tbody = ""
    for inp in sorted(inp_names):
        # calculate median and MAD
        value_list = list()
        for tname in tester_names:
            val = tester_values[tname].get(inp, None)
            if val:
                value_list.append(float(val))
        values: Any = np.array(value_list)
        median_iterp = np.median(values)  # takes midpoint if len(values)%2==0
        median_idx = (np.abs(values - median_iterp)).argmin()  # find closest point
        median = values[median_idx]
        norm = median + 1.0 * (median == 0.0)
        rel_diff = abs((values - median) / norm)
        mad = np.amax(rel_diff)  # Maximum Absolute Deviation
        outlier = list(rel_diff > test_defs[inp].tolerance)

        # output table-row
        tbody += '<tr align="right">\n'
        style = 'bgcolor="#FF9900"' if any(outlier) else ""
        tbody += '<th align="left" %s>%s</th>\n' % (style, inp)
        ttype_num = test_defs[inp].test_type
        ttype_re = test_types[test_defs[inp].test_type].split("!")[0].strip()
        tbody += '<td title="%s" >%s</td>\n' % (ttype_re, ttype_num)
        tbody += "<td>%.1e</td>\n" % test_defs[inp].tolerance
        tbody += "<td>%.1e</td>\n" % mad
        tol_mad = test_defs[inp].tolerance / max(1e-14, mad)
        if tol_mad < 1.0:
            tbody += '<td bgcolor="#FF9900">%.1e</td>\n' % tol_mad
        elif tol_mad > 10.0:
            tbody += '<td bgcolor="#99FF00">%.1e</td>\n' % tol_mad
        else:
            tbody += "<td>%.1e</td>\n" % tol_mad
        tbody += "<td>%s</td>\n" % test_defs[inp].ref_value
        tbody += "<td>%.17g</td>\n" % median
        tbody += "<td>%i</td>\n" % np.sum(outlier)
        for tname in tester_names:
            val = tester_values[tname].get(inp, None)
            if not val:
                tbody += "<td></td>\n"
                tester_nskipped[tname] += 1
            elif outlier.pop(0):
                tbody += '<td bgcolor="#FF9900">%s</td>' % val
                tester_nfailed[tname] += 1
            else:
                tbody += "<td>%s</td>" % val
        tbody += '<th align="left" %s>%s</th>\n' % (style, inp)
        tbody += "</tr>\n"

    # table-header
    theader = '<tr align="center"><th>Name</th><th>Type</th><th>Tolerance</th>'
    theader += '<th><abbr title="Maximum Absolute Deviation">MAD</abbr></th>'
    theader += "<th>Tol. / MAD</th><th>Reference</th><th>Median</th>"
    theader += '<th><abbr title="Failures if ref. were at median.">#failed</abbr></th>'
    for tname in tester_names:
        theader += '<th><span class="nowrap">%s</span>' % config.get(tname, "name")
        theader += "<br>#failed: %d" % tester_nfailed[tname]
        theader += "<br>#skipped: %d" % tester_nskipped[tname]
        theader += "</th>\n"
    theader += "<th>Name</th>"
    theader += "</tr>\n"

    # assemble table
    output += '<div class="waitmsg">Sorting, please wait...</div>\n'
    output += "<p>Click on table header to sort by column.</p>\n"
    output += '<table border="1" cellpadding="5">\n'
    output += "<thead>" + theader + "</thead>"
    output += "<tfoot>" + theader + "</tfoot>"
    output += "<tbody>" + tbody + "</tbody>"
    output += "</table>\n"

    # html-footer
    now = datetime.utcnow().replace(microsecond=0)
    output += "<p><small>Page last updated: %s</small></p>\n" % now.isoformat()
    output += "</body></html>"

    # write output file
    fn = outdir + "regtest_survey.html"
    f = open(fn, "w", encoding="utf8")
    f.write(output)
    f.close()
    print("Wrote: " + fn)


# ======================================================================================
def parse_test_files() -> Dict[str, TestDef]:
    test_defs = dict()

    tests_root = "../../tests/"
    test_dir_lines = open(tests_root + "TEST_DIRS", encoding="utf8").readlines()
    for dline in test_dir_lines:
        dline = dline.strip()
        if len(dline) == 0:
            continue
        if dline.startswith("#"):
            continue
        d = dline.split()[0]
        flags = dline.split()[1:]  # flags requiremented by this test_dir
        fn = tests_root + d + "/TEST_FILES"
        content = open(fn, encoding="utf8").read()
        for line in content.strip().split("\n"):
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            name = d + "/" + parts[0]
            test_type = int(parts[1])
            if len(parts) == 2:
                tolerance = 1.0e-14  # default
                ref_value = ""
            elif len(parts) == 3:
                tolerance = float(parts[2])
                ref_value = ""
            elif len(parts) == 4:
                tolerance = float(parts[2])
                ref_value = parts[3]  # do not parse float
            else:
                raise (Exception("Found strange line in: " + fn))
            test_defs[name] = TestDef(test_type, flags, tolerance, ref_value)

    return test_defs


# ======================================================================================
def parse_test_types() -> List[str]:
    test_types = [""]
    lines = open("../../tests/TEST_TYPES", encoding="utf8").readlines()
    ntypes = int(lines[0])
    for i in range(1, ntypes + 1):
        test_types.append(lines[i])
    return test_types


# ======================================================================================
def parse_report(report_url: str) -> Optional[Report]:
    print("Parsing: " + report_url)
    data = urlopen(report_url, timeout=5).read()
    report_txt = gzip.decompress(data).decode("utf-8", errors="replace")

    if "Keepalive:" not in report_txt:
        print("Could not recognize as regtests report - skipping.")
        return None

    base_dir = None
    curr_dir = None
    values: Report = cast(Report, {})
    for line in report_txt.split("\n"):
        if line.startswith("Work base dir:"):
            base_dir = line.split(":", 1)[1].strip()
        elif "/UNIT/" in line:
            continue  # ignore unit-tests
        elif line.startswith(">>>"):
            prefix = f">>> {base_dir}/"
            assert line.startswith(prefix)
            curr_dir = line[len(prefix) :]
        elif line.startswith("<<<"):
            curr_dir = None
        elif curr_dir:
            # Found an actual result line.
            if "OK" not in line and "WRONG RESULT" not in line:
                continue  # ignore failed tests
            parts = line.split()
            if parts[0].rsplit(".", 1)[1] not in ("inp", "restart"):
                print("Found strange line:\n" + line)
                continue
            test_name = curr_dir + "/" + parts[0]
            if parts[1] == "-":
                continue  # test without numeric check
            try:
                float(parts[1])  # try parsing float...
                values[test_name] = parts[1]  # ... but pass on the original string
            except ValueError:
                pass  # ignore values which can not be parsed

    return values


# ======================================================================================
main()

# EOF
