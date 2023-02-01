#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
import ast
from os import path
from argparse import RawTextHelpFormatter  # enable newline

blas_re = re.compile(
    r"[SDCZ]"
    r"(ROTG|ROTMG|ROT|ROTM|SWAP|SCAL|COPY|AXPY|DOT|DOTU|DOTC"  # level 1
    r"|GEMV|GBMV|HEMV|HBMV|HPMV|SYMV|SBMV|SPMV|TRMV|TBMV|TPMV|TRSV"  # level 2
    r"|TBSV|TPSV|GER|GERU|GERC|HER|HPR|HER2|HPR2|SYR|SPR|SYR2|SPR2"
    r"|GEMM|SYMM|HEMM|SYRK|HERK|SYR2K|HER2K|TRMM|TRSM"  # level 3
    r"|LANGE|LARNV|LAMCH|NRM2|CNRM2|ZNRM2|DSCAL)"
)  # aux

lapack_re = re.compile(
    r"ILAENV|"
    r"([SDCZ]" + "(BD|DI|GB|GE|GG|GT|HB|HE|HG|HP|HS|OP"
    r"|OR|PB|PO|PP|PT|SB|SP|ST|SY|TB|TG|TP|TR|TZ|UN|UP)"
    r"(BAK|BAL|BRD|CON|EBZ|EDC|EIN|EQR|EGR|EQU|EQZ|ERF|EVC"
    r"|EXC|GBR|GHR|GLQ|GQL|GQR|GRQ|GST|GTR|HRD|LQF|MBR|MHR"
    r"|MLQ|MQL|MQR|MRQ|MRZ|MTR|QLF|QPF|QRF|RFS|RQF|RZF|SDC"
    r"|SEN|SJA|SNA|SQR|SVP|SYL|TRD|TRF|TRI|TRS"
    r"|SDD|EV|GV|SV|BS2D|BR2D|LS))"
)

warning_re = re.compile(r".*[Ww]arning: (.*)")

IGNORED_WARNINGS = (
    "-Wrealloc-lhs",
    "-Wdo-subscript",
    "-Winteger-division",
    "-Wmaybe-uninitialized",
    "-Wfunction-elimination",
    "-Wunused-function",
    "Creating array temporary",
    "quality comparison",
    "Removing call to function",
    "Conversion from",
    "Non-significant digits in 'REAL(8)'",
    "POINTER-valued function appears on right-hand side",
    "style of line directive is a GCC extension",
)

# ======================================================================================
def check_warnings(fhandle):
    loc = loc_short = ""

    for line in fhandle:
        line = line.rstrip("\n")

        if not line:
            continue

        # line directives issues by gcc directly
        # we ignore non-absolut paths, because they point to indermediate fypp output
        if line.startswith("/") and line.endswith(":"):
            loc = line.rsplit(":")[0].strip()
            loc_short = path.basename(loc)
            if not path.exists(loc):
                return  # source file gone - skipping
            continue

        # fypp line directives that leaked through as part of warning messages
        if line.startswith(' # 1 "'):
            loc = line.split()[2].strip('"')
            loc_short = path.basename(loc)
            if not path.exists(loc):
                return  # source file gone - skipping
            continue

        if loc.endswith("include/fftw3.f"):
            continue  # an external file

        m = warning_re.match(line)
        if not m:
            continue  # we are only looking for warnings
        warning = m.group(1)

        warning = re.sub(r"[‘’]", "'", warning)  # replace unicode apostrophes
        warning = re.sub(r"'\d+'", "*", warning)  # replace numbers with *

        if any(iw in warning for iw in IGNORED_WARNINGS):
            continue

        if ("CHARACTER expression" in warning) and ("truncated" in warning):
            continue

        if "called with an implicit interface" in warning:
            parts = warning.split()
            assert parts[0] == "Procedure"
            routine = parts[1].strip("'").upper()
            if may_call_implicit(loc, routine):
                continue

        if ".ubound' is used uninitialized" in warning:
            continue
        if ".stride' is used uninitialized" in warning:
            continue

        # ok this warning we should handle
        print("%s: %s" % (loc_short, warning))


# ======================================================================================
def may_call_implicit(loc, routine):
    if blas_re.match(routine):
        return True  # BLAS calls are allowed everywhere

    if lapack_re.match(routine):
        return True  # Lapack calls are allowed everywhere

    pkg = path.dirname(loc)
    manifest_fn = pkg + "/PACKAGE"
    with open(manifest_fn, encoding="utf8") as fhandle:
        manifest = ast.literal_eval(fhandle.read())

    if "implicit" not in manifest:
        return False

    return re.match(manifest["implicit"], routine)


# ======================================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Checks given stderr output files from gfortran for violations of the coding conventions",
        epilog="""\
For generating the warn-files run gfortran with all warning flags and redirect the output to a file.
This can be achieved by putting
    FCLOGPIPE   = 2>&1 | tee $(notdir $<).warn
in the cp2k arch-file.
""",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "files",
        metavar="<warn-file>",
        type=str,
        nargs="+",
        help="files containing the compiler warnings",
    )
    args = parser.parse_args()

    for fn in args.files:
        with open(fn, encoding="utf8") as fhandle:
            check_warnings(fhandle)

# EOF
