import argparse
from GROSeqPL_v3_updated import extract_rpkm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GRO-Seq Pipeline')
    parser.add_argument("-a", dest = "annotation", type = str, required = False, help = "Genome annotation (.gtf/.gtf.tz file)")
    parser.add_argument("-o", dest="output", type=str, required=True, help="Output name")
    args = parser.parse_args()

    print("RPKM extracted")

    extract_rpkm(args)