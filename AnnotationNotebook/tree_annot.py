import pandas as pd
import argparse
import tree_code as tc
import sys
import random
import numpy as np

def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", help="Path to tree", required=True)
    parser.add_argument("-a", "--align", help="Path to alignment", required=True)
    parser.add_argument("-c", "--csv", help="Path to csv", required=True)
    parser.add_argument("-col", "--col", help="Column to annotate based on", required=True)
    parser.add_argument("-r", "--random_seed", help="Set the random seed to an integer to make colour selection "
                                                    "reproducible")

    parser.add_argument("-o", "--outpath", help="Outpath", default="tree_annot.png")
    parser.add_argument("-m", "--match_from", help="Match from a specific column in the annotations file",
                        default="Name")

    return parser.parse_args(args)


if __name__ == "__main__":

    print ("\nRunning tree_annot")

    # Parse the arguments
    parser = parse_args(sys.argv[1:])

    df = pd.read_csv(parser.csv)


    annot_dict = dict(zip(df[parser.match_from], df[parser.col]))

    for key, val in annot_dict.items():
        if pd.isnull(val):
            annot_dict[key] = None



    # else:
    # annot_dict = df.to_dict(orient='index')

    print (annot_dict)

    col = parser.col.strip() if parser.col else df.columns[0]

    random_seed = parser.random_seed if parser.random_seed else random.randint(0, 999)


    print ("\nTo reproduce this colour scheme set the random seed as " + str(random_seed))

    color_dict = tc.get_color_dict(annot_dict, random_seed)

    # print (color_dict)

    # Load tree
    tree = tc.load_tree(parser.tree, parser.align)

    # Manually set the Clade (Mischko) colours to match Mischko et al. 2018
    if col == 'Clade (Mischko)':

        color_dict[1] = "#f88485"
        color_dict[2] = "#96b9da"
        color_dict[3] = "#c4e0a4"
        color_dict[4] = "#ffdf80"
        color_dict[5] = "#d17de8"
        color_dict[6] = "#FFA533"
        color_dict[7] = "#905E22"


    tree, ts = tc.get_example_tree(tree = tree, color_dict = color_dict, annot_dict=annot_dict, col=col)

    # Write to out path

    print ("\nTree image has been written to " + parser.outpath)
    tree.render(parser.outpath, tree_style=ts, dpi=300)
