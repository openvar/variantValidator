#! /usr/bin/env python

from VariantValidator import mysql_refSeqGene_noMissmatch, compile_lrg_data

if __name__ == '__main__':

    # Update refSeqGene Primary assembly alignment data
    mysql_refSeqGene_noMissmatch.update()
    # Update LRG records
    compile_lrg_data.update()