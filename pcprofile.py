from __future__ import print_function
import sys
import numpy as np

class _pcprofile(object):

    def __init__(self,pcmat):
        self.pcmat = pcmat
        self.common_pc = pcmat*pcmat.transpose()

    def print_common_pc(self):
        print("Common PC")
        print(self.common_pc)

    def 

def read_csv(csv_file):
    '''
    read in csv file { protein, contig, pcs, vcs }
    '''
    data = np.genfromtxt(csv_file,delimiter=',',dtype=None)
    contigs = data['f1'][1:]
    pcs = data['f2'][1:]
    #vcs = data['f3'][1:]
    return [pcs, contigs]


def matrix_builder(pcs,contigs):
    '''
    build index matrix and get singletons
    '''
    uniq_pcs = np.unique(pcs)
    n = uniq_pcs.size
    uniq_pcs = uniq_pcs[ uniq_pcs != np.array('') ]
    if n > uniq_pcs.size:
        print("PCs: %d vanished element(s) are removed" % (n-uniq_pcs.size))
    uniq_contigs = np.unique(contigs)
    n = uniq_contigs.size
    if n > uniq_contigs.size:
        print("Contigs: %d vanished element(s) are removed" % (n-uniq_contigs.size))

    pcmat = np.empty((0,uniq_pcs.size),int)
    #print(uniq_contigs)
    #print(uniq_pcs)
    for contig in np.nditer(uniq_contigs):
        where_contig = np.where(np.in1d(contigs,contig))
        where_pc = np.in1d(uniq_pcs,pcs[where_contig])+[0]
        pcmat = np.append(pcmat,np.array([where_pc]),axis=0) 

    return np.matrix(pcmat)

if __name__ == '__main__':
    pcs, contigs = read_csv("test.csv")
    pcmat = matrix_builder(pcs,contigs)
    pcprofile = _pcprofile(pcmat)
    pcprofile.print_common_pc()

