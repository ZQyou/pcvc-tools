from __future__ import print_function
import sys
import numpy as np

class _pcprofile(object):

    def __init__(self,pcmat,singletons=[],contigs=[]):
        self.pcmat = pcmat
        self.common_pc = pcmat*pcmat.transpose()
        self.singletons = singletons
        self.contigs = contigs

    def print_common_pc(self,print_mat=False):
        if print_mat:
            print("Common PC matrix:\n",self.common_pc)

    def print_percentage(self):
        print("Percentage (%):")
        n = self.contigs.size
        pct_mat = np.zeros((n,n),int)
        for i in xrange(0,n-1):
            for j in xrange(i+1,n):
                pct_mat[i,j] = 100*self.common_pc[i,j]/(self.common_pc[i,i] + [ 1 if singletons[i] > 0 else 0 ])
                pct_mat[j,i] = 100*self.common_pc[i,j]/(self.common_pc[j,j] + [ 1 if singletons[j] > 0 else 0 ])
        print(pct_mat)

        pct_file = open('pct_mat.csv','w')
        pct_file.write(",%s\n" % (",".join(self.contigs)))
        for i in xrange(0,n):
            pct_file.write("%s,%s\n" % (self.contigs[i], ",".join(str(x) for x in pct_mat[i,:])))
        pct_file.close();

def read_csv(csv_file):
    '''
    read in csv file { protein, contig, keywords, pc_id, vc_id }
    '''
    data = np.genfromtxt(csv_file,delimiter=',',dtype=None)
    #print(data.shape)
    #print(data)
    contigs = data[1:,1]
    pcs = data[1:,3]
    #vcs = data[1:,4]
    return [pcs, contigs]


def matrix_builder(pcs,contigs):
    '''
    build index matrix and get singletons
    '''
    singletons_b = pcs == np.array('')
    num_singletons = np.sum(singletons_b + [0])
    contigs_singletons = contigs[singletons_b]
    uniq_pcs = np.unique(pcs)
    uniq_pcs = uniq_pcs[ uniq_pcs != np.array('') ]
    uniq_contigs = np.unique(contigs)

    npc = uniq_pcs.size
    nv = uniq_contigs.size
    print("Number of Viruses = %d" % nv)
    print("Number of PCs = %d" % npc)
    print("Number of Singletons = %d" % num_singletons)

    pcmat = np.empty((0,npc),int)
    singletons = np.zeros(nv,int)
    contig_pc_mat = np.zeros((nv,npc+1),int)
    for v, contig in enumerate(np.nditer(uniq_contigs)):
        where_contig = np.where(np.in1d(contigs,contig))
        where_pc = np.in1d(uniq_pcs,pcs[where_contig])+[0]
        contig_pc_mat[v,:-1] = where_pc
        singletons[v] = sum((contigs_singletons == np.array(contig))+[0])
#       contig_pc_mat[v,npc] = 1 if singletons[v] > 0 else 0
        contig_pc_mat[v,npc] = singletons[v]
        pcmat = np.append(pcmat,np.array([where_pc]),axis=0) 

    print("Viruses:\n",uniq_contigs)
    print("PCs:\n",uniq_pcs)
    print("Singletons:\n",singletons)

    print("Contig-PC Matrix:\n",contig_pc_mat)
    f = open('contig_pc_mat.csv','w')
    f.write(",%s,Singletons\n" % (",".join(uniq_pcs)))
    for i in xrange(0,nv):
        f.write("%s,%s\n" % (uniq_contigs[i], ",".join(str(x) for x in contig_pc_mat[i,:])))
    f.close()


    return [np.matrix(pcmat),uniq_pcs,uniq_contigs,singletons]

if __name__ == '__main__':
    pcs_file = sys.argv[1];
    pcs, contigs = read_csv(pcs_file);
    pcmat, pcs, contigs, singletons = matrix_builder(pcs,contigs)
    pcprofile = _pcprofile(pcmat,singletons,contigs)
    pcprofile.print_common_pc(True)
    pcprofile.print_percentage()

