

import itertools 

global chars
chars=["a","c","g","t"]



######################################################################################
#******function to compute max char in score matrix
def computeMax(l,start): 
    max_ind=0
    max_char=chars[0]
    max_val=l[0]
    for u in range(start,len(chars)):
        if l[u]>max_val:
            max_val=l[u]
            max_char=chars[u]
            max_ind=u
    

    return max_val,max_char,max_ind

##################################################################################################
#******function takes a sequence of dna and computes its score matrix and return suggested motif 
def scoreMat(dna_seq):
    
    #tanspose for easier computations
    dna_seq=[[dna_seq[j][i] for j in range(len(dna_seq))] for i in range(len(dna_seq[0]))]
    
    a_count=0
    c_count=0
    g_count=0
    t_count=0

    list=[]
    for i in dna_seq:
        for j in i:
            if j=="a":
                a_count+=1
            if j=="c":
                c_count+=1
            if j=="g":
                g_count+=1
            if j=="t":
                t_count+=1
        list.append(str(a_count)+str(c_count)+str(g_count)+str(t_count))
        a_count=0
        c_count=0
        g_count=0
        t_count=0
    count_matrix=[[list[j][i] for j in range(len(list))] for i in range(len(list[0]))]
    
    j=0
    sum_list=[]

    sg_motif=""
    cons_score=0
    for m in range(k):
        for i in range(len(count_matrix)):
            sum_list.append(int(count_matrix[i][j]))             
        same_val=[]
        max_val,max_char,max_index=computeMax(sum_list,0)   
        for i in range(max_index+1,len(sum_list)):
            if sum_list[i]==max_val:
                st=sg_motif
                st+=chars[sum_list[i]]
                same_val.append(st)    
        sg_motif+=max_char
        cons_score+=max_val
        sum_list=[]
        j+=1
    sg_motifs=[]
    #to handle motifs with same value in two or more chars ex: A&T both have same value in count matrix
    for j in range(len(same_val)):
        same_val[j]+=sg_motif[len(same_val[j]):]
        sg_motifs.append(same_val[j])
    cons_score=round((cons_score/(t*k))*100)
    sg_motifs.append(sg_motif)

    return count_matrix,sg_motifs,cons_score
    


##################################################################################
#*****takes a DNA, num of seqs in dna, length of dna seq, kmer length,and returns 
#***** list of its motifs, their score, and their starting posions

def motifSearch(dna,t,n,k):
    indexes=itertools.product(range(n-k), repeat=t)
    score_arr=[]
    motif_arr=[]
    startPos_arr=[]
    for i in indexes:
        v=0
        k_list=[]
        for j in i:
            p=dna[v][j:j+k]
            k_list.append(p)
            v+=1
        par_score=0
        sug_motifs=[]
        _,sug_motifs,par_score=scoreMat(k_list)
        #adding only unique motifs 
        for j in sug_motifs:
            if j not in motif_arr: 
                startPos_arr.append(i)
                score_arr.append(par_score)
                motif_arr.append(j)
       
    return score_arr,motif_arr,startPos_arr

###################################################################################
#**********MAIN*******************

DNA=[
    "ggctatccaggtactt",
     "gatccatacgtacacc",
    "aaacgttagtgcaccc"]
#lower bound to start our motifs selection
lower_bound=60

#length of the motif
k=4

#number of sequences in the DNA
t=len(DNA)

#length of each sequence in the DNA
n=len(DNA[0])

#arrays of all scores , motifs and starting postions
score_arr=[]
motif_arr=[]
s_arr=[]


score_arr,motif_arr,s_arr= motifSearch(DNA,t,n,k)
final_motifs=[]
final_scores=[]
final_s=[]
for i in range(len(score_arr)):
    if score_arr[i]>=lower_bound:
        final_scores.append(score_arr[i])
        final_motifs.append(motif_arr[i])
        final_s.append(s_arr[i])
for i in range(len(final_motifs)):
    print("motif : ",final_motifs[i], " with score: ",final_scores[i],"% ")
    print("starting positions : ",final_s[i])
    






