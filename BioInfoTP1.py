# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 09:36:57 2023

@author: Alexandre de Carvalho Assunção.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path


BLOSUM62_FILENAME = Path(__file__).with_name('blosum62.csv')
BLOSUM62 = pd.read_csv(BLOSUM62_FILENAME)

OUTPUT_FOLDER = Path(__file__).with_name('output')

'''
Constants for Gap penalties.
The values use were suggested by chatGPT as the standard ones in the
literature.
'''
GAP_OPENING_PENALTY = -10
GAP_WIDENING_PENALTY = -2


'''
AUXILIARY INPUT-OUTPUT FUNCTIONS:
'''

'''Reads first fasta file in the input folder'''    
def ParseFasta(mode):
    '''
    Input:
       ...
    Output:
        Names(String List) -> list containg the names of the proteins in file
        Descs(String List) -> list containg metadata of the proteins in file
        Sequences(String List) -> list containg the amino acid sequences of the
        proteins in file
    '''
    os.chdir(Path(__file__).with_name('input'))
    files = os.listdir(os.getcwd())
    
    Names = []
    Descs = []
    Sequences = []
    seq = ""
    
    data = open(files[0],'r')
    
    number_of_proteins = 0
        
    for line in data:
        
        if line.startswith('>'):
            
            if len(seq) > 0:
                Sequences.append("".join(seq.split('\n')))
                seq = ""
                
            metadata = line.split('|') 
            Names.append(metadata[1])
            Descs.append(metadata[2])
            number_of_proteins += 1
            
        elif line[0] in 'ACDEFGHIKLMNPQRSTVWXY':
            seq = seq + line
            if mode == 'P' and number_of_proteins > 2:
                Names.pop()
                Descs.pop()
                break
        else:
            print("IO Error.")
            
        

    Sequences.append("".join(seq.split('\n')))    
    
    data.close()
    os.chdir('..')

    return Names,Descs,Sequences

'''Generates fasta file with aligned sequences'''
def GenFasta(names,descs,seqs):
    
    '''
    Input:
        names(String List) -> list containg the names of the proteins to be printed
        descs(String List) -> list containg metadata of the proteins to be printed
        seqs(String List) -> list containg the amino acid sequences of the proteins
                             to be printed.
    Output:
        None
    '''
    os.chdir(Path(__file__).with_name('output'))
    
    f = open('Aligned_Proteins.fasta','w')
    
    for i in range(len(names)):
        f.write('>sp|'+names[i]+'|'+descs[i])
        f.write(seqs[i])
        f.write('\n')
    
    f.close()
    os.chdir('..')

'''
CLASSES USED BY THE PROJECT:
     --> Protein
     --> Alignment
'''
class Protein(object):
    
    def __init__(self,P_name,P_desc,P_seq):
        
        '''
        Initializes a Protein object
        
        A Protein object has two attributes:
        name(String) -> the protein's name
        desc(String) -> protein metadata from the fasta file
        seq(String) -> the protein's amino acid sequence abreviated with capital letters.
        size(Int) -> number os amino acids in the protein
        '''
        self.name = P_name
        self.desc = P_desc
        self.seq = P_seq
        self.size = len(P_seq)
        
    def Set_Seq(self,s):
        
        self.seq = s
        
    def get_ProteinName(self):
        
        '''
        Used to safely access the protein's name
        Returns self.P_name
        '''
        return self.name
    
    
    def get_ProteinSequence(self):
        
        '''
        Used to safely access the protein's amino acid sequence
        Returns self.P_Seq
        '''
        return self.seq
    
    
    def __repr__(self):
        
       return f'Protein(\'{self.name}\', {self.seq})'
           
class PairwiseGlobalAlignment(Protein):
    
    def __init__(self,P1,P2):
        
        self.P1 = P1
        self.P2 = P2
        self.score = 0
        self.alignment = ("","")
        
    def Set_Score(self):
        
        self.score,self.alignment = Align_Protein(self.P1, self.P2)
        
    def get_AlignedP1(self):
            
        return self.alignment[0]
        
    def get_AlignedP2(self):
            
        return self.alignment[1]
    
    def Get_Protein1(self):
        
        return self.P1.get_ProteinName()
    
    def Get_Protein2(self):
        
        return self.P2.get_ProteinName()
    
    def Get_Score(self):
        
        return self.score
    
    def __str__(self):
        
        print('Alignment between '+self.P1.get_ProteinName()+' and '+self.P2.get_ProteinName()+':')
        print('Score = '+str(self.Get_Score()))
        return 'Aligned '+ self.P1.get_ProteinName()+': '+self.get_AlignedP1()+'\nAligned '+self.P2.get_ProteinName()+': '+self.get_AlignedP2()
       
'''
MAIN FUNCTION OF THE PROJECT:
    THE NEEDLEMAN-WUNSCH DYNAMIC PROGRAMING ALGORITHM
'''
def Align_Protein(P1,P2):
      
    '''
    Simple implementation of the Needleman-Wunsch algorithm
    to align to sequences of amino acid residues in a polypeptide
    chain
    
    INPUTS:
        P1(Protein Object)
        P2(Protein Object)
    
    OUTPUTS:
        max_score(Int) --> Optimal aligment score using BLOSUM62.
        best_path(Int List) --> List of integers representing the optimal path.
    '''
    
    #Step1 --> Initialization of scores matrix and paths matrix.
    
    m,n = P1.size+1,P2.size+1
    
    scores = np.zeros((m,n), dtype=int)
    path = np.zeros((m,n), dtype=int)
    
    scores[0][0] = 0
    path[0][0] = -1
    for i in range(1,m):
        scores[i][0] = GAP_OPENING_PENALTY + (i-1)*(GAP_WIDENING_PENALTY)
        path[i][0] = 1
    for i in range(1,n):
        scores[0][i] = GAP_OPENING_PENALTY + (i-1)*(GAP_WIDENING_PENALTY)
        path[0][i] = 2
    
    #Step2 --> Recursevely fill the scores matrix and the paths matrix
    
    for i in range(1,m):
       for j in range(1,n):
           
           #---------------------------------------------------------------------------------------
           #Value of score[i][j] if we're coming from the diagonal score[i-1][j-1]
           #---------------------------------------------------------------------------------------
           Match_Mismatch = scores[i-1][j-1] + BLOSUM62[P1.seq[i-1]][P2.seq[j-1]]
           #---------------------------------------------------------------------------------------
           #----Value of score[i][j] if we either widen the gap by coming from the left (i,j-1) or
           #----open a new gap by coming from above: i-1,j.
           #---------------------------------------------------------------------------------------
           Pro1_Del = max(scores[i][j-1] + GAP_WIDENING_PENALTY, scores[i-1][j] + GAP_OPENING_PENALTY)
           #---------------------------------------------------------------------------------------
           #----Value of score[i][j] if we either widen the gap by coming from above (i-1,j) or
           #----open a new gap by coming from above: i,j-1.
           #---------------------------------------------------------------------------------------
           Pro2_Del = max(scores[i-1][j] + GAP_WIDENING_PENALTY, scores[i][j-1] + GAP_OPENING_PENALTY)
           #---------------------------------------------------------------------------------------
           #----Put the above values into a list and then pick the maximum score.
           #----Record the path taken as well.
           #---------------------------------------------------------------------------------------
           Record_Step = [Match_Mismatch, Pro1_Del, Pro2_Del]
           scores[i][j] = max(Record_Step)
           path[i][j] = np.argmax(Record_Step) 
           
    #-----------------------------------------------------#
    #---Recovers the best path via the path matrix--------#
    #-----------------------------------------------------#
    
    i = m-1
    j = n-1
    best_path = []
    
    while i>0 or j>0:
        
        best_path.append(path[i][j])
        
        if best_path[-1] == 0:
            i = i-1
            j = j-1
        elif best_path[-1] == 1 and j > 0:
            j = j-1
        elif best_path[-1] == 2 and i > 0:
            i = i-1
        else:
            print("Invalid Path")
            return
        
    #-----------------------------------------------------#
    #---Sets the highest score----------------------------#
    #-----------------------------------------------------#
    
    score = scores[P1.size][P2.size] 
    
    #-----------------------------------------------------#
    #---Defines the best alignment------------------------#
    #-----------------------------------------------------#
    s1 = ""
    s2 = ""
    
    rev_path = best_path[::-1]
    i = 0 
    j = 0
    
    for A in rev_path:
        
        if A == 0:
            s1 = s1 + str(P1.get_ProteinSequence()[i])
            s2 = s2 + str(P2.get_ProteinSequence()[j])
            i = i+1
            j = j+1
        elif A == 1:
            s1 = s1 + '-'
            s2 = s2 + str(P2.get_ProteinSequence()[j])
            j = j+1
        elif A == 2:
           s1 = s1 + str(P1.get_ProteinSequence()[i])
           s2 = s2 + '-'
           i = i+1
        else:
            print("Invalid Alignment Path")
            
    alignment = (s1,s2)
    
    return score,alignment
    
'''
FUNCTIONS FOR THE MULTIPLE ALIGNMENT HEURISTIC
'''  
def Is_Adjacent(E1,E2):
    
   '''
   Checks if edges E1 and E2 share one AND ONLY one node
   '''
   L1 = [E1.Get_Protein1(),E1.Get_Protein2()]
   L2 = [E2.Get_Protein1(),E2.Get_Protein2()]
   
   count = 0
   
   for n in L1:
       if n in L2:
           count += 1
           
   if count == 1:
       return True
   else:
       return False
      
def Forms_Loop(P,A):
    
    if A.Get_Protein1() in P and A.Get_Protein2() in P:
        return True
    else:
        return False
    
def GuidingTree(N,G):
    
  
    Path = []
    Nodes_in_path = []
    
    Node_Visits = {}
    
    for node in N:
        Node_Visits[node] = 0
    
    Current_min = G[0]
    for edge in G[1:]:
        if edge.Get_Score() > Current_min.Get_Score():
            Current_min = edge
            
   
    Path.append(Current_min)
    
    
    Nodes_in_path.append(Current_min.Get_Protein1())
    Nodes_in_path.append(Current_min.Get_Protein2())
    
    Node_Visits[Current_min.Get_Protein1()] += 1 
    Node_Visits[Current_min.Get_Protein2()] += 1 
   
    '''
    Main Loop of the algorithm:
    Create an Adjacency list of the current edge we have, and then greedily select the one
    with the best score. This will be done until an edge has an empty adjacency list.
    '''
    
    Next = Current_min
    
    while(True):
        
        
        Curr = Next
        Adjacency = []
        for edge in G:
            if Is_Adjacent(Curr, edge) and not Forms_Loop(Nodes_in_path,edge):
                if Node_Visits[edge.Get_Protein1()] < 2 and Node_Visits[edge.Get_Protein2()] < 2:
                    Adjacency.append(edge)
         
        
        '''
        (iii) Use a greedy strategy to select the next edge to add to path
        '''
        Next = Adjacency[0]
        for edge in Adjacency[1:]:
            if edge.Get_Score() > Curr.Get_Score():
                Next = edge
                
              
        if Next.Get_Protein1() not in Nodes_in_path:
            Nodes_in_path.append(Next.Get_Protein1())
            Node_Visits[Next.Get_Protein2()] += 1
            
            
        if Next.Get_Protein2() not in Nodes_in_path:
            Nodes_in_path.append(Next.Get_Protein2())
            Node_Visits[Next.Get_Protein1()] += 1 
            
        '''
        Add the edge to the circuit
        '''
        Path.append(Next)
        '''
        Check if circuit is complete
        '''
        if len(Nodes_in_path) == len(N):
            break
      
    return(Nodes_in_path)
    
def Eulerian_Multiple_Alignment():
    
    '''
    INPUT:
        fasta_in(String) -> Path to a fasta file with the proteins to align.
        fasta_out(String) -> Path to where you want to save the fasta file with the aligned proteins.
    OUTPUT:
    '''
    
    n,d,s = ParseFasta('M')
    
    
    Pro_ids = range(len(n))
    Proteins = []
    Pro_Dict = {}
    Pairwise_Al = []
    
    
    for i in Pro_ids:
        Proteins.append(Protein(n[i],d[i],s[i]))
        Pro_Dict[n[i]] = (d[i],s[i])
        
    for i in Pro_ids:
        for j in Pro_ids[i:]:   
            if i == j:
                continue
            else:
                Pairwise_Al.append(PairwiseGlobalAlignment(Proteins[i],Proteins[j]))
                
    for A in Pairwise_Al:
        print("aligning " + A.Get_Protein1() + " and " + A.Get_Protein2())
        A.Set_Score()
        
    Ordering = GuidingTree(n,Pairwise_Al)
    
    Ordered_Proteins = []    
        
    Ordered_Proteins.append(Protein(Ordering[0],Pro_Dict[Ordering[0]][0],Pro_Dict[Ordering[0]][1]))
    
    for i in range(1,len(Ordering),2):
        
        P_Temp1 = Protein(Ordering[i-1],Pro_Dict[Ordering[i-1]][0],Pro_Dict[Ordering[i-1]][1])
        P_Temp2 = Protein(Ordering[i],Pro_Dict[Ordering[i]][0],Pro_Dict[Ordering[i]][1])
        
        print('Realigning '+P_Temp1.get_ProteinName()+' and '+P_Temp2.get_ProteinName())
        s,a = Align_Protein(P_Temp1, P_Temp2)
        
        Ordered_Proteins.append(P_Temp1)
        Ordered_Proteins[-1].Set_Seq(a[0])
        Ordered_Proteins.append(P_Temp2)
        Ordered_Proteins[-1].Set_Seq(a[1])
    
    out_seq = []
    for i in Ordered_Proteins:
        out_seq.append(i.get_ProteinSequence())
    
    GenFasta(n,d,out_seq)
    
    

'''CORE OF THE PROGRAM:

Opens a fasta file in the input folder that is in the same
directory as this file.

IMPORTANT: If there are multiple fasta files in the input directory/
           the program will run the first element of that list.
'''

def AlignSequenceFromFasta():

    #  n <-- names of the proteins in the fasta file
    #  d <-- description of the proteins in the fasta file
    #  s <-- amino acid sequences of the proteins in the fasta file
    n,d,s = ParseFasta('P')
    
    '''
    Gets input from the user to decide whether to use the standard NW pairwise/
    algorithm or the multiple alignment sequence
    '''
    
    print('Hello!')
    
    mode = input('Type P for pairwise alignment or M for multiple alignment:\n')
    
    if mode == 'P' or mode == 'p':
        
        P1 = Protein(n[0],d[0],s[0])
        P2 = Protein(n[1],d[1],s[1])
        P_Alignment = PairwiseGlobalAlignment(P1, P2)
        print("\nAligning "+n[0]+" and "+n[1])
        P_Alignment.Set_Score()
        
        '''Generating the output file'''
        if not os.path.exists(OUTPUT_FOLDER):
            os.makedir(OUTPUT_FOLDER)
        GenFasta(n,d,[P_Alignment.get_AlignedP1(),P_Alignment.get_AlignedP2()])
        
        print("\nA fasta file containg the aligned sequences (gaps included) of "+n[0]+" and "+n[1]+" is in the output folder.\n")
        print_flag = input("If you wish to print the NW score and sequences in the prompt, type 'Y'. If not, type any other character.\n")
        if print_flag == 'Y' or print_flag == 'y':
            print(P_Alignment)
            
    if mode == 'M' or mode == 'm':
        
        Eulerian_Multiple_Alignment()



if __name__ == '__main__':
    
    AlignSequenceFromFasta()
    
    print("\nDo you Wish to see the alignment in biodash app?\n")
    print("WARNING: Requires dashbio, dash, and an internet browser!\n")
    view_app = input("Type 'Y' if yes or any other character for no:")
    print("\n")
    
    if view_app == 'Y' or view_app == 'y':
        
       
        from dash import Dash, html, Input, Output, callback
        import dash_bio as dashbio
        
        os.chdir(Path(__file__).with_name('output'))
        
        f = open('Aligned_Proteins.fasta','r')
        d = f.read()
        
        f.close()
        os.chdir('..')
        
        app = Dash(__name__)
        
        app.layout = html.Div([
            dashbio.AlignmentChart(
                id='my-default-alignment-viewer',
                data=d,
                height=900,
                tilewidth=30,
            ),
            html.Div(id='default-alignment-viewer-output')
        ])

        @callback(
            Output('default-alignment-viewer-output', 'children'),
            Input('my-default-alignment-viewer', 'eventDatum')
        )
        def update_output(value):
            if value is None:
                return 'No data.'
            return str(value)
        
        app.run(debug=True, port=8051)
        print("\nTo see the alignment, open your brower and go to:\n--> http://127.0.0.1:8051/")
