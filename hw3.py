import math
from operator import itemgetter 

import plotly
plotly.tools.set_credentials_file(username='moell162', api_key='XWeDQVjHzfvEN01WsE7p')

import plotly.plotly as py
import plotly.graph_objs as go

# Create random data with numpy
import numpy as np

class Pair:
    def __init__(self, proteinA, proteinB):
        self.proteinA = proteinA
        self.proteinB = proteinB  


def q2(pairs):
    # Measure the degree of each protein with at least 1 interaction in the network (exclude selfinteractions)
    degs = {}

    for p in pairs:

        #exclude self interactions
        if(p.proteinA != p.proteinB):
            #Increment proteinA count in degs if found, otherwise add it
            if (p.proteinA not in degs):
                degs[p.proteinA] = 1
            else:
                degs[p.proteinA] = degs[p.proteinA] + 1

            #Increment proteinB count in degs if found, otherwise add it
            if (p.proteinB not in degs):
                degs[p.proteinB] = 1
            else:
                degs[p.proteinB] = degs[p.proteinB] + 1        

    for k, v in degs.items():
        print(k, v)

    # Plot the degree distribution of the protein-protein interaction network (a histogram is fine)

def main():
    print("In main")
    #Reading in data from txt file
    with open('Human_PPI.csv','r') as f:
        data = f.readlines()

    with open("Lit_degrees.csv", "r") as f2:
        data2 = f2.readlines()

    count = 0
    pairs = []
    for line in data:
        line = line.strip()
        if (count == 0):
            print("did nothing")
        if(count >= 1 ):
            pairterms = line.split('\t')
            pair = Pair(pairterms[0], pairterms[1])
            pairs.append(pair)

        count += 1

    count = 0
    degrees = {}
    for line in data2:
        line = line.strip()
        if (count == 0):
            print("did nothing")
        if(count >= 1 ):
            degreeLine = line.split('\t')
            degrees[degreeLine[0]] = degreeLine[1]

        count += 1

    print("First line ", degrees["A1CF"])

    q2(pairs)

if __name__ == '__main__':
    main()
