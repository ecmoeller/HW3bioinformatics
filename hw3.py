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


def main():
    print("In main")
    #Reading in data from txt file
    with open('Human_PPI.csv','r') as f:
        data = f.readlines()


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


if __name__ == '__main__':
    main()
