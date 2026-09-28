"""tabkit: A toolkit for sequence similarity analysis based on tabular data"""

import os
import sys
from Bio import SeqIO
import multiprocessing
import argparse
import itertools
from progress.bar import ChargingBar
import re
import subprocess
from pathlib import Path
import numpy as np
import shutil
import pandas as pd
import networkx as nx
import math
import datetime
from .tabtk import tabtkMain

def main():
    try:
        tabtkMain()            
    except KeyboardInterrupt:
        print("Execulted halted")

if __name__=="__main__":
        main()
