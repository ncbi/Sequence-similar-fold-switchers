# Sequence-similar-fold-switchers
Code and data for predicting mutation-induced protein fold switching

Dependencies: clustal omega: http://www.clustal.org/omega/
              python3: https://www.python.org/downloads/
              pandas module for python3: https://www.python.org/downloads/
              
Scripts included: jpredapi
                  multibatch_jpredapi
                  
NOTES:
    Before running, open the source code (compare_preds.py) and:
        ->Change paths for clustalo, Jpred4 API, and scheduler as appropariate
        
usage: compare_test.py [-h] [-t THRES] [-np N_PRED] [-e email] [-c FASTA_A FASTA_B] [--skip] [--clean]

Compares secondary structure prediction discrepancies between two protein families.

optional arguments:
  -h, --help          show this help message and exit
  -t THRES            threshold for discrepancy (default: 0.1)
  -np N_PRED          number of residues to use for prediction (default: 200)
  -e email            email address for Jpred4 submission (required) (default: None)
  -c FASTA_A FASTA_B  fasta files for secondary structure comparison (required) (default: None)
  --skip              skip scheduler (default: None)
  --clean             clean temporary folders and files (default: None)
        
        

        
    
        
                  
