#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

##########################################################################
#This code predicts whether similar sequences switch folds
#Authors: Allen Kim and Lauren Porter 
#Contact: Lauren Porter (lauren.porter@nih.gov)
#Inputs: Your email address (-e flag)
#        2 fasta files containing sequences to be compared.
#        The headers of the fastas should be formatted as:
#        > title|accession name|information about protein
#        Examples: GA.txt, GB.txt, RfaH_whole(or frag).txt
#Outputs: summary.txt
#Notes: The last field in summary.txt is the average secondary structure
#       discrepancy between the inputted seq (accession number in field 2)
#       and secondary structure predictions of all seqs in the file not
#       containing the inputted seq
###########################################################################

import warnings
warnings.filterwarnings('ignore', 'Could not import the lzma module.')
import argparse
import functools
import os
import pandas as pd
import re
import subprocess
import textwrap

# important information:
# * fasta file needs to be in the format given by uniprot

## path to programs and scripts
path_clustalo = '/Users/lauren.lea.porter/data/prediction_pipeline/clustalo' # path to clustalo binary
path_jpred4 = 'jpredapi' # path to jpred4 script
path_scheduler = 'multibatch_jpredapi.pl' # path to jpred scheduler


## default parameters
t_def = 0.1 # default threshold
# temp files should be name.temp.ext for clean to work properly
unaln_temp = 'unaln.temp.fasta'
aln_temp = 'aln.temp.fasta'
acc_temp = 'accessions.temp.txt'
jpred_temp = 'jpred_temp'
dld_temp = 'downloads_temp'


## internal use
## commands to execute
r_eh = re.compile('[EH]')
r_fasta = re.compile('>(\w+)\\|(\w+)\\|\w+')
r_input = re.compile('(?:.*/)*(.+)\\..+')
r_jout = re.compile('(.+).tar.gz')
r_jacc = re.compile('(\w+)\\.name')
r_jnet = re.compile('jnetpred:(.+)')
cmd_clustalo = path_clustalo+' --force --outfmt=fa -i \'{fasta:s}\' -o \'{fasta_aln:s}\''


# perform multiple sequence alignment
def run_msa(df):
        write_fasta(unaln_temp, df)
        run_clustalo(unaln_temp, aln_temp)
        df_aln = load_fasta(aln_temp, None)[['acc', 'seq']] \
                .rename(columns={'seq': 'seq_aln'}) \
                .set_index('acc')
        df = df.join(df_aln, on='acc')
        return df


# run clustal omega
def run_clustalo(fd_in, fd_out):
        cmd = cmd_clustalo.format(fasta=fd_in, fasta_aln=fd_out)
        subprocess.run(cmd, shell=True)


# write fasta file from dataframe
def write_fasta(fd, df):
        f = open(fd, 'w')
        for row in df.iterrows():
                f.write(row[1]['raw_h']+'\n')
                f.write(textwrap.fill(row[1]['seq'], width=60))
                f.write('\n')
        f.close()


# run pred scheduler
def run_scheduler(email):
        cmd = ('perl {scheduler:s} -j {jpred:s} -i {temp_dir:s} ' +\
                '-d {temp_dld:s} -f {accs:s} -e {email:s}') \
                .format(scheduler=path_scheduler, jpred=path_jpred4, \
                temp_dir=jpred_temp, temp_dld=dld_temp, accs=acc_temp, email=email)
        print(cmd)
        subprocess.call(cmd, shell=True)


# load uniprot-format fasta files onto dataframe
def load_fasta(fd, category):
        f = open(fd, 'r')
        seq = ''
        seqs = []
        accs = []
        raw_hs = []
        names = []
        for line in f:
                m = line.split('|')
                if len(m) > 1:
                        seqs.append(seq)
                        seq = ''
                        accs.append(m[1].strip())
                        names.append(m[2].strip())
                        raw_hs.append(line.strip())
                else:
                        seq += line.strip()
        f.close()
        seqs.append(seq)
        seqs = seqs[1:]
        df = pd.DataFrame({'acc': accs, 'name': names, 'category': category, \
                'seq': seqs, 'raw_h': raw_hs})
        return df


def prep_input_files(names):
        m = (r_input.match(names[0]), r_input.match(names[1]))

        if None in m:
                raise ValueError('File names must have extensions.')
        elif m[0].group(1) == m[1].group(1):
                raise ValueError('File names must be different.')

        # load fasta files
        fasta_a = load_fasta(names[0], m[0].group(1))
        fasta_a['category'] = m[0].group(1)
        fasta_b = load_fasta(names[1], m[1].group(1))
        fasta_b['category'] = m[1].group(1)

        return fasta_a.append(fasta_b).reset_index(drop=True)


def populate_seq_dir(df, n):
        # create fasta files
        try: 
                os.mkdir(jpred_temp)
        except FileExistsError:
                pass

        for row in df.iterrows():
                fd = jpred_temp+'/'+row[1]['acc']+'.fasta'
                f = open(fd, 'w')
                f.write(row[1]['raw_h']+'\n')
                if n != None:
                        f.write(textwrap.fill(row[1]['seq'][n:], width=60))
                else:
                        f.write(textwrap.fill(row[1]['seq'], width=60))
                f.close()


def generate_acc_file(df):
        f = open(acc_temp, 'w')
        for row in df.iterrows():
                f.write(row[1]['acc']+'\n')
        f.close()


def clean_temp():
        subprocess.run('rm -rf *_temp', shell=True)
        subprocess.run('rm -rf *.temp.*', shell=True)


def process_jpred_out(df):
        files = os.listdir(dld_temp)
        
        accs = []
        jids = []
        jpreds = []
        for i in files:
                m = r_jout.match(i)
                # create directory and unzip contents into it
                if (m == None):
                        continue
                try:
                        os.mkdir(dld_temp+'/'+m.group(1))
                except FileExistsError:
                        pass
                subprocess.run ('tar xfz {:s} -C {:s}'.format(  \
                        dld_temp+'/'+i, dld_temp+'/'+m.group(1)), \
                        shell=True)
        
                # get map between file name and list
                acc = map(r_jacc.match, os.listdir(dld_temp+'/'+m.group(1)))
                acc = functools.reduce(lambda x, y: x if y == None else y, acc)
                accs.append(acc.group(1))
                jids.append(m.group(1))

                # load predictions
                f_pred = open('{folder:s}/{jnetid:s}/{jnetid:s}.jnet' \
                        .format(folder=dld_temp, jnetid=m.group(1)))
                for line in f_pred:
                        m_jnet = r_jnet.match(line)
                        if m_jnet != None:
                                jpreds.append(m_jnet.group(1).replace(',', ''))
                f_pred.close()

        df_map = pd.DataFrame({'acc': accs, 'jpred_id': jids, \
                'jpred': jpreds}).set_index('acc')
        df = df.join(df_map, on='acc')
        return df

def alignSS(SS,aligned_seq):

    j = 0
    alignedSS = ''

    for i in range(len(aligned_seq)):
        
        if aligned_seq[i] != '-':
 
            if j >= len(SS):
                alignedSS += '-'
            
            elif SS[j] == '-':
                alignedSS += 'C'
            else:
                alignedSS += SS[j]
            j += 1

        else:
            alignedSS += '-'

    return alignedSS


def align_jpred(df, np):
        jpreds_aln = []
        for row in df.iterrows():
                jpred = row[1]['jpred']
                seq_aln = row[1]['seq_aln']

                # determine which direction to go
                #if np < 0:
                        #seq_aln = seq_aln[::-1]
                        #jpred = jpred[::-1]

                jpred_aln = alignSS(jpred,seq_aln)

                print(seq_aln)
                print(jpred_aln)

                # reverse if negative
                #if np < 0:
                        #jpred_aln = jpred_aln[::-1]
                jpreds_aln.append(jpred_aln)
        df['jpred_aln'] = jpreds_aln
        return df


def count_mismatches(seq_a, seq_b, np):

        if np < 0:
                _a = seq_a[::-1]
                _b = seq_b[::-1]
        i = 0
        j = 0
        s_disc = 0
        while i < len(_a) and j < len(_b):
                if (r_eh.match(_a[i]) != None) and (r_eh.match(_b[j]) != None):
                        if _a[i] != _b[j]:
                                s_disc += 1
                i += 1
                j += 1
        return s_disc

def calc_difference(ss1,ss2):

    l1 = len([x for x in ss1 if x in ['H','E']])
    l2 = len([x for x in ss2 if x in ['H','E']])

    diffs = [i for i in range(len(ss1)) if \
             ss1[i] != ss2[i] and \
             ss1[i] in ['H','E'] and ss2[i] in ['H','E']]

    return float(len(diffs))/min(l1,l2)


def calculate_stats(df, np):
        ss_diff = []
        ss_same = []
        for i in range(len(df)):
                s_diff = 0
                s_same = 0
                n_diff = 0
                n_same = 0
                for j in range(len(df)):
                        if i == j:
                                continue
                        mis = calc_difference(df.iloc[i]['jpred_aln'], \
                                df.iloc[j]['jpred_aln'])
                        if df.iloc[i]['category'] != df.iloc[j]['category']:
                                s_diff += mis
                                n_diff += 1
                        else:
                                s_same += mis
                                n_same += 1
                if n_diff == 0:
                        s_diff = 0
                else:
                        s_diff /= n_diff
                if n_same == 0:
                        s_same = 0
                else:
                        s_same /= n_same
                ss_diff.append(s_diff)
                ss_same.append(s_same)
        df['s_diff'] = ss_diff
        #df['s_same'] = ss_same
        return df


# main execution sequence
if __name__ == '__main__':
        # argument parser
        parser = argparse.ArgumentParser(description= \
                'Compares secondary structure prediction discrepancies '+ \
                'between two protein families.', \
                formatter_class=argparse.ArgumentDefaultsHelpFormatter )
        parser.add_argument('-t', type=float, metavar='THRES', \
                default=t_def, help='threshold for discrepancy')
        parser.add_argument('-np', type=int, metavar='N_PRED', \
                            default=None,help='number of residues to use for prediction'\
                            +'inputting sequences of desired lengths is recommended')
        parser.add_argument('-e', type=str, metavar='email', \
                help='email address for Jpred4 submission (required)')
        parser.add_argument('-c', nargs=2, type=str, \

                            metavar=('FASTA_A', 'FASTA_B'), \
                help='fasta files for secondary structure comparison (required)')
        parser.add_argument('--skip', action='store_const', const=True, \
                help ='skip scheduler')
        parser.add_argument('--clean', action='store_const', const=True, \
                help='clean temporary folders and files')
        args = parser.parse_args()

        # if no arguments are passed then...
        if args.clean == True:
                clean_temp()
                print('Temporary folders deleted.')
                if args.c == None:
                        exit()
        if args.c == None:
                parser.print_help()
                exit()

        # if we have arguments
        df = prep_input_files(args.c) # make sequence files for jpred inputs
        
        populate_seq_dir(df, args.np) # populate sequence directory
        generate_acc_file(df) # generate accession list for scheduler
        if args.skip != True:
                run_scheduler(args.e)
        df = process_jpred_out(df)
        df = run_msa(df)
        align_jpred(df, args.np)
        df = calculate_stats(df, args.np)
        df.to_csv('summary.txt')
