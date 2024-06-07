"""

#################### TO RUN EPITOPE PREDICTION ONLY ####################

python3 MMPRED.py -q QUERY -a ALLELE
    - QUERY (mandatory) is the fasta file to which prediction are applied
    - ALLELE (mandatory) is a 1-column txt file with the identifier of the alleles
    this script applies the CNNPEPPRED prediction for the alleles specified in ALLELE to the protein sequences in FASTA 


python3 MMPRED.py.py -q ... -a ... -n NETMHCIIPAN_PATH
    - NETMHCIIPAN_PATH (optional) is the relative or absolute path for the folder "mhc_ii" of the
       iedb prediction program. metti link and version ###############

python3 MHCIIPRED.py ... -m MODE
    - MODE (optional) can be either "protein" or "peptide". It specifies if the QUERY file contains full length
      protein sequences or small peptide/epitopes. 
      If MODE is set to "protein" the program will identity a 9-mer core for each window of size W of the sequences in QUERY
      If MODE is set to "peptide" the program predict one 9-mer core for each sequence in QUERY
      DEFAULT = protein

python3 MMPRED.py.py ... -m protein -w W
    - W (optional) is the window size for the prediction, default = 15. Is used only when MODE = protein


python3 MHCIIPRED.py ... -r RESULTS_FOLDER
    - RESULTS_FOLDER (optional) is the name of the resutls folder



##################### TO RUN ALIGNMENT AND EPITOPE PREDICTION #####################

python3 MMPRED.py.py -b BLAST_PATH -q QUERY -t TARGET -a ...
    - BLAST_PATH (mandatory) relative or absolute path to the blast folder "ncbi-blast-2.12.0+"
    - QUERY (mandatory) fasta file
    - TARGET (mandatory) fasta file 
    QUERY is aligned against TARGET. The TARGET's sequences that show a significant alignment 
    are then predicted has epitopes



python3 MMPRED.py -b BLAST_PATH -q QUERY -t TARGET -a ... -afp ALGN_FILT_PAR -afv ALGN_FILT_VAL
    - ALGN_FILT_PAR (optional) parameter choosen to filter alignment, possible values: "evalue", "bitscore". DEFAULT = evalue
    - ALGN_FILT_VAL (optional) is the cutoff to filter the alignments. DEFAULT = 0.05

python3 MMPRED.py ... -alg_mode ALG_MODE
	- ALGN_MODE (optional) can be either "blastp" or "psiblast"
	
python3 MMPRED.py ... -alg_mode psiblast -pssm_comp_db EPITOPE_DB
	- EPITOPE_DB (mandatory if psiblast is used) is the fasta file of the epitope dataset with which the pssm is computed


python3 MMPRED.py ... -n_core N
	- N is the maximum number of CPU to be used, parallel computation is used only when psiblast is applied
	


##################### TO RUN THE PIPELINE FROM THE PARAMETER FILE #####################

python3 MMPRED.py.py -getPF PARAM_FILE_NAME
    - PARAM_FILE_NAME: name of the empty parameter file
    An empty parameter named PARAM_FILE_NAME file is generated using this script. Istruction on how to use it are in the file itself

python3 MMPRED.py.py -PF PARAM_FILE_NAME
    runs the pipeline with the parameters specified in PARAM_FILE_NAME


"""






import sys
import os 
sys.path.append(os.getcwd()+'/bin/')
import shutil
import re
import pandas as pd
import numpy as np
from getOptPar import getOptPar
import subprocess
from model_initializer2 import CNNPepPred
from compute_hla_ranks import apply_ranks
import numpy as np
#from parseAlignmentOut import pars_alignment as get_human_hr_alignment
from getFastaFromAlignment import filter_alignment as getFastaFromAlignment
from fasta_utilities import read_fasta
from eliminateRedundancyFromFasta import eliminateRedundancy
from concurrent.futures import ThreadPoolExecutor


class MHCIIPRED:

    def get_input(self, iargs):
        idict = {  '-q':False, 
                   '-t':False, 
                   '-a':False, 
                   '-m':'protein', 
                   '-w':15, 
                   '-b':False, 
                   '-n':False,
                   '-PF':False,
                   '-getPF':False,
                   '-r':'./results/',
                   '-afp':'evalue',
                   '-afv':0.05, 
                   '-alg_mode':'blastp',
                   '-pssm_comp_db':False,
                   '-n_core':1}

        
        #check if input format is right
        if not(all([iargs[i][0] == '-' for i in range(1,len(iargs),2)])) or iargs[-1][0] == '-':
            print('error in command line format')
            return False


        # read the inputs
        for key in iargs:
            if key in ['-q', '-m', '-t', '-s', '-n', '-r', '-a', '-w', '-b', '-afp', '-afv', '-PF', '-getPF', '-alg_mode', '-pssm_comp_db', '-n_core']:
                if key in idict:
                    idict[key] = iargs[iargs.index(key)+1]
                else:
                    print(key+' is not a valid option')
                    return False


        self.idict = idict


        if self.idict['-getPF']:
            param_filename = self.idict['-getPF']
            shutil.copy('./bin/parameter_file.txt', param_filename)
            return False

        # - if the the parameter file is passed in input "-PF" parse it and extract idict
        if self.idict['-PF']:
            parameter_file = self.idict['-PF']
            self.parse_parameters_file(parameter_file)




        #print(self.idict)


        return True



  
    def run_analysis(self):



        
        # 1 - parsing the arguments
        wd = os.getcwd()
        query     = os.path.realpath(self.idict['-q'])
        
        target  = self.idict['-t']
        if target: 
            target = os.path.realpath(target)
            self.target = target

        # intersect the alleles in input with the one available with netmhciipan
        alleles_file = self.idict['-a']
        if not alleles_file:
            print('error: alleles must be specified ith the option -a')
            return
        alleles = [hla.rstrip().replace('HLA-', '').replace(' ','') for hla in open(alleles_file).readlines()]
        netmhciipan_modelled_alleles_file = './bin/netmhciipan_alleles.txt'
        modelled_alleles = [hla.rstrip() for hla in open(netmhciipan_modelled_alleles_file).readlines()]
        alleles_netmhciipan = list(set(alleles) & set(modelled_alleles))
        self.alleles_netmhciipan = alleles_netmhciipan



 
        
        # intersect the alleles in input with the one available with cnnpeppred
        alleles = [re.sub('\*|/|:', '_', hla.rstrip()).replace('HLA-', '') for hla in open(alleles_file).readlines()] 
        cnnpeppred_modelled_alleles_file = './bin/HLA_MODELS_TL/modelled_alleles.txt'
        modelled_alleles = [hla.rstrip() for hla in open(cnnpeppred_modelled_alleles_file)]
        alleles_cnnpeppred = list(set(alleles) & set(modelled_alleles))
        self.alleles_cnnpeppred = alleles_cnnpeppred








        results_path = self.idict['-r']
        results_path = os.path.realpath(results_path)
        if os.path.exists(results_path):
            shutil.rmtree(results_path)
        os.mkdir(results_path)
        self.results_path = results_path



        # 2 - running the analysis

        # if a blast path is specified, run the alignment 
        blast_path = self.idict['-b']
        if not os.path.exists(blast_path):
            raise IOError('path to blast does not exists')
        self.run_alignment = True
        blast_path=os.path.realpath(blast_path)
        print('RUNNING ALIGNMENT')
        if blast_path and self.idict['-alg_mode'] == 'blastp':
            query, any_seq_matching = self.runAlignment(blast_path, query, target, results_path)

        elif blast_path and self.idict['-alg_mode'] == 'psiblast':
            query, any_seq_matching = self.runAlignmentPssm(blast_path, query, target, self.idict['-pssm_comp_db'], results_path)


        if not any_seq_matching:
            print('No alignement satisified the filter: '+self.idict['-afp']+','+str(self.idict['-afv']))
            return
        pred_mode = self.idict['-m']
        W         = int(self.idict['-w'])


        pred_out = []


   

        netmhciipan_path = self.idict['-n']
        if netmhciipan_path:
            self.run_netmhciipan = True
            if os.path.exists(netmhciipan_path):
                self.alleles_netmhciipan = alleles_netmhciipan
                print('RUNNING_NETMHCII_BA')
                pred_out_netmhciipan_ba = self.runNetmhciiPanPred('netmhciipan_ba', alleles_netmhciipan, query, pred_mode, W, results_path)
                pred_out.extend(pred_out_netmhciipan_ba)
                
                print('RUNNING_NETMHCII_EL')
                pred_out_netmhciipan_el = self.runNetmhciiPanPred('netmhciipan_el', alleles_netmhciipan, query, pred_mode, W, results_path)
                pred_out.extend(pred_out_netmhciipan_el)
            else:
                print('path to netmhciipan does not exists')
                return
        else:
            self.run_netmhciipan = False
            
        #running cnnpeppred
        
        
        models_path = os.path.realpath(wd+'/bin/HLA_MODELS_TL/')
        rank_dist_comp = os.path.realpath(wd+'/bin/HLA_MODELS_TL/uniref50_rand_ep_score.txt')

        pred_out_cnnpeppred = self.runCnnPepPred(query, alleles_cnnpeppred, models_path, pred_mode, W, results_path, rank_dist_comp, gen_out_files=False)
        pred_out.extend(pred_out_cnnpeppred)
        

        self.pred_out = pred_out

        self.generate_summary()

        #generate summary

        return


    def parse_parameters_file(self, parameter_file):
        f = open(parameter_file)
        for line in f.readlines():
            line = line.rstrip()
            if line != '' and line[0] != '#':
                line = line.replace('"', '')
                [PAR, VAL] = line.split('=')


                if VAL != '':
                    if PAR == 'QUERY':
                        self.idict['-q'] = os.path.abspath(VAL)
                    elif PAR == 'ALLELES':
                        self.idict['-a'] = os.path.abspath(VAL)
                    elif PAR == 'PREDICTION_MODE':
                        self.idict['-m'] = VAL
                    elif PAR == 'W':
                        self.idict['-w'] = int(VAL)
                    elif PAR == 'NETMHCIIPANPATH':
                        self.idict['-n'] = os.path.abspath(VAL)
                    elif PAR == 'RESULTS_FOLDER':
                        self.idict['-r'] = VAL
                    elif PAR == 'TARGET':
                        self.idict['-t'] = VAL
                    elif PAR == 'BLASTPATH':
                        self.idict['-b'] = os.path.abspath(VAL)
                    elif PAR == 'ALGN_FILT_PAR':
                        self.idict['-afp'] = VAL
                    elif PAR == 'ALGN_FILT_THR':
                        self.idict['-afv'] = float(VAL)
                    elif PAR == 'ALGN_MODE':
                        self.idict['-alg_mode'] = VAL
                    elif PAR == 'PSSM_COMP_DB':
                        self.idict['-pssm_comp_db'] = os.path.abspath(VAL)
                    elif PAR == 'N_CORE':
                        self.idict['-n_core'] = int(VAL)
                        
                    else:
                        print([PAR, VAL])


        print(self.idict)

   

        return


    def runCnnPepPred(self, applyData, alleles, models_path, applyMode, W, savePath, randDistRankComp, gen_out_files=True):

        savePathTable = os.path.realpath(savePath)
        applyDataName, applyData = read_fasta(applyData, 'lists')

        if applyMode == 'protein':
            uniqLengths = [W]
            applyDataSubsets = [applyData]
            applyDataNameSubsets = [applyDataName]
            applyDataIdxSubsets = [[idx for idx in range(len(applyData))]]
            seqL = np.array([len(s) for s in applyData])
        elif applyMode == 'peptide':
            seqL = np.array([len(s) for s in applyData])
            uniqLengths = np.unique(seqL).tolist()
            uniqLengthsSeq = []
            applyDataSubsets = []
            applyDataNameSubsets = []
            applyDataIdxSubsets = []
            for L in uniqLengths:
                idxs = [idx for idx,l in enumerate(seqL) if l == L]
                applyDataSubsets.append([applyData[idx] for idx in idxs])
                applyDataNameSubsets.append([applyDataName[idx] for idx in idxs])
                applyDataIdxSubsets.append(idxs)
            

        
        #apply the models 

        doTraining = False
        trainingData = None
        trainingOutcome = None
        doLogoSeq = False
        doCV = False
        cvPart = None
        kFold = None
        doApplyData = True
        parametersFile = None
        


        print(alleles)
        colNames = ['Peptide_Source','Start','End','Peptide']
        colNames.extend(alleles)


        if applyMode == 'peptide':
            output_predscore = pd.DataFrame(index = range(len(applyData)), columns = colNames)
            output_core = pd.DataFrame(index = range(len(applyData)), columns = colNames)
        elif applyMode == 'protein':
            nRows = np.array([l-W+1 for l in seqL]).sum()
            output_predscore = pd.DataFrame(index = range(nRows), columns = colNames)
            output_core = pd.DataFrame(index = range(nRows), columns = colNames)



        model_names = os.listdir(models_path)


        for aln, al in enumerate(alleles):

            ###########

            print('predicting for allele '+al+' ('+str(aln+1)+'/'+str(len(alleles))+')')#+' '+str(aln+1)+'/'+str(len(alleles))+ ' done')
            
            allele = 'HLA_' + al

            curr_model = [m for m in model_names if m[6:m.find('_TL_')]==allele]        
            curr_model_path = models_path+'/'+curr_model[0]


            for idxL, curr_pep_len in enumerate(uniqLengths):
                applyDataCurrLen = applyDataSubsets[idxL]
                applyDataNameCurrLen = applyDataNameSubsets[idxL]




                modelCNN = CNNPepPred(allele,\
                                savePath,\
                                doTraining,\
                                trainingData,\
                                trainingOutcome,\
                                doLogoSeq,\
                                doCV,\
                                cvPart,\
                                kFold,\
                                doApplyData,\
                                curr_model_path,\
                                applyDataCurrLen,\
                                applyDataNameCurrLen,\
                                curr_pep_len,\
                                parametersFile)



                # ------------- # 2
                if curr_pep_len >= modelCNN.maxL:
                    #print(['--------------->', curr_pep_len, modelCNN.maxL])
                    break

                sIntApply,sApplyName = modelCNN.seq2Lmer(modelCNN.aa2int(modelCNN.applyData),L=None,nameSeq=modelCNN.applyDataName,saveLmer = True,takeUniqueLmer=False)[0:2]
                sIntApply = modelCNN.addEmptyPositions(sIntApply) 
                modelCNN.feedForwardAndGetScore(sIntApply,saveOutcome = True)
                modelCNN.getCoreBinder(modelCNN.int2aa(sIntApply),modelCNN.contributionScore,sApplyName,saveCoreBinders = True)
                modelCNN.printApplyOutcome()
                sApply = modelCNN.int2aa(modelCNN.applyDataSeq)
                nameApply = modelCNN.applyDataSeqName
                yhatApply = modelCNN.predictedOutcomeWithContributionScore
                yhatApply = yhatApply.round(3)
                coreBinders = modelCNN.coreBinders
                posStart = modelCNN.applyDataPositionStart 
                posEnd = modelCNN.applyDataPositionEnd
                modelCNN.save_object()




                # ------------- # 3



            
                if applyMode == 'peptide':
                    for idx, si in enumerate(applyDataIdxSubsets[idxL]):
                        output_predscore.loc[output_predscore.index[si], 'Peptide_Source'] = nameApply[idx]
                        output_predscore.loc[output_predscore.index[si], 'Start'] = posStart[idx]
                        output_predscore.loc[output_predscore.index[si], 'End'] = posEnd[idx]
                        output_predscore.loc[output_predscore.index[si], 'Peptide'] = sApply[idx]
                        output_predscore.loc[output_predscore.index[si], al] = yhatApply[idx]

                        output_core.loc[output_core.index[si], 'Peptide_Source'] = nameApply[idx]
                        output_core.loc[output_core.index[si], 'Start'] = posStart[idx]
                        output_core.loc[output_core.index[si], 'End'] = posEnd[idx]
                        output_core.loc[output_core.index[si], 'Peptide'] = sApply[idx]
                        output_core.loc[output_core.index[si], al] = coreBinders[idx]


                        
                
                elif applyMode == 'protein':
                    for si in range(nRows):
                        output_predscore.loc[output_predscore.index[si], 'Peptide_Source'] = nameApply[si]
                        output_predscore.loc[output_predscore.index[si], 'Start'] = posStart[si]
                        output_predscore.loc[output_predscore.index[si], 'End'] = posEnd[si]
                        output_predscore.loc[output_predscore.index[si], 'Peptide'] = sApply[si]
                        output_predscore.loc[output_predscore.index[si], al] = yhatApply[si]

                        output_core.loc[output_core.index[si], 'Peptide_Source'] = nameApply[si]
                        output_core.loc[output_core.index[si], 'Start'] = posStart[si]
                        output_core.loc[output_core.index[si], 'End'] = posEnd[si]
                        output_core.loc[output_core.index[si], 'Peptide'] = sApply[si]
                        output_core.loc[output_core.index[si], al] = coreBinders[si]
                        
                subprocess.Popen(('rm -r '+savePath+'/HLA_*').replace('//', '/'), shell=True).wait()

            


        # computing ranks
        # 1 reading random score distribution
        f = open(randDistRankComp)
        hlaRandScoreDict = {}
        for line in f.readlines():
            line=line.split(',')
            hla = line[0]
            randscores = np.array([float(x) for x in line[1:None]])
            hlaRandScoreDict[hla] = randscores


        # 2 generating the rank dataframe
        output_scorerank = pd.DataFrame(index = range(len(output_predscore[output_predscore.columns[0]])), columns = output_predscore.columns)
        for col_name in output_predscore.columns[0:4]:
            output_scorerank[col_name] = output_predscore[col_name]

        for hla in output_predscore.columns[4:None]:
            predscores = np.array(output_predscore[hla])
            ranks = np.zeros(len(predscores))
            randDist = hlaRandScoreDict[hla]
            for i,s in enumerate(predscores):
                if not np.isnan(s):
                    ranks[i] = np.round(np.sum(randDist >= s)/len(randDist)*100,3)
            output_scorerank[hla] = ranks



        #generating the output ########## RIVEDI BENE


        
        pred_out = []
        for ei in range(len(output_predscore[output_predscore.columns[0]])):
            for aj in range(4, len(output_predscore.columns)):
                
                seq_name = output_core['Peptide_Source'][ei]
                epitope = output_core['Peptide'][ei]
                start = int(output_core['Start'][ei])


                core = output_core[output_predscore.columns[aj]][ei]
                score = output_predscore[output_predscore.columns[aj]][ei]


                if applyMode == 'protein':
                    start = start + epitope.find(core)
                    end = start+8
                elif applyMode == 'peptide':
                    start = epitope.find(core)+1
                    end = start+8

                if not np.isnan(score):
                    rank = output_scorerank[output_predscore.columns[aj]][ei]
                    allele = output_core.columns[aj]
                    pred_out.append([core, seq_name, start, end, 'CNNPEPRED', score, rank, allele])

                
        return pred_out


    def runNetmhciiPanPred(self, method, alleles, fasta, predMode='peptide', W=15, savePath=None):


        if method not in ['netmhciipan_el', 'netmhciipan_ba']:
            raise Exception('method must be equal to netmhciipan_el or netmhciipan_ba')


        if type(alleles) == list:
            alleles = ','.join(alleles)


        seq_names = []
        seq_seq = []
        for line in open(fasta).readlines():
            if line[0]=='>':
                seq_names.append(line.rstrip())
                seq_seq.append('')
            else:
                seq_seq[-1] = seq_seq[-1]+line.rstrip()
        seq_len = [len(seq) for seq in seq_seq]





        pred_out = []
        netmhciipan_path = (self.idict['-n']+'/').replace('//', '/')
        tmp_out_fname = (savePath+'/netmhciipan_tmpout.txt').replace('//', '/')
        tmp_in_fname = (savePath+'/netmhciipan_tmpin.fasta').replace('//', '/')
        

        if predMode == 'peptide':
            uniqLen = list(set(seq_len))
            for leni, L in enumerate(uniqLen):
                # generate a netmhciipan_in_tmp.fasta with all sequences with current length L
                idx_curr_len =  [idx for idx in range(len(seq_len)) if seq_len[idx] == L]
                currL_fasta_cont = [seq_names[idx]+'\n'+seq_seq[idx] for idx in idx_curr_len]
                currL_fasta_cont = '\n'.join(currL_fasta_cont)+'\n'
                seq_names_currL = [seq_names[idx] for idx in idx_curr_len]
                open(tmp_in_fname, 'w').write(currL_fasta_cont)
                fasta = os.path.abspath(tmp_in_fname)

                

                #command = 'python3 {}mhc_II_binding.py {} {} {} {} > netmhciipan_tmpout'.format(netmhciipan_path, method, alleles, fasta, str(W))
                #process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                
                command = 'python3 {}mhc_II_binding.py {} {} {} {}'.format(netmhciipan_path, method, alleles, fasta, str(W))
                tmp_out_file = open(tmp_out_fname, 'w')
                process = subprocess.Popen(command, shell=True, stdout=tmp_out_file, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL).wait()
                

                for line in open(tmp_out_fname).readlines()[1:None]:
                    line = line.rstrip().split()
                    core = line[5]
                    peptide = line[6]
                    seq_name = seq_names_currL[int(line[1])-1]


                    start = peptide.find(core)+1
                    end = str(start+8)
                    start = str(start)

                    score = line[7]
                    rank = line[8]
                    allele = line[0]
                    pred_out.append([core, seq_name, start, end, method.upper(), score, rank, allele])
                os.remove(tmp_in_fname)
                os.remove(tmp_out_fname)
                
     

        elif predMode == 'protein':
            
            fasta = os.path.abspath(fasta)

            command = 'python3 {}mhc_II_binding.py {} {} {} {}'.format(netmhciipan_path, method, alleles, fasta, str(W))
            print(command)
            tmp_out_file = open(tmp_out_fname, 'w')
            process = subprocess.Popen(command, shell=True, stdout=tmp_out_file, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL).wait()

            #extracting output
            for line in open(tmp_out_fname).readlines()[1:None]:
                line = line.rstrip().split()

                core = line[5]
                peptide = line[6]
                start = int(line[2])+peptide.find(core)
                end = start+8




                seq_name = seq_names[int(line[1])-1]
                score = line[7]
                rank = line[8]
                allele = line[0]
                pred_out.append([core, seq_name, start, end, method.upper(), score, rank, allele])
            os.remove(tmp_out_fname)




        return pred_out        


    def runAlignment(self, blast_path, query, target, results_path):
        print([blast_path, query,  target])
        align_respath = results_path+'/alignment/'.replace('//','/')
        if os.path.exists(align_respath):
            shutil.rmtree(align_respath)
        os.mkdir(align_respath)
        blast_path = blast_path+'/bin/'.replace('//','/')


        #running the alignment
        command = '{}blastp -out {}.query_vs_target_alignment_fmt11.csv -outfmt 11 \
                -query {} -subject {} -ungapped -comp_based_stats F \
                -task blastp-short'.format(blast_path, align_respath, query, target)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()


        #reformatting the output 1
        command = '{}blast_formatter -archive {}.query_vs_target_alignment_fmt11.csv \
                -outfmt "10 delim=, qaccver qseq qstart qend saccver sseq sstart send bitscore pident length mismatch gapopen evalue"\
                -out {}query_vs_target_alignment.csv'.format(blast_path, align_respath, align_respath)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()



        #this function produce good ooutput only if the query is a fasta o 9-mer cores (maintain or selecte cores/noncores?)
        # metti in if insieme a repformatting input 2
        
        """
        #reformatting the output 2
        command = '{}blast_formatter -archive {}.query_vs_target_alignment_fmt11.csv \
                -outfmt 0\
                -out {}query_vs_target_alignment_hr.raw'.format(blast_path, align_respath, align_respath)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        get_human_hr_alignment('{}query_vs_target_alignment_hr.raw'.format(align_respath), \
                            '{}query_vs_target_alignment_hr.txt'.format(align_respath))
        """


        filter_par = self.idict['-afp']
        filter_value = float(self.idict['-afv'])
        any_seq_matching = getFastaFromAlignment('{}query_vs_target_alignment.csv'.format(align_respath),\
                                '{}target_seq_matching_query.fasta'.format(align_respath),\
                                 target, filter_par, filter_value, self.idict['-w'])


        #eliminate redundancy from matching sequences

        if any_seq_matching:

            self.alignment_redundaant_seq_dict = \
            eliminateRedundancy('{}target_seq_matching_query.fasta'.format(align_respath), 
                                '{}target_seq_matching_query_nr.fasta'.format(align_respath),
                                '{}fasta_redundancy_db.txt'.format(align_respath))

            query = '{}target_seq_matching_query_nr.fasta'.format(align_respath)

            #query = '{}target_seq_matching_query.fasta'.format(align_respath)
        return query, any_seq_matching




    def runAlignmentPssm(self, blast_path, query, target, epitopes_db, results_path):




        align_respath = results_path+'/alignment/'.replace('//','/')
        if os.path.exists(align_respath):
            shutil.rmtree(align_respath)
        os.mkdir(align_respath)
        blast_path = blast_path+'/bin/'.replace('//','/')
        
        query, any_match = self.pssm_ep_align(blast_path, query, target, epitopes_db)



        return query, any_match

    def pssm_ep_align(self, blast_path, query_fasta, target_fasta, epitope_db):


        # 1 - compute_pssm
        pssm_folder = self.compute_all_pssms(blast_path, query_fasta, target_fasta, epitope_db)
        alg_res_name = (self.idict['-r']+'/alignment/pssm_algn_query_vs_target.csv').replace('//', '/')

	files = glob.glob(f"{pssm_folder}*.csv")
        with open(alg_res_name, 'w') as outfile:
            for file in files:
                with open(file, 'r') as infile:
                    for line in infile:
                        if line.strip() and 'CONVERGED' not in line.upper():
                            outfile.write(line)


        filter_par = self.idict['-afp']
        filter_value = float(self.idict['-afv'])
        matching_seqs_file =(self.idict['-r']+'/alignment/target_seq_matching_query.fasta').replace('//', '/')
        any_seq_matching = getFastaFromAlignment(alg_res_name, matching_seqs_file,
                                 target_fasta, filter_par, filter_value, self.idict['-w'])


        #eliminate redundancy from matching sequences
        query_nr = matching_seqs_file.replace('_query.fasta', '_query_nr.fasta')
        if any_seq_matching:
            
            self.alignment_redundaant_seq_dict = \
            eliminateRedundancy(matching_seqs_file, 
                                query_nr,
                                self.idict['-r']+'/alignment/fasta_redundancy_db.txt')

            

        return query_nr, any_seq_matching

    def compute_all_pssms(self, blast_path, query_fasta, target_fasta, epitope_db):


        max_threads=self.idict['-n_core']


        pssm_folder =  (self.idict['-r']+'/alignment/PSSMs/').replace('//', '/')
        if os.path.exists(pssm_folder):
            shutil.rmtree(pssm_folder)
        os.makedirs(pssm_folder)


        self.sh_process_wrapper(f"{blast_path}makeblastdb -in {epitope_db} -dbtype prot")
        self.sh_process_wrapper(f"{blast_path}makeblastdb -in {target_fasta} -dbtype prot")
        query_eps_dict = self.read_fasta_to_dict(query_fasta)





        with ThreadPoolExecutor(max_workers=max_threads) as executor:
            futures = []
            for ep_name, ep_seq in query_eps_dict.items():
                futures.append(executor.submit(self.compute_pssm_n_align, blast_path, pssm_folder, ep_name, ep_seq, epitope_db, target_fasta))

            # Wait for all tasks to complete
            for future in futures:
                future.result()

        return pssm_folder



    def compute_pssm_n_align(self, blast_path, pssm_folder, ep_name, ep_seq, epitope_db, target_fasta):
        

        pssm_1_asn = os.path.join(pssm_folder, (ep_name+'_pssm1.asn').replace('|', '_'))
        pssm_1_ascii = pssm_1_asn.replace('_pssm1.asn', '_pssm1.ascii')
        currep_fasta_in = os.path.join(pssm_folder, (ep_name+'.fasta').replace('|', '_'))
        
        with open(currep_fasta_in, 'w') as f:
            f.write(f'>{ep_name}\n{ep_seq}\n')

        # 1 -running alignment against epitope_db
        
        command = f"{blast_path}psiblast -query {currep_fasta_in} -db {epitope_db} -inclusion_ethresh 0.05 \
                    -num_iterations 5 -out_pssm {pssm_1_asn}\
                    -pseudocount 1 \
                    -out_ascii_pssm {pssm_1_ascii}\
                    -save_pssm_after_last_round"
        stdout,_ = self.sh_process_wrapper(command)


 
        # 2 running alignment against target 


        command = f"{blast_path}psiblast -in_pssm {pssm_1_asn} -db {target_fasta}\
                    -num_iterations 1 -inclusion_ethresh 10\
                    -outfmt \"10 delim=, qaccver qseq qstart qend saccver sseq sstart  \
                    send bitscore pident length mismatch gapopen evalue\" \
                    -out {pssm_1_asn.replace('_pssm1.asn', '.csv')}"
        stdout,_  = self.sh_process_wrapper(command)










    def generate_summary(self):

        
        ranks = [float(x[6]) for x in self.pred_out]
        o = np.argsort(ranks)
        self.pred_out = [self.pred_out[i] for i in o[::-1]]

        eps = [x[0] for x in self.pred_out]
        o = np.argsort(eps)
        self.pred_out = [self.pred_out[i] for i in o]   


        f_out = open(self.results_path+'/PRED_SUMMARY.csv'.replace('//', '/'), 'w')
        if self.run_alignment:

            
            #read redundancy db
            id_redundant_ids_dict  = {}
            for line in open(self.results_path+'/alignment/fasta_redundancy_db.txt'.replace('//', '/')).readlines():
                line  = line.rstrip().split()
                id = line[0]
                fastas_ids = line[1].split(',')
                id_redundant_ids_dict[id] = fastas_ids


            f_out.write(','.join(('core', 'target_id', 'target_window_start', 'target_window_end', 'target_alg_start', 'target_alg_end', 'target_core_start', 'target_core_end', 
                            'ident', 'evalue', 'bitscore', 'target_aligned_seq', 'query_id', 'query_alg_start', 'query_alg_end', 'query_aligned_seq',
                             'method', 'score', 'rank', 'allele\n')))
            for pred in self.pred_out:
                core = pred[0]
                target_core_start = str(pred[2])
                target_core_end = str(pred[3])
                method, score, rank, allele = [str(x) for x in pred[4:8]]
                if float(rank) > 10:
                    continue
                allele = self.convert_allele_format(allele)
                for red_desc in id_redundant_ids_dict[pred[1][1::]]:
                    red_desc = red_desc.split('@@@')
                    target_id = red_desc[0].replace('>','')
                    target_window_start, target_window_end = target_id.split('|')[-2::]
                    target_id = '|'.join(target_id.split('|')[0:-2:])
                    ident, evalue, bitscore, target_aligned_seq, target_alg_start, target_alg_end, query_aligned_seq = red_desc[2].split('|')
                    query_id = red_desc[1]
                    query_alg_start, query_alg_end = query_id.split('|')[-2::]
                    query_id = '|'.join(query_id.split('|')[0:-2:])
                    
                    f_out.write(','.join((core, target_id, target_window_start, target_window_end, target_alg_start, target_alg_end, target_core_start, target_core_end, 
                            ident, evalue, bitscore, target_aligned_seq, query_id, query_alg_start, query_alg_end, query_aligned_seq,
                             method, score, rank, allele))+'\n')

        else:
            f_out.write(','.join(['core', 'query', 'start', 'end', 'method', 'score', 'rank', 'allele'])+'\n')
            for pred in self.pred_out:
                pred = [str(x) for x in pred]
                pred[1] = pred[1].replace('>', '')
                pred[7] = self.convert_allele_format(pred[7])
                f_out.write(','.join(pred)+'\n')
                               
                    





        
    def convert_allele_format(self, allele):
        allele = allele.replace('HLA-', '')
        symbs = '*:/*:'
        idxs = [4,7,10,15,18]
        allele = allele.replace('_',symbs[0],1)
        allele = allele.replace('_',symbs[1],1)
        allele = allele.replace('_',symbs[2],1)
        allele = allele.replace('_',symbs[3],1)
        allele = allele.replace('_',symbs[4],1)

        return allele



    def sh_process_wrapper(self, command, stdout_return=False, stderr_return=False):
        #print(command)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL)

        # Read the standard output and standard error
        stdout, stderr = process.communicate()

        # Decode the byte strings to UTF-8 and print them
        if stdout_return and stdout:
            print("Standard Output:")
            print(stdout.decode('utf-8'))
        if stderr_return and stderr:
            print("Standard Error:")
            print(stderr.decode('utf-8'))
            raise IOError(stderr.decode('utf-8'))
        process.wait()


        return stdout.decode('utf-8'), stderr.decode('utf-8')


    def read_fasta_to_dict(self, file_path):
        fasta_dict = {}
        with open(file_path, 'r') as fasta_file:
            seq_id = None
            sequence = ''
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    # If sequence ID is not None, store the previous sequence
                    if seq_id is not None:
                        fasta_dict[seq_id] = sequence
                    # Extract sequence ID
                    seq_id = line[1:]
                    sequence = ''
                else:
                    sequence += line
            # Store the last sequence
            if seq_id is not None:
                fasta_dict[seq_id] = sequence
        return fasta_dict




if __name__ == '__main__':
    iargs = sys.argv

    X = MHCIIPRED()
    if X.get_input(iargs):
        X.run_analysis()
