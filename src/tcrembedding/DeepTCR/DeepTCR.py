import sys
sys.path.append('../')
from DeepTCR.functions.Layers import *
from DeepTCR.functions.utils_s import *
from DeepTCR.functions.act_fun import *
from DeepTCR.functions.plot_func import *
import glob
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, MultiLabelBinarizer
from multiprocessing import Pool
import pickle
import shutil

class DeepTCR_base(object):

    def __init__(self,Name,max_length=40,device=0,tf_verbosity=3):
        """
        # Initialize Training Object.
        Initializes object and sets initial parameters.

        All DeepTCR algorithms begin with initializing a training object. This object will contain all methods, data, and results during the training process. One can extract learned features, per-sequence predictions, among other outputs from DeepTCR and use those in their own analyses as well.

        This method is included in the three main DeepTCR objects:

        - DeepTCR_U (unsupervised)
        - DeepTCR_SS (supervised sequence classifier/regressor)
        - DeepTCR_WF (supervised repertoire classifier/regressor)


        Args:
            Name (str): Name of the object. This name will be used to create folders with results as well as a folder with parameters and specifications for any models built/trained.

            max_length (int): maximum length of CDR3 sequence.

            device (int): In the case user is using tensorflow-gpu, one can specify the particular device to build the graphs on. This selects which GPU the user wants to put the graph and train on.

            tf_verbosity (str): determines how much tensorflow log output to display while training.
            0 = all messages are logged (default behavior)
            1 = INFO messages are not printed
            2 = INFO and WARNING messages are not printed
            3 = INFO, WARNING, and ERROR messages are not printed

        """

        #Assign parameters
        self.Name = Name
        self.max_length = max_length
        self.use_beta = False
        self.use_alpha = False
        self.device = 'cpu'
        #self.device = '/device:GPU:'+str(device)
        self.use_v_beta = False
        self.use_d_beta = False
        self.use_j_beta = False
        self.use_v_alpha = False
        self.use_j_alpha = False
        self.use_hla = False
        self.use_hla_sup = False
        self.keep_non_supertype_alleles = False
        self.regression = False
        self.use_w = False
        self.ind = None
        self.unknown_str = '__unknown__'

        #Create dataframes for assigning AA to ints
        aa_idx, aa_mat = make_aa_df()
        aa_idx_inv = {v: k for k, v in aa_idx.items()}
        self.aa_idx = aa_idx
        self.aa_idx_inv = aa_idx_inv

        tf.compat.v1.disable_eager_execution()
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = str(tf_verbosity)

        # Create directory for results of analysis
        directory = os.path.join(self.Name, 'results')
        self.directory_results = directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Create directory for any temporary files
        directory = self.Name
        if not os.path.exists(directory):
            os.makedirs(directory)

        tf.compat.v1.disable_eager_execution()
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = str(tf_verbosity)

    def Get_Data(self,directory,Load_Prev_Data=False,classes=None,type_of_data_cut='Fraction_Response',data_cut=1.0,n_jobs=1,
                    aa_column_alpha = None,aa_column_beta = None, count_column = None,sep='\t',aggregate_by_aa=True,
                    v_alpha_column=None,j_alpha_column=None,
                    v_beta_column=None,j_beta_column=None,d_beta_column=None,
                 p=None,hla=None,use_hla_supertype=False,keep_non_supertype_alleles=False):
        """
        # Get Data for DeepTCR
        Parse Data into appropriate inputs for neural network from directories where data is stored.

        This method can be used when your data is stored in directories and you want to load it from directoreis into DeepTCR. This method takes care of all pre-processing of the data including:

         - Combining all CDR3 sequences with the same nucleotide sequence (optional).
         - Removing any sequences with non-IUPAC characters.
         - Removing any sequences that are longer than the max_length set when initializing the training object.
         - Determining how much of the data per file to use (type_of_data_cut)
         - Whether to use HLA/HLA-supertypes during training.

        This method is included in the three main DeepTCR objects:

        - DeepTCR_U (unsupervised)
        - DeepTCR_SS (supervised sequence classifier/regressor)
        - DeepTCR_WF (supervised repertoire classifier/regressor)

        Args:
            directory (str): Path to directory with folders with tsv/csv files are present for analysis. Folders names become           labels for files within them. If the directory contains the TCRSeq files not organized into classes/labels,             DeepTCR will load all files within that directory.

            Load_Prev_Data (bool): Loads Previous Data. This allows the user to run the method once, and then set this parameter to True to reload the data from a local pickle file.

            classes (list): Optional selection of input of which sub-directories to use for analysis.

            type_of_data_cut (str): Method by which one wants to sample from the TCRSeq File.

                ###
                Options are:

                - Fraction_Response: A fraction (0 - 1) that samples the top fraction of the file by reads. For example, if one wants to sample the top 25% of reads, one would use this threshold with a data_cut = 0.25. The idea of this sampling is akin to sampling a fraction of cells from the file.

                - Frequency_Cut: If one wants to select clones above a given frequency threshold, one would use this threshold. For example, if one wanted to only use clones about 1%, one would enter a data_cut value of 0.01.

                - Num_Seq: If one wants to take the top N number of clones, one would use this threshold. For example, if one wanted to select the top 10 amino acid clones from each file, they would enter a data_cut value of 10.

                - Read_Cut: If one wants to take amino acid clones with at least a certain number of reads, one would use this threshold. For example, if one wanted to only use clones with at least 10 reads,they would enter a data_cut value of 10.

                - Read_Sum: IF one wants to take a given number of reads from each file, one would use this threshold. For example, if one wants to use the sequences comprising the top 100 reads of hte file, they would enter a data_cut value of 100.

            data_cut (float or int): Value  associated with type_of_data_cut parameter.

            n_jobs (int): Number of processes to use for parallelized operations.

            aa_column_alpha (int): Column where alpha chain amino acid data is stored. (0-indexed).

            aa_column_beta (int): Column where beta chain amino acid data is stored.(0-indexed)

            count_column (int): Column where counts are stored.

            sep (str): Type of delimiter used in file with TCRSeq data.

            aggregate_by_aa (bool): Choose to aggregate sequences by unique amino-acid. Defaults to True. If set to False, will allow duplicates of the same amino acid sequence given it comes from different nucleotide clones.

            v_alpha_column (int): Column where v_alpha gene information is stored.

            j_alpha_column (int): Column where j_alpha gene information is stored.

            v_beta_column (int): Column where v_beta gene information is stored.

            d_beta_column (int): Column where d_beta gene information is stored.

            j_beta_column (int): Column where j_beta gene information is stored.

            p (multiprocessing pool object): For parellelized operations, one can pass a multiprocessing pool object to this method.

            hla (str): In order to use HLA information as part of the TCR-seq representation, one can provide a csv file where the first column is the file name and the remaining columns hold HLA alleles for each file. By including HLA information for each repertoire being analyzed, one is able to find a representation of TCR-Seq that is more meaningful across repertoires with different HLA backgrounds.

            use_hla_supertype (bool): Given the diversity of the HLA-loci, training with a full allele may cause over-fitting. And while individuals may have different HLA alleles, these different allelees may bind peptide in a functionality similar way. This idea of supertypes of HLA is a method by which assignments of HLA genes can be aggregated to 6 HLA-A and 6 HLA-B supertypes. In roder to convert input of HLA-allele genes to supertypes, a more biologically functional representation, one can se this parameter to True and if the alleles provided are of one of 945 alleles found in the reference below, it will be assigned to a known supertype.

                - For this method to work, alleles must be provided in the following format: A0101 where the first letter of the designation is the HLA loci (A or B) and then the 4 digit gene designation. HLA supertypes only exist for HLA-A and HLA-B. All other alleles will be dropped from the analysis.

                - Sidney, J., Peters, B., Frahm, N., Brander, C., & Sette, A. (2008). HLA class I supertypes: a revised and updated classification. BMC immunology, 9(1), 1.

            keep_non_supertype_alleles (bool): If assigning supertypes to HLA alleles, one can choose to keep HLA-alleles that do not have a known supertype (i.e. HLA-C alleles or certain HLA-A or HLA-B alleles) or discard them for the analysis. In order to keep these alleles, one should set this parameter to True. Default is False and non HLA-A or B alleles will be discarded.

        Returns:
            variables into training object

            - self.alpha_sequences (ndarray): array with alpha sequences (if provided)
            - self.beta_sequences (ndarray): array with beta sequences (if provided)
            - self.class_id (ndarray): array with sequence class labels
            - self.sample_id (ndarray): array with sequence file labels
            - self.freq (ndarray): array with sequence frequencies from samples
            - self.counts (ndarray): array with sequence counts from samples
            - self.(v/d/j)_(alpha/beta) (ndarray): array with sequence (v/d/j)-(alpha/beta) usage

        """

        if Load_Prev_Data is False:

            if aa_column_alpha is not None:
                self.use_alpha = True

            if aa_column_beta is not None:
                self.use_beta = True

            if v_alpha_column is not None:
                self.use_v_alpha = True

            if j_alpha_column is not None:
                self.use_j_alpha = True

            if v_beta_column is not None:
                self.use_v_beta = True

            if d_beta_column is not None:
                self.use_d_beta = True

            if j_beta_column is not None:
                self.use_j_beta = True


            #Determine classes based on directory names
            data_in_dirs = True
            if classes is None:
                classes = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory,d))]
                classes = [f for f in classes if not f.startswith('.')]
                if not classes:
                    classes = ['None']
                    data_in_dirs = False


            self.lb = LabelEncoder()
            self.lb.fit(classes)
            self.classes = self.lb.classes_

            if p is None:
                p_ = Pool(n_jobs)
            else:
                p_ = p

            if sep == '\t':
                ext = '*.tsv'
            elif sep == ',':
                ext = '*.csv'
            else:
                print('Not Valid Delimiter')
                return

            #Get data from tcr-seq files
            alpha_sequences = []
            beta_sequences = []
            v_beta = []
            d_beta = []
            j_beta = []
            v_alpha = []
            j_alpha = []
            label_id = []
            file_id = []
            freq = []
            counts=[]
            file_list = []
            seq_index = []
            print('Loading Data...')
            for type in self.classes:
                if data_in_dirs:
                    files_read = glob.glob(os.path.join(directory, type, ext))
                else:
                    files_read = glob.glob(os.path.join(directory,ext))
                num_ins = len(files_read)
                args = list(zip(files_read,
                                [type_of_data_cut] * num_ins,
                                [data_cut] * num_ins,
                                [aa_column_alpha] * num_ins,
                                [aa_column_beta] * num_ins,
                                [count_column] * num_ins,
                                [sep] * num_ins,
                                [self.max_length]*num_ins,
                                [aggregate_by_aa]*num_ins,
                                [v_beta_column]*num_ins,
                                [d_beta_column]*num_ins,
                                [j_beta_column]*num_ins,
                                [v_alpha_column]*num_ins,
                                [j_alpha_column]*num_ins))

                DF = p_.starmap(Get_DF_Data, args)

                DF_temp = []
                files_read_temp = []
                for df,file in zip(DF,files_read):
                    if df.empty is False:
                        DF_temp.append(df)
                        files_read_temp.append(file)

                DF = DF_temp
                files_read = files_read_temp

                for df, file in zip(DF, files_read):
                    if aa_column_alpha is not None:
                        alpha_sequences += df['alpha'].tolist()
                    if aa_column_beta is not None:
                        beta_sequences += df['beta'].tolist()

                    if v_alpha_column is not None:
                        v_alpha += df['v_alpha'].tolist()

                    if j_alpha_column is not None:
                        j_alpha += df['j_alpha'].tolist()

                    if v_beta_column is not None:
                        v_beta += df['v_beta'].tolist()

                    if d_beta_column is not None:
                        d_beta += df['d_beta'].tolist()

                    if j_beta_column is not None:
                        j_beta += df['j_beta'].tolist()

                    label_id += [type] * len(df)
                    file_id += [file.split('/')[-1]] * len(df)
                    file_list.append(file.split('/')[-1])
                    freq += df['Frequency'].tolist()
                    counts += df['counts'].tolist()
                    seq_index += df.index.tolist()

            alpha_sequences = np.asarray(alpha_sequences)
            beta_sequences = np.asarray(beta_sequences)
            v_beta = np.asarray(v_beta)
            d_beta = np.asarray(d_beta)
            j_beta = np.asarray(j_beta)
            v_alpha = np.asarray(v_alpha)
            j_alpha = np.asarray(j_alpha)
            label_id = np.asarray(label_id)
            file_id = np.asarray(file_id)
            freq = np.asarray(freq)
            counts = np.asarray(counts)
            seq_index = np.asarray(seq_index)

            Y = self.lb.transform(label_id)
            OH = OneHotEncoder(sparse=False,categories='auto')
            Y = OH.fit_transform(Y.reshape(-1,1))

            print('Embedding Sequences...')
            #transform sequences into numerical space
            if aa_column_alpha is not None:
                args = list(zip(alpha_sequences, [self.aa_idx] * len(alpha_sequences), [self.max_length] * len(alpha_sequences)))
                result = p_.starmap(Embed_Seq_Num, args)
                sequences_num = np.vstack(result)
                X_Seq_alpha = np.expand_dims(sequences_num, 1)

            if aa_column_beta is not None:
                args = list(zip(beta_sequences, [self.aa_idx] * len(beta_sequences),  [self.max_length] * len(beta_sequences)))
                result = p_.starmap(Embed_Seq_Num, args)
                sequences_num = np.vstack(result)
                X_Seq_beta = np.expand_dims(sequences_num, 1)


            if p is None:
                p_.close()
                p_.join()

            if self.use_alpha is False:
                X_Seq_alpha = np.zeros(shape=[len(label_id)])
                alpha_sequences = np.asarray([None]*len(label_id))

            if self.use_beta is False:
                X_Seq_beta = np.zeros(shape=[len(label_id)])
                beta_sequences = np.asarray([None]*len(label_id))

            #transform v/d/j genes into categorical space
            num_seq = X_Seq_alpha.shape[0]
            if self.use_v_beta is True:
                self.lb_v_beta = LabelEncoder()
                self.lb_v_beta.classes_ = np.insert(np.unique(v_beta), 0, self.unknown_str)
                v_beta_num = self.lb_v_beta.transform(v_beta)
            else:
                self.lb_v_beta = LabelEncoder()
                v_beta_num = np.zeros(shape=[num_seq])
                v_beta = np.asarray([None]*len(label_id))

            if self.use_d_beta is True:
                self.lb_d_beta = LabelEncoder()
                self.lb_d_beta.classes_ = np.insert(np.unique(d_beta), 0, self.unknown_str)
                d_beta_num = self.lb_d_beta.transform(d_beta)
            else:
                self.lb_d_beta = LabelEncoder()
                d_beta_num = np.zeros(shape=[num_seq])
                d_beta = np.asarray([None]*len(label_id))

            if self.use_j_beta is True:
                self.lb_j_beta = LabelEncoder()
                self.lb_j_beta.classes_ = np.insert(np.unique(j_beta), 0, self.unknown_str)
                j_beta_num = self.lb_j_beta.transform(j_beta)
            else:
                self.lb_j_beta = LabelEncoder()
                j_beta_num = np.zeros(shape=[num_seq])
                j_beta = np.asarray([None]*len(label_id))

            if self.use_v_alpha is True:
                self.lb_v_alpha = LabelEncoder()
                self.lb_v_alpha.classes_ = np.insert(np.unique(v_alpha), 0, self.unknown_str)
                v_alpha_num = self.lb_v_alpha.transform(v_alpha)
            else:
                self.lb_v_alpha = LabelEncoder()
                v_alpha_num = np.zeros(shape=[num_seq])
                v_alpha = np.asarray([None]*len(label_id))

            if self.use_j_alpha is True:
                self.lb_j_alpha = LabelEncoder()
                self.lb_j_alpha.classes_ = np.insert(np.unique(j_alpha), 0, self.unknown_str)
                j_alpha_num = self.lb_j_alpha.transform(j_alpha)
            else:
                self.lb_j_alpha = LabelEncoder()
                j_alpha_num = np.zeros(shape=[num_seq])
                j_alpha = np.asarray([None]*len(label_id))

            if hla is not None:
                self.use_hla = True
                hla_df = pd.read_csv(hla)
                if use_hla_supertype:
                    hla_df = supertype_conv(hla_df,keep_non_supertype_alleles)
                    self.use_hla_sup = True
                    self.keep_non_supertype_alleles = keep_non_supertype_alleles
                hla_df = hla_df.set_index(hla_df.columns[0])
                hla_id = []
                hla_data = []
                for i in hla_df.iterrows():
                    hla_id.append(i[0])
                    temp = np.asarray(i[1].dropna().tolist())
                    hla_data.append(temp)

                hla_id = np.asarray(hla_id)
                hla_data = np.asarray(hla_data)

                keep,idx_1,idx_2 = np.intersect1d(file_list,hla_id,return_indices=True)
                file_list = keep
                hla_data = hla_data[idx_2]

                self.lb_hla = MultiLabelBinarizer()
                hla_data_num = self.lb_hla.fit_transform(hla_data)

                hla_data_seq_num = np.zeros(shape=[file_id.shape[0],hla_data_num.shape[1]])
                for file,h in zip(file_list,hla_data_num):
                    hla_data_seq_num[file_id==file] = h
                hla_data_seq_num = hla_data_seq_num.astype(int)
                hla_data_seq = np.asarray(self.lb_hla.inverse_transform(hla_data_seq_num))

                #remove sequences with no hla information
                idx_keep = np.sum(hla_data_seq_num,-1)>0
                X_Seq_alpha = X_Seq_alpha[idx_keep]
                X_Seq_beta = X_Seq_beta[idx_keep]
                Y = Y[idx_keep]
                alpha_sequences = alpha_sequences[idx_keep]
                beta_sequences = beta_sequences[idx_keep]
                label_id = label_id[idx_keep]
                file_id = file_id[idx_keep]
                freq = freq[idx_keep]
                counts = counts[idx_keep]
                seq_index = seq_index[idx_keep]
                v_beta = v_beta[idx_keep]
                d_beta = d_beta[idx_keep]
                j_beta = j_beta[idx_keep]
                v_alpha = v_alpha[idx_keep]
                j_alpha = j_alpha[idx_keep]
                v_beta_num = v_beta_num[idx_keep]
                d_beta_num = d_beta_num[idx_keep]
                j_beta_num = j_beta_num[idx_keep]
                v_alpha_num = v_alpha_num[idx_keep]
                j_alpha_num = j_alpha_num[idx_keep]
                hla_data_seq = hla_data_seq[idx_keep]
                hla_data_seq_num = hla_data_seq_num[idx_keep]

            else:
                self.lb_hla = MultiLabelBinarizer()
                file_list = np.asarray(file_list)
                hla_data = np.asarray(['None']*len(file_list))
                hla_data_num = np.asarray(['None']*len(file_list))
                hla_data_seq = np.asarray(['None']*len(file_id))
                hla_data_seq_num = np.asarray(['None']*len(file_id))

            with open(os.path.join(self.Name,'Data.pkl'), 'wb') as f:
                pickle.dump([X_Seq_alpha,X_Seq_beta,Y, alpha_sequences,beta_sequences, label_id, file_id, freq,counts,seq_index,
                             self.lb,file_list,self.use_alpha,self.use_beta,
                             self.lb_v_beta, self.lb_d_beta, self.lb_j_beta,self.lb_v_alpha,self.lb_j_alpha,
                             v_beta, d_beta,j_beta,v_alpha,j_alpha,
                             v_beta_num, d_beta_num, j_beta_num,v_alpha_num,j_alpha_num,
                             self.use_v_beta,self.use_d_beta,self.use_j_beta,self.use_v_alpha,self.use_j_alpha,
                             self.lb_hla, hla_data, hla_data_num,hla_data_seq,hla_data_seq_num,
                             self.use_hla,self.use_hla_sup,self.keep_non_supertype_alleles],f,protocol=4)

        else:
            with open(os.path.join(self.Name,'Data.pkl'), 'rb') as f:
                X_Seq_alpha,X_Seq_beta,Y, alpha_sequences,beta_sequences, label_id, file_id, freq,counts,seq_index,\
                self.lb,file_list,self.use_alpha,self.use_beta,\
                    self.lb_v_beta, self.lb_d_beta, self.lb_j_beta,self.lb_v_alpha,self.lb_j_alpha,\
                    v_beta, d_beta,j_beta,v_alpha,j_alpha,\
                    v_beta_num, d_beta_num, j_beta_num,v_alpha_num,j_alpha_num,\
                    self.use_v_beta,self.use_d_beta,self.use_j_beta,self.use_v_alpha,self.use_j_alpha,\
                    self.lb_hla, hla_data,hla_data_num,hla_data_seq,hla_data_seq_num,\
                self.use_hla,self.use_hla_sup,self.keep_non_supertype_alleles = pickle.load(f)

        self.X_Seq_alpha = X_Seq_alpha
        self.X_Seq_beta = X_Seq_beta
        self.Y = Y
        self.alpha_sequences = alpha_sequences
        self.beta_sequences = beta_sequences
        self.class_id = label_id
        self.sample_id = file_id
        self.freq = freq
        self.counts = counts
        self.sample_list = file_list
        self.v_beta = v_beta
        self.v_beta_num = v_beta_num
        self.d_beta = d_beta
        self.d_beta_num = d_beta_num
        self.j_beta = j_beta
        self.j_beta_num = j_beta_num
        self.v_alpha = v_alpha
        self.v_alpha_num = v_alpha_num
        self.j_alpha = j_alpha
        self.j_alpha_num = j_alpha_num
        self.seq_index = np.asarray(list(range(len(self.Y))))
        self.predicted = np.zeros((len(self.Y),len(self.lb.classes_)))
        self.hla_data_seq = hla_data_seq
        self.hla_data_seq_num = hla_data_seq_num
        self.w = np.ones(len(self.seq_index))
        #self.seq_index_j = seq_index
        print('Data Loaded')

    def Load_Data(self,alpha_sequences=None,beta_sequences=None,v_beta=None,d_beta=None,j_beta=None,
                  v_alpha=None,j_alpha=None,class_labels=None,sample_labels=None,freq=None,counts=None,Y=None,
                  p=None,hla=None,use_hla_supertype=False,keep_non_supertype_alleles=False,w=None):
        """
        # Load Data programatically into DeepTCR.

        DeepTCR allows direct user input of sequence data for DeepTCR analysis. By using this method,
        a user can load numpy arrays with relevant TCRSeq data for analysis.

        Tip: One can load data with the Get_Data command from directories and then reload it into another DeepTCR object with the Load_Data command. This can be useful, for example, if you have different labels you want to train to, and you need to change the label programatically between training each model. In this case, one can load the data first with the Get_Data method and then assign the labels pythonically before feeding them into the DeepTCR object with the Load_Data method.

        Of note, this method DOES NOT combine sequences with the same amino acid sequence. Therefore, if one wants this, one should first do it programatically before feeding the data into DeepTCR with this method.

        Another special use case of this method would be for any type of regression task (sequence or repertoire models). In the case that a per-sequence value is fed into DeepTCR (with Y), this value either becomes the per-sequence regression value or the average of all Y over a sample becomes the per-sample regression value. This is another case where one might want to load data with the Get_Data method and then reload it into DeepTCR with regression values.

        This method is included in the three main DeepTCR objects:

        - DeepTCR_U (unsupervised)
        - DeepTCR_SS (supervised sequence classifier/regressor)
        - DeepTCR_WF (supervised repertoire classifier/regressor)

        Args:

            alpha_sequences (ndarray of strings): A 1d array with the sequences for inference for the alpha chain.

            beta_sequences (ndarray of strings): A 1d array with the sequences for inference for the beta chain.

            v_beta (ndarray of strings): A 1d array with the v-beta genes for inference.

            d_beta (ndarray of strings): A 1d array with the d-beta genes for inference.

            j_beta (ndarray of strings): A 1d array with the j-beta genes for inference.

            v_alpha (ndarray of strings): A 1d array with the v-alpha genes for inference.

            j_alpha (ndarray of strings): A 1d array with the j-alpha genes for inference.

            class_labels (ndarray of strings): A 1d array with class labels for the sequence (i.e. antigen-specificities)

            sample_labels (ndarray of strings): A 1d array with sample labels for the sequence. (i.e. when loading data from different samples)

            counts (ndarray of ints): A 1d array with the counts for each sequence, in the case they come from samples.

            freq (ndarray of float values): A 1d array with the frequencies for each sequence, in the case they come from samples.

            Y (ndarray of float values): In the case one wants to regress TCR sequences or repertoires against a numerical label, one can provide these numerical values for this input. For the TCR sequence regressor, each sequence will be regressed to the value denoted for each sequence. For the TCR repertoire regressor, the average of all instance level values will be used to regress the sample. Therefore, if there is one sample level value for regression, one would just repeat that same value for all the instances/sequences of the sample.

            hla (ndarray of tuples/arrays): To input the hla context for each sequence fed into DeepTCR, this will need to formatted as an ndarray that is (N,) for each sequence where each entry is a tuple or array of strings referring to the alleles seen for that sequence. ('A*01:01', 'A*11:01', 'B*35:01', 'B*35:02', 'C*04:01')

            use_hla_supertype (bool): Given the diversity of the HLA-loci, training with a full allele may cause over-fitting. And while individuals may have different HLA alleles, these different allelees may bind peptide in a functionality similar way. This idea of supertypes of HLA is a method by which assignments of HLA genes can be aggregated to 6 HLA-A and 6 HLA-B supertypes. In roder to convert input of HLA-allele genes to supertypes, a more biologically functional representation, one can se this parameter to True and if the alleles provided are of one of 945 alleles found in the reference below, it will be assigned to a known supertype.

                - For this method to work, alleles must be provided in the following format: A0101 where the first letter of the designation is the HLA loci (A or B) and then the 4 digit gene designation. HLA supertypes only exist for HLA-A and HLA-B. All other alleles will be dropped from the analysis.

                - Sidney, J., Peters, B., Frahm, N., Brander, C., & Sette, A. (2008). HLA class I supertypes: a revised and updated classification. BMC immunology, 9(1), 1.

            keep_non_supertype_alleles (bool): If assigning supertypes to HLA alleles, one can choose to keep HLA-alleles that do not have a known supertype (i.e. HLA-C alleles or certain HLA-A or HLA-B alleles) or discard them for the analysis. In order to keep these alleles, one should set this parameter to True. Default is False and non HLA-A or B alleles will be discarded.

            p (multiprocessing pool object): a pre-formed pool object can be passed to method for multiprocessing tasks.

            w (ndarray): optional set of weights for training of autoencoder

        Returns:
            variables into training object

            - self.alpha_sequences (ndarray): array with alpha sequences (if provided)
            - self.beta_sequences (ndarray): array with beta sequences (if provided)
            - self.label_id (ndarray): array with sequence class labels
            - self.file_id (ndarray): array with sequence file labels
            - self.freq (ndarray): array with sequence frequencies from samples
            - self.counts (ndarray): array with sequence counts from samples
            - self.(v/d/j)_(alpha/beta) (ndarray):array with sequence (v/d/j)-(alpha/beta) usage

        ---------------------------------------

        """

        inputs = [alpha_sequences,beta_sequences,v_beta,d_beta,j_beta,v_alpha,j_alpha,
                  class_labels,sample_labels,counts,freq,Y,hla,w]

        for i in inputs:
            if i is not None:
                assert isinstance(i,np.ndarray),'Inputs into DeepTCR must come in as numpy arrays!'

        inputs = [alpha_sequences,beta_sequences,v_beta,d_beta,j_beta,v_alpha,j_alpha]
        for i in inputs:
            if i is not None:
                len_input = len(i)
                break

        if p is None:
            #p_ = Pool(40)
            p_ = Pool(1)
        else:
            p_ = p

        if alpha_sequences is not None:
            self.alpha_sequences = alpha_sequences
            args = list(zip(alpha_sequences, [self.aa_idx] * len(alpha_sequences), [self.max_length] * len(alpha_sequences)))
            result = p_.starmap(Embed_Seq_Num, args)
            sequences_num = np.vstack(result)
            self.X_Seq_alpha = np.expand_dims(sequences_num, 1)
            self.use_alpha = True
        else:
            self.X_Seq_alpha = np.zeros(shape=[len_input])
            self.alpha_sequences = np.asarray([None] * len_input)

        if beta_sequences is not None:
            self.beta_sequences = beta_sequences
            args = list(zip(beta_sequences, [self.aa_idx] * len(beta_sequences), [self.max_length] * len(beta_sequences)))
            result = p_.starmap(Embed_Seq_Num, args)
            sequences_num = np.vstack(result)
            self.X_Seq_beta = np.expand_dims(sequences_num, 1)
            self.use_beta = True
        else:
            self.X_Seq_beta = np.zeros(shape=[len_input])
            self.beta_sequences = np.asarray([None] * len_input)

        if v_beta is not None:
            self.v_beta = v_beta
            self.lb_v_beta = LabelEncoder()
            self.lb_v_beta.classes_ = np.insert(np.unique(v_beta), 0, self.unknown_str)
            self.v_beta_num = self.lb_v_beta.transform(v_beta)
            self.use_v_beta = True
        else:
            self.lb_v_beta = LabelEncoder()
            self.v_beta_num = np.zeros(shape=[len_input])
            self.v_beta = np.asarray([None] * len_input)

        if d_beta is not None:
            self.d_beta = d_beta
            self.lb_d_beta = LabelEncoder()
            self.lb_d_beta.classes_ = np.insert(np.unique(d_beta), 0, self.unknown_str)
            self.d_beta_num = self.lb_d_beta.transform(d_beta)
            self.use_d_beta = True
        else:
            self.lb_d_beta = LabelEncoder()
            self.d_beta_num = np.zeros(shape=[len_input])
            self.d_beta = np.asarray([None] * len_input)

        if j_beta is not None:
            self.j_beta = j_beta
            self.lb_j_beta = LabelEncoder()
            self.lb_j_beta.classes_ = np.insert(np.unique(j_beta), 0, self.unknown_str)
            self.j_beta_num = self.lb_j_beta.transform(j_beta)
            self.use_j_beta = True
        else:
            self.lb_j_beta = LabelEncoder()
            self.j_beta_num = np.zeros(shape=[len_input])
            self.j_beta = np.asarray([None] * len_input)

        if v_alpha is not None:
            self.v_alpha = v_alpha
            self.lb_v_alpha = LabelEncoder()
            self.lb_v_alpha.classes_ = np.insert(np.unique(v_alpha), 0, self.unknown_str)
            self.v_alpha_num = self.lb_v_alpha.transform(v_alpha)
            self.use_v_alpha = True
        else:
            self.lb_v_alpha = LabelEncoder()
            self.v_alpha_num = np.zeros(shape=[len_input])
            self.v_alpha = np.asarray([None] * len_input)

        if j_alpha is not None:
            self.j_alpha = j_alpha
            self.lb_j_alpha = LabelEncoder()
            self.lb_j_alpha.classes_ = np.insert(np.unique(j_alpha), 0, self.unknown_str)
            self.j_alpha_num = self.lb_j_alpha.transform(j_alpha)
            self.use_j_alpha = True
        else:
            self.lb_j_alpha = LabelEncoder()
            self.j_alpha_num = np.zeros(shape=[len_input])
            self.j_alpha = np.asarray([None] * len_input)

        if p is None:
            p_.close()
            p_.join()

        if counts is not None:
            if sample_labels is not None:
                count_dict={}
                for s in np.unique(sample_labels):
                    idx = sample_labels==s
                    count_dict[s]=np.sum(counts[idx])

                freq = []
                for c,n in zip(counts,sample_labels):
                    freq.append(c/count_dict[n])
                freq = np.asarray(freq)
                self.counts = counts
            else:
                print('Counts need to be provided with sample labels')
                return

        if freq is not None:
            self.freq = freq

        if sample_labels is not None:
            self.sample_id = sample_labels
        else:
            self.sample_id = np.asarray(['None']*len_input)

        if class_labels is not None:
            self.class_id = class_labels
        else:
            self.class_id = np.asarray(['None']*len_input)

        if (counts is None) & (freq is None):
            counts = np.ones(shape=len_input)
            count_dict = {}
            for s in np.unique(self.sample_id):
                idx = self.sample_id == s
                count_dict[s] = int(np.sum(counts[idx]))

            freq = []
            for c, n in zip(counts, self.sample_id):
                freq.append(c / count_dict[n])
            freq = np.asarray(freq)
            self.counts = counts
            self.freq = freq

        if hla is not None:
            if use_hla_supertype:
                hla = supertype_conv_op(hla,keep_non_supertype_alleles)
                self.use_hla_sup = True
                self.keep_non_supertype_alleles = keep_non_supertype_alleles

            self.lb_hla = MultiLabelBinarizer()
            self.hla_data_seq_num = self.lb_hla.fit_transform(hla)
            self.hla_data_seq = hla
            self.use_hla = True
        else:
            self.lb_hla = MultiLabelBinarizer()
            self.hla_data_seq_num = np.zeros([len_input,1])
            self.hla_data_seq = np.zeros(len_input)

        if Y is not None:
            if Y.ndim == 1:
                Y = np.expand_dims(Y,-1)
            self.Y = Y
            self.lb = LabelEncoder()
            self.regression = True
        else:
            self.lb = LabelEncoder()
            Y = self.lb.fit_transform(self.class_id)
            OH = OneHotEncoder(sparse=False, categories='auto')
            Y = OH.fit_transform(Y.reshape(-1, 1))
            self.Y = Y

        if w is not None:
            self.use_w = True
            self.w = w
        else:
            self.w = np.ones(len_input)

        self.seq_index = np.asarray(list(range(len(self.Y))))
        if self.regression is False:
            self.predicted = np.zeros((len(self.Y),len(self.lb.classes_)))
        else:
            self.predicted = np.zeros([len(self.Y),1])
        self.sample_list = np.unique(self.sample_id)
        print('Data Loaded')

class DeepTCR_U(DeepTCR_base):

    def _reset_models(self):
        self.models_dir = os.path.join(self.Name, 'models')
        if os.path.exists(self.models_dir):
            shutil.rmtree(self.models_dir)
        os.makedirs(self.models_dir)

    def embedding(self, latent_dim=256, kernel = 5, trainable_embedding=True, embedding_dim_aa = 64,embedding_dim_genes = 48,embedding_dim_hla=12,
                  use_only_seq=True,use_only_gene=False,use_only_hla=False,size_of_net='medium',latent_alpha=1e-3,sparsity_alpha=None,var_explained=None,graph_seed=None,
                  batch_size=10000, epochs_min=0,stop_criterion=0.01,stop_criterion_window=30, accuracy_min=None,
                  suppress_output = False,learning_rate=0.001,split_seed=None,Load_Prev_Data=False):
        #print(self.X_Seq_beta)
        if Load_Prev_Data is False:
            GO = graph_object()
            GO.size_of_net = size_of_net
            GO.embedding_dim_genes = embedding_dim_genes
            GO.embedding_dim_aa = embedding_dim_aa
            GO.embedding_dim_hla = embedding_dim_hla
            GO.l2_reg = 0.0

            graph_model_AE = tf.Graph()
            with graph_model_AE.device(self.device):
                with graph_model_AE.as_default():
                    if graph_seed is not None:
                        tf.compat.v1.set_random_seed(graph_seed)

                    GO.net = 'ae'
                    if self.use_w:
                        GO.w = tf.compat.v1.placeholder(tf.float32, shape=[None])
                    GO.Features = Conv_Model(GO, self, trainable_embedding, kernel, use_only_seq, use_only_gene, use_only_hla)

            with tf.compat.v1.Session(graph=graph_model_AE) as sess:
                sess.run(tf.compat.v1.global_variables_initializer())
                features_tensor = GO.Features

                # 运行张量，获取数据
                features_data = sess.run(features_tensor, feed_dict={GO.X_Seq_beta: self.X_Seq_beta})

            # 在 features_data 中获取 GO.Features 的数据
            print(features_data)

        return GO.Features
