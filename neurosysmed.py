import umap
import hdbscan
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor


class AC:    
    """
    This class contains the method for clustering anomalies of a masscytometry data. 
    """
    def __init__(self, embedded):
        """
        PARAMETERS
        ----------
        
        
        """
        self.df = embedded
    
    def fit(self, df, columns=None, embedd=True, train_size=.15, save_to_file=None):
        """
        PARAMETERS
        ----------
        df: DataFrame
            The entire cytometry data.
        columns: list of string
            Columns to choose from df for the entire process.
            If None, then the all columns are included.
        embedd: Bool
            choice to embedd or not.
        train_size: float between 0 and 1
            the training size for fitting umap.
        
        """
        assert(0<train_size<=1), "train size must be between 0 and 1!"
        
        if columns:
            self.df = df[columns]
            self.df["group"] = df.group
            self.df["patient"] = df.patient
        else: 
            self.df = df

        if embedd:
            raw = df[df.columns.difference(["group","patient"])]
            embedder = self.train_embedder(raw.sample(frac=train_size))
            embedding = embedder.transform(raw)
            
            self.df = pd.DataFrame(embedding, columns=["fst","snd"])
            self.df["group"] = df.group
            self.df["patient"] = df.patient
        
        if save_to_file:
            self.df.to_csv(save_to_file, index=False)

    def train_embedder(self, df, embedding_method="umap"):
        """
        PARAMETERS
        ----------
        df: DataFrame

        RETURNS
        -------
        embedder: a trained embedder
        """
        up = umap.UMAP(n_components=2, n_neighbors=10, random_state=0)
        self.trained_embedder = up.fit(df)
        return self.trained_embedder

    def collectAnomalies(self, group=None, patient_id=None, kind = "if", sample_size=50, nr_samples=20):
        """
        PARAMETERS
        ----------
        train_df: DataFrame
        
        df: DataFrame

        sample_size: int

        nr_samples: int

        embedding_method: int
        
        use_embedded: 
        
        kind : str, optional
                Type of anomaly detection algorithm to use.
                Defualt = "if", for isolation forest. "lf" for local outlier factor.
        

        RETURNS
        -------
        anoms: DataFrame 
            df with an added column "AnomCount" indecating the number of times a cells is considered an anomaly.

        """
        # select a group to collect anomlies from.
        if group:
            df = self.df[self.df.group==group]
        elif patient_id:
            df = self.df[self.df.patient==patient_id]
        else:
            df = self.df

       
        
        
        # collect the anomalies
        df = df.assign(AnomCount = 0)
        
        # take the raw data.
        raw = df[df.columns.difference(["group","patient"])]
        
        if kind == "if":
            possible_anomalies = self.find_anomalies(raw)
            anomaly_indecies = raw[possible_anomalies == -1].index
            df.loc[anomaly_indecies,"AnomCount"] = df.loc[anomaly_indecies,"AnomCount"]+1
        else:
            for sample_i in range(nr_samples):
                sample = raw.sample(sample_size)
                possible_anomalies = self.find_anomalies(sample)
                anomaly_indecies = sample[possible_anomalies == -1].index
                df.loc[anomaly_indecies,"AnomCount"] = df.loc[anomaly_indecies,"AnomCount"]+1

        return df
    
   
    def find_anomalies(self, df, kind="if"):
        """
        PARAMETERS
        ----------
        df : DataFrame 
        kind : str, optional
                Type of anomaly detection algorithm to use.
                Defualt = "if", for isolation forest. "lf" for local outlier factor.

        RETURNS
        -------
        anomalies: nparray of shape (nr_samples, )
            -1 for anomalies, 1 for outliers
        """
        assert kind in ["if","lf"], "unknown anomaly detector!"

        def iso_for(df): # Isolation forest run subsampling anyways. Do not need to subsample for it.
            clf = IsolationForest(max_samples=.3, random_state=np.random.RandomState(42))
            clf.fit(df)
            return clf.predict(df)

        def loc_fac(df):
            clf = LocalOutlierFactor(n_neighbors=2)
            return clf.fit_predict(df)

        if kind=="if":
            return iso_for(df)
        elif kind == "lf":
            return loc_fac(df)

   
    def cluster_anomalies(self, df, threshold=0, cluster_method="HDBSCAN"):
        """
        PARAMETERS
        ----------
        df:  DataFrame
        cluster_method: String
        threshold: int

        RETURNS
        -------
        clusterer: A fitted clusterer 

        """
        df = df.assign(cluster_id = np.NaN)
        anoms = df[df.AnomCount>threshold]
        anoms_raw = anoms[["fst", "snd"]]
        
        assert anoms.size >0, "There are no anomalies"
        
        clusterer = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=True)
        self.clusterer = clusterer.fit(anoms_raw)
        df.loc[anoms.index,"cluster_id"] = self.clusterer.labels_
        return df

import re
import fcsparser
from os import listdir
from os.path import isfile, join
from sklearn.preprocessing import StandardScaler

class Reader:
    """
    This class is used to read the necessary files.
    """
    def __init__(self ):
        pass
    
    def read_fcs_file(self, filepath, asinh_factor=5, asinh_transform=True, verbose=False):
        """
        PARAMETERS
        ----------
        
        RETURNS
        -------
        """
        if verbose:
            print(filepath)
        meta, data = fcsparser.parse(filepath, reformat_meta=True)
        # pick out markers that have been used and rename columns
        Pchannels = [str for str in meta.keys() if
                     (str.endswith('S') and str.startswith('$P'))]
        Pchannel_numbers = [
            [int(s) for s in re.findall(r'\d+', channel)][0]
            for channel in Pchannels]
        names = [name for name in meta['_channel_names_']]
        for k, n in enumerate(Pchannel_numbers):
            names[n-1] = meta[Pchannels[k]]
        data.columns = names
        asinh_markers = list(meta[name] for name in Pchannels)
        if asinh_transform:
            data[asinh_markers] = np.arcsinh(data[asinh_markers] / 5)
        data = data[np.sort([x for x in data.columns])]
        return data


    def get_cells_of_patient(self, patient_nr, group="pre"):
        """
        PARAMETERS
        ---------
        control_nr: int
            the index of control patients
        group: String
            'pre' for pre medicaiton group
            'post' for post medicaiton group

        RETURNS
        -------
        control_cells: DataFrame
            the cells of the patient. The number of cells varies.
        """
        assert 1 <= patient_nr <= 40

        if(patient_nr in [i for i in range(1,10)]):
            patient_nr = '0'+ str(patient_nr)
        else:
            patient_nr =str(patient_nr)

        if group == "pre" :
            # Get filename
            p_file = self.patientFileIndex[self.patientFileIndex.ID == 'Pt_'+ patient_nr].Pre_FCS_file.values[0]
        else: 
            p_file = self.patientFileIndex[self.patientFileIndex.ID == 'Pt_'+ patient_nr].Post_FCS_file.values[0]

        # Parse file
        df = self.read_fcs_file(join("DATA", p_file))
        df["patient"] = patient_nr
        return df #fcsparser.parse(join("DATA",p_file), reformat_meta=True) 


    def get_cells_of_control(self, control_nr):
        """
        PARAMETERS
        ---------
        control_nr: int
            the index of control patients

        RETURNS
        -------
        control_cells: DataFrame
            the cells of the patient. The number of cells varies.
        """
        assert 1 <= control_nr <= 10
        if(control_nr in [i for i in range(1,10)]):
            control_nr = '0'+ str(control_nr)
        else:
            control_nr =str(control_nr)

        # Get filename
        c_file = self.controlFileIndex[self.controlFileIndex.ID == 'Ctr_MS_'+ control_nr].FCS_file.values[0]

        # Parse file
        return fcsparser.parse(join("DATA",c_file), reformat_meta=True)

    def get_patients(self):
        """
        PARAMETERS
        ----------
        group: String
            'all' for all patients
        RETURNS
        -------
        data: DataFrame
            the cells of a group of patient
        """
        
        
        data = self.get_medication_group('pre')
        data["group"]="pre"
        data_2 = self.get_medication_group('post')
        data_2["group"]="post"
        data = data.append(data_2,ignore_index=True, sort=True)
        
        return data
     
    def get_medication_group(self, group="pre"):
        """
        PARAMETERS
        ----------
        group: String
            'pre' for pre medicaiton group
            'post' for post medicaiton group
        RETURNS
        -------
        data: DataFrame
            the cells of a group of patient
        """
        data = self.get_cells_of_patient(1, group)
        for i in range(2,41):
            data = data.append(self.get_cells_of_patient(i, group), sort=True)
        return data
    
    def scale(self,df):
        """
        PARAMETERS
        ----------
        df: DataFrame, all values must contain numbers.

        """
        scaler = StandardScaler()
        df[:-5] = np.arcsinh(df[:-5]/5)
        df_ = scaler.fit_transform(df)
        return pd.DataFrame(df_, columns = df.columns)

import seaborn as sns
import matplotlib.pyplot as plt
class Plotter:

    
    def __init__(self):
        pass
        
    def plot_clusters(self, df, title, just_anomalies=True, remove_legend=False):
        """
        PARAMETERS
        ----------
        """
        # choose clusters to plot
        if just_anomalies :
            df = df[df.cluster_id != np.NaN]
            title = "Cluster of anomalies" + title
        else :
            non_anos = df.cluster_id == np.NaN
            df.loc[non_anos.index, "cluster_id"] = df.cluster_id.unique().max()+1
            title = "Cluster of anomlies + non anomalies" + title

        # plot the clusters
        plt.figure(figsize= (15,15))
        sns.scatterplot(x = df.fst, y = df.snd, hue = df.cluster_id, alpha = .7,
                        palette=sns.color_palette("Set1", df.cluster_id.nunique()))
        plt.title(title)
        if remove_legend:
            plt.legend([],[], frameon=False)

    def plot_clusters_wcount(self, df, title, just_anomalies=True, remove_legend=False) :
        """
        PARAMETERS
        ----------
        
        Useless with isolation forrest is used for anomaly detection. Just run plot_clusters in that case.
        """
        # choose clusters to plot
        if just_anomalies :
            df = df[df.cluster_id != np.NaN]
            title = "Cluster of anomalies" + title
        else :
            non_anos = df.cluster_id == np.NaN
            df.loc[non_anos.index, "cluster_id"] = df.cluster_id.unique().max()+1
            title = "Cluster of anomlies + non anomalies" + title

        # plot the clusters
        plt.figure(figsize= (15,15))
        sns.scatterplot(x = df.fst, y = df.snd, hue = df.cluster_id, alpha = .7, size=df.AnomCount,
                        palette=sns.color_palette("Set1", df.cluster_id.nunique()))
        plt.title(title)
        if remove_legend:
            plt.legend([],[], frameon=False)

    def plot_heatmap_largest_clusters(self, df, whole_df, markers, clusters=3):
        """
        PARAMETERS
        ----------
        """
        # count size, ignore non anomalies
        cluster_sizes = df[df.cluster_id.notna()].groupby("cluster_id").count()

        largest_clusters = cluster_sizes.nlargest(clusters, ["fst"]).index


        # ignore the outliers during clustering.
        largest_clusters = largest_clusters[largest_clusters != -1]

        # selected clusters
        selected_clusters = df[df.cluster_id.isin(largest_clusters)]

        # the unembedded version of selected clusters, with all markers
        selected_clusters_markers = whole_df.loc[selected_clusters.index]

        # each cluster grouped 
        cluster_groups = selected_clusters_markers.groupby(selected_clusters.cluster_id).mean()

        # ordere by size
        cluster_groups =cluster_groups.reindex(largest_clusters)

        plt.figure(figsize= (len(markers)/2,clusters/2.5))
        g = sns.heatmap(cluster_groups[markers], cmap="Spectral_r", linewidths=2, cbar_kws={'label': 'mean'})
        plt.title("Relative mean of each marker in each cluster; the clusters are ordered by size")
        plt.xlabel("Markers")

    def plot_kde_cluster(self, cluster):
        """
        PARAMETERS
        ----------
        """
        f, axes = plt.subplots(15,5,figsize = (30,65), sharex=False)
        for axi, cl in zip(axes.flat, cluster.columns):
            sns.kdeplot(data=cluster[cl], ax=axi)
        f.tight_layout