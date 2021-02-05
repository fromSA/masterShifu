import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import pickle

class GmmModel():

    def __init__(self):
        pass

    
    def fit(self, data):
        """
        PARAMETERS
        ----------
        
        data : DataFrame
            must contain the two groups, diseased and control
        
        """
        
        c = data[data.group == "control"]
        d = data[data.group == "diseased"]
        
        c_ = c[c.columns.difference(["id","group"])]
        d_ = d[d.columns.difference(["id","group"])]
        
        self.model_control = GaussianMixture(n_components=8).fit(X = c_, y = None) 
        self.model_diseased = GaussianMixture(n_components=8).fit(X = d_, y = None) 


    def save(self, filepath):
        # save the model
        with open(filepath + "/control", 'wb') as f:
            pickle.dump(self.model_control, f)
        with open(filepath + "/diseased", 'wb') as f:
            pickle.dump(self.model_diseased, f)

    def load(self, filepath):
        
        """
        PARAMETERS
        ----------
        filepath : str
            should point to a directory containing the two models of control and diseased
        
        """
        # load the model
        # set self.model
        
        
        # TODO if filepath exists.
        with open(filepath + "/control", 'rb') as f:
            self.model_control = pickle.load(f)
        with open(filepath + "/diseased", 'rb') as f:
            self.model_diseased = pickle.load(f)
        
        
    def generate_patients(self, nr_markers=12, nr_cells = 20000, nr_patients = 20, column_names=None, group=None):
        """
        PARAMETERS:
        ----------
        nr_markers : int
            nr of markers

        nr_cells : int
            number of cells per patient
        
        nr_patients : int
            number of patients to generate
        
        column_names : 
            Dataframe.columns, names for the markers
        
        group : str
            "control" or "diseased"

        RETURNS:
        -------
        patients : dataframe
            list of patients, each patient with `sample_size` cells.
        """
        
        assert (group in ["control", "diseased"]), "group must be one of [control, diseased]"
        
        
        if group == "control":
            assert (self.model_control), "load the models first"
            model = self.model_control
        else :
            assert (self.model_diseased), "load the models first"
            model = self.model_diseased
            
        if model:

            patients = np.empty(shape=(nr_patients* nr_cells, nr_markers))
            p_id = np.empty(shape=nr_patients*nr_cells, dtype="int32")

            for i in range(nr_patients):
                p_id[nr_cells*i : nr_cells*(i+1)] = np.full(shape=(nr_cells), fill_value=i+1, dtype="int32")
                patients[nr_cells*i : nr_cells*(i+1)] = model.sample(nr_cells)[0]

            patients_df = pd.DataFrame(patients, columns=column_names)
            patients_df["id"] = p_id
            patients_df["group"] = group

        return patients_df
        