import os
import pandas as pd
import numpy as np
import numpy.ma as ma
from keras.models import load_model
from keras.preprocessing.text import Tokenizer
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import MinMaxScaler

def preprocess_seq(seq_string):

    err = None
    seq = seq_string.split()

    if len(seq) > 0 and len(seq) % 2 == 0:
        seq = [ seq[i] for i in range(1, len(seq), 2) ]
    else:
        err = 'Invalid data format!'

    # Validate input
    # Check if length = 200
    if not err:
        for i in range(len(seq)):
            if len(seq[i]) != 200:
                err = 'Incorrect sequence length!'
                break

    # Check if A, T, C, G
    if not err:
        bases = set('ACTG')
        for i in range(len(seq)):
            if set(seq[i]) != bases:
                err = 'Invalidate sequences!'
                break

    if err:
        return (None, err)
    else:
        df = pd.DataFrame()
        df['sequence'] = seq
        # X = my_pipeline.fit_transform(df)
        X = my_pipeline.transform(df)
        return (X, None)


def identify_enhancer(X):

    # layer 1 predictions
    path = [os.path.join('models-layer1', 'pcp-2x16gru1x16dense-dropout010101-wr41-weights', i)
            for i in ['model_wgts_cyc1000.h5', 'model_wgts_cyc1200.h5', 'model_wgts_cyc1400.h5',
                      'model_wgts_cyc1600.h5', 'model_wgts_cyc1800.h5']]
    model = load_model(os.path.join('models-layer1', 'pcp-2x16gru1x16dense-dropout010101-wr41.best-epch1283.h5'))
    results_df = pd.DataFrame()
    for idx, a_file in enumerate(path):
        model.load_weights(filepath=a_file, by_name=False)
        results_df['model' + str(idx)] = pd.Series(model.predict_classes(X, batch_size=None, verbose=1).flatten())
    threshold = len(results_df.columns) // 2
    results_df['layer1'] = results_df.apply(lambda x: 'enhancer' if x.sum() > threshold else 'non-enhancer', axis=1)

    # layer 2 predictions
    string_arr = np.empty((len(results_df)), dtype='U12')
    string_arr[:] = 'non-enhancer'
    results_df['layer2'] = string_arr
    print('initialisation = ', results_df['layer2'])


    print(results_df['layer1'])
    mask = (results_df['layer1'] == 'enhancer').values
    print(mask)

    # mask = np.empty(X.shape, dtype=bool)
    # mask[row_mask] = True
    # X_layer2 = ma.masked_array(X, mask)

    X_layer2 = X[np.where(mask > 0)]
    print(X_layer2)
    # X_layer2 = X

    if X_layer2.size > 0:
        print('X_layer2 is not empty!!!')

        path_layer2 = [os.path.join('models-layer2', 'layer2-2x16gru1x16dense-dropout010101-wr01-weights', i)
                       for i in ['model_wgts_cyc0400.h5', 'model_wgts_cyc0600.h5', 'model_wgts_cyc0800.h5']]
        model_layer2 = load_model(os.path.join('models-layer2',
                                               'layer2-2x16gru1x16dense-dropout010101-wr01.best-epch228.h5'))

        results2_df = pd.DataFrame()
        for idx, a_file in enumerate(path_layer2):
            model_layer2.load_weights(filepath=a_file, by_name=False)
            results2_df['model' + str(idx)] = pd.Series(model_layer2.predict_classes(X_layer2, batch_size=None, verbose=1).flatten())
        threshold2 = len(results2_df.columns) // 2
        results2_df['ensemble'] = results2_df.apply(lambda x: 'strong' if x.sum() > threshold2 else 'weak', axis=1)

        print('layer2 ensemble\n', results2_df['ensemble'])

        print('layer1\n', results_df['layer1'])

        # results_df.loc[mask, 'layer2'] = results2_df['ensemble']  # to fix !!!
        results_df.loc[results_df['layer1'] == 'enhancer', 'layer2'] = results2_df['ensemble'].values  # to fix !!!

        print('layer2\n', results_df['layer2'])

    return list(zip(results_df['layer1'], results_df['layer2']))

# This class selects the desired attributes and drops the rest.
class DataFrameSelector(BaseEstimator, TransformerMixin):

    def __init__(self, attribute_names):
        self.attribute_names = attribute_names

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return X[self.attribute_names]

# This class converts a nucleotide base (A, C, G, T) to one-hot-encoding.
class one_hot_encoder(BaseEstimator, TransformerMixin):

    def __init__(self):
        self.tokenizer = Tokenizer(num_words=4, lower=False, char_level=True)

    def fit(self, X, y=None):
        # Note that X is a data frame.
        # Fit the tokenizer on the 1st sequence in the dataset.
        self.tokenizer.fit_on_texts(X.iloc[0, 0])
        self.len_sequence = len(X.iloc[0, 0])
        return self

    def transform(self, X):
        # Note that X is a data frame.
        # one_hot_X = X.iloc[:, 0].map(lambda x: self.tokenizer.texts_to_matrix(x, mode='binary')).values
        one_hot_X = X.iloc[:, 0].map(lambda x: self.tokenizer.texts_to_matrix(x, mode='binary')).values
        one_hot_X = np.concatenate(one_hot_X)
        one_hot_X = np.reshape(one_hot_X, (-1, self.len_sequence, 4))
        return one_hot_X

# This class converts a sequence of nucleotide bases (A, C, G, T) to a sequence of dinucleotides and then to a sequence
# of pysiochemical properties of each dinucleotide.
class pcp_encoder(BaseEstimator, TransformerMixin):

    def __init__(self, pcp_df):
        self.pcp_df = pcp_df

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        # Note that X is a data frame.
        dinuc_seq = X.iloc[:, 0].map(lambda x: [ x[i:i+2] for i in range(len(x) - 1) ])
        pcp_seq = dinuc_seq.map(lambda x: [ self.pcp_df[j][i] for i in x for j in self.pcp_df.columns.tolist() ])
        # Pad with -1 for last element of sequence; it does not have an associated di-nucleotide
        pcp_seq = pcp_seq.map(
            lambda x: np.array(
                x + [-1. for i in range(len(self.pcp_df.columns))]).reshape((len(X.iloc[0, 0]),
                                                                             len(self.pcp_df.columns)))).values
        # pandas values returns a 1-D array of objects; use numpy stack to reshape it to a multi-dimensional array
        return np.stack(pcp_seq)


# This class shapes a numpy array.
class Array_Shaper(BaseEstimator, TransformerMixin):

    def __init__(self, shape):
        self.shape = shape

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return X.reshape(self.shape)

pcp_df = pd.read_csv(os.path.join('data', 'S2.csv'), index_col=0)
scaler = MinMaxScaler()
pcp_df.loc[:, :] = scaler.fit_transform(pcp_df.values)

attrbs = ['sequence']
num_bases = 4  # number of nucleotide bases
num_pcp = 6  # number of di-nucleotide physiochemical properties
# len_seq = len(all_data_df['sequence'][0])
len_seq = 200
one_hot_pipeline = Pipeline([
    ('selector', DataFrameSelector(attrbs)),
    ('one_hot_encoder', one_hot_encoder()),
    ('array_shaper2D', Array_Shaper((-1, num_bases)))
])
pcp_pipeline = Pipeline([
    ('selector', DataFrameSelector(attrbs)),
    ('pcp_encoder', pcp_encoder(pcp_df)),
    ('array_shaper2D', Array_Shaper((-1, num_pcp)))
])
union_pipeline = FeatureUnion(transformer_list=[
    ("one_hot_pipeline", one_hot_pipeline),
    ("pcp_pipeline", pcp_pipeline)
])
my_pipeline = Pipeline([
    ('feature_combiner', union_pipeline),
    ('array_shaper3D', Array_Shaper((-1, len_seq, num_bases + num_pcp)))
])

enhancer_df = pd.read_csv(os.path.join('data', 'enhancer.csv'), nrows=1)
my_pipeline.fit(enhancer_df)