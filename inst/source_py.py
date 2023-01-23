import numpy as np
import tensorly as tl
from tensorly.decomposition import parafac
from sklearn.model_selection import KFold

def parafac_svd(G, K, random_state):
  K = int(K)
  random_state = int(random_state)
  np.random.seed(random_state)
  Gparafac = parafac(tensor = G, rank = K, n_iter_max = 10000, 
  init = 'svd', random_state = random_state)[1]
  return Gparafac

def parafac_random(G, K, random_state):
  K = int(K)
  random_state = int(random_state)
  np.random.seed(random_state)
  Gparafac = parafac(tensor = G, rank = K, n_iter_max = 100, 
  init = 'random', random_state = random_state)[1]
  return Gparafac

def cutfoldid_py(n, nfolds, random_state):
  n = int(n)
  nfolds = int(nfolds)
  random_state = int(random_state)
  k_fold = KFold(nfolds, shuffle = True, random_state = random_state)
  subject_ids = np.arange(n)
  split_id_obj = k_fold.split(subject_ids)
  train_ids = []
  test_ids = []
  for train_index, test_index in split_id_obj:
      train_ids.append(train_index.copy())
      test_ids.append(test_index.copy())
  return train_ids, test_ids
