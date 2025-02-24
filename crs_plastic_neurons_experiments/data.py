import numpy as np
from sklearn.model_selection import train_test_split
import torch
from torch.utils.data import TensorDataset, DataLoader

def resize_vector(vec, target_size=65598):
    if vec.shape[0] < target_size:
        return np.pad(vec, (0, target_size - len(vec)), mode='constant')
    else:
        return vec

def load_data(base_dir):
  my_listdir = os.listdir(base_dir)

  X_list = []
  Y_list = []

  for index, timestep in enumerate(my_listdir):
    path = os.path.join(base_dir, timestep)
    path = os.path.join(path, 'cMatrix.npz')
    arr = np.load(path)
    x = arr['data']
    y = arr['indices']
    x = resize_vector(x)
    y = resize_vector(y)
    X_list.append(x)
    Y_list.append(y)

  X = np.stack(X_list, axis=1)
  print(X.shape)
  Y = np.stack(Y_list, axis=1)
  print(Y.shape)

  x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

  x_train_tensor = torch.tensor(x_train, dtype=torch.float32)
  y_train_tensor = torch.tensor(y_train, dtype=torch.long)

  x_test_tensor = torch.tensor(x_test, dtype=torch.float32)
  y_test_tensor = torch.tensor(y_test, dtype=torch.long)

  train_dataset = TensorDataset(x_train_tensor, y_train_tensor)
  test_dataset = TensorDataset(x_test_tensor, y_test_tensor)

  train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
  test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

  return train_loader, test_loader