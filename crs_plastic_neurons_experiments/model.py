import torch
import torch.nn as nn

class ThreeLayerNet(nn.Module):
    def __init__(self):
        super(ThreeLayerNet, self).__init__()
        self.fc1 = nn.Linear(257, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 257)
    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x