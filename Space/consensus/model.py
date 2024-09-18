import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module

class Encoder_S(torch.nn.Module):
    def __init__(self,  n_spot):
        super(Encoder_S, self).__init__()
        self.n_spot = n_spot
        self.M = Parameter(torch.FloatTensor(self.n_spot, self.n_spot))
        self.reset_parameters()
    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.M)
    def forward(self):
        x = self.M
        return x
