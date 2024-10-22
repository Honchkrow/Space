from .model import Encoder_S
from tqdm import tqdm
from torch import nn
import torch.nn.functional as F
import torch
from .utils import set_seed
import numpy as np

from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import SpectralClustering

class Space():
    def __init__(self,
            matrices_dict,
            pos_similarity,
            gt,
            device='cuda:0',
            epochs=200,
            alpha=1,
            beta = 1,
            learning_rate=0.001,
            weight_decay=1e-6,
            k=20,
            seed =123,
            ):
        self.n_methods= len(matrices_dict)
        self.device = device
        self.n_spots = list(matrices_dict.values())[0].shape[0]
        self.epochs=epochs
        self.alpha = alpha
        self.beta = beta
        self.learning_rate=learning_rate
        self.weight_decay=weight_decay
        self.matrices_dict=matrices_dict
        self.k = k
        self.pos_similarity = torch.FloatTensor(pos_similarity).to(self.device)
        self.gt = gt
        self.seed =seed

    def train(self):
        set_seed(self.seed)
        self.models = Encoder_S(self.n_spots).to(self.device)
        self.optimizer = torch.optim.Adam(self.models.parameters(), self.learning_rate,
                                          weight_decay=self.weight_decay)
        self.models.train()
        for epoch in range(self.epochs):
            self.models.train()
            m_raw = self.models()
            self.matrixs = m_raw + m_raw.t()
            loss = self.loss()
            print(epoch,loss.item())

            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
        m = self.matrixs.cpu().detach().numpy()
        return (m - np.min(m)) / (np.max(m) - np.min(m))

    def loss(self):
        loss_mse=0
        for key, value in self.matrices_dict.items():
            relust = torch.FloatTensor(value).to(self.device)
            loss_temp = F.mse_loss(self.matrixs, relust, reduction='mean')
            loss_mse += loss_temp
        loss_mse= loss_mse/self.n_methods
        loss_pos = F.mse_loss(self.matrixs, self.pos_similarity, reduction='mean')

        loss_norm = self.approximate_nuclear_norm()

        return loss_mse + self.alpha * loss_pos + self.beta * loss_norm

    def approximate_nuclear_norm(self):

        M = self.matrixs.float()
        Omega = torch.randn(self.n_spots, self.k, device=M.device, dtype=M.dtype)
        Y = torch.matmul(M, Omega)
        Q, R = torch.linalg.qr(Y)
        _, s, _ = torch.svd(Q)
        approximate_nuclear_norm = s.sum()

        return approximate_nuclear_norm
