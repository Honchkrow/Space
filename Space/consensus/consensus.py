from .model import Encoder_S
from tqdm import tqdm
from torch import nn
import torch.nn.functional as F
import torch
from .utils import set_seed
import numpy as np

from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import SpectralClustering

class consensus():
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
        # self.optimizer = torch.optim.SGD(self.models.parameters(), lr=self.learning_rate, momentum=0.9)
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
            if (epoch+1)%10==0:
            # if epoch > 11:
                sc = SpectralClustering(n_clusters=self.k, affinity='precomputed',n_jobs=128,random_state=self.seed)
                m=self.matrixs.cpu().detach().numpy()
                labels = sc.fit_predict(m)
                ari = adjusted_rand_score(labels, self.gt)
                print("ari:",ari)

        return self.matrixs.cpu().detach().numpy()

    def loss(self):
        loss_mse=0
        for key, value in self.matrices_dict.items():
            relust = torch.FloatTensor(value).to(self.device)
            loss_temp = F.mse_loss(self.matrixs, relust, reduction='mean')
            loss_mse += loss_temp
        loss_mse= loss_mse/self.n_methods
        loss_pos = F.mse_loss(self.matrixs, self.pos_similarity, reduction='mean')

        loss_norm = self.approximate_nuclear_norm()

        # print(loss_mse.item(),loss_pos.item(),loss_norm.item())

        return loss_mse + self.alpha * loss_pos + self.beta * loss_norm

    def approximate_nuclear_norm(self):
        """
        使用PyTorch估计矩阵M的核范数（奇异值之和）的近似值，
        包含一个QR分解步骤以改善数值稳定性。

        参数:
        - M: 原始矩阵，PyTorch张量。
        - k: 随机投影的维度，应小于M的列数。

        返回:
        - 核范数的近似值。
        """
        # 确保M是浮点张量，以便进行SVD
        M = self.matrixs.float()

        # 生成随机矩阵Omega。这个随机矩阵用于将原始矩阵M投影到一个较低的维度空间，减少计算量。
        Omega = torch.randn(self.n_spots, self.k, device=M.device, dtype=M.dtype)

        # 计算投影矩阵Y = MOmega。这一步实质上是在用随机矩阵Omega压缩M，减少其列数，从而降维。
        Y = torch.matmul(M, Omega)

        # 对Y进行QR分解。QR分解可以将矩阵分解为一个正交矩阵Q和一个上三角矩阵R。
        # 这一步是为了改善数值稳定性，通过预先正交化来减少后续SVD的数值误差。
        Q, R = torch.linalg.qr(Y)

        # 对Q进行SVD，只计算奇异值。由于Y=QR，而Q是正交的，因此Y的奇异值与R的奇异值相同。
        # 这里直接对Q进行SVD是为了利用Y的数值稳定性改进。
        _, s, _ = torch.svd(Q)

        # 计算核范数的近似值（奇异值之和）。这个值是原始矩阵M通过随机投影和数值稳定化处理后的奇异值之和，用作核范数的近似。
        approximate_nuclear_norm = s.sum()

        return approximate_nuclear_norm
