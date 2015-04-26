__author__ = 'kikoc_000'

import numpy as np

"""
Class allows filtering
"""


class Kalman(object):
    isCor = True
    i = 0
    x = np.zeros((6, 100))
    xF = np.zeros((6, 100))
    z = np.zeros((6, 100))
    sigmaR = 0.1
    sigmaTH = 0.1
    T = 1
    tspine = np.arange(0, 101 * T, T)
    fi = np.pi / 4
    Rad = 2

    sigmaX = (sigmaR * np.cos(fi)) ** 2 + (Rad * sigmaTH * np.sin(fi)) ** 2
    sigmaY = (sigmaR * np.sin(fi)) ** 2 + (Rad * sigmaTH * np.cos(fi)) ** 2
    sigmaXY = (np.sin(2 * fi) * sigmaR ** 2) / 2 - (np.sin(2 * fi) * (Rad * sigmaTH) ** 2) / 2

    PF = np.array([[sigmaX, 3 * sigmaX / (2 * T), sigmaX / (T * T), sigmaXY, 3 * sigmaXY / (2 * T), sigmaXY / (T * T)],
                   [3 * sigmaX / (2 * T), 13 * sigmaX / (2 * T * T), 6 * sigmaX / (T * T * T), 3 * sigmaXY / (2 * T),
                    13 * sigmaXY / (2 * T * T), 6 * sigmaXY / (T * T * T)],
                   [sigmaX / (T * T), 6 * sigmaX / (T * T * T), 6 * sigmaX / (T * T * T * T), sigmaX / (T * T),
                    6 * sigmaX / (T * T * T), 6 * sigmaX / (T * T * T * T)],
                   [sigmaXY, 3 * sigmaXY / (2 * T), sigmaXY / (T * T), sigmaY, 3 * sigmaY / (2 * T), sigmaY / (T * T)],
                   [3 * sigmaXY / (2 * T), 13 * sigmaXY / (2 * T * T), 6 * sigmaXY / (T * T * T), 3 * sigmaY / (2 * T),
                    13 * sigmaY / (2 * T * T), 6 * sigmaY / (T * T * T)],
                   [sigmaX / (T * T), 6 * sigmaX / (T * T * T), 6 * sigmaX / (T * T * T * T), sigmaY / (T * T),
                    6 * sigmaY / (T * T * T), 6 * sigmaY / (T * T * T * T)]])

    P = np.array([[sigmaR ** 2, 3 * sigmaR ** 2 / (2 * T), sigmaR ** 2 / (T * T), 0, 0, 0],
                  [3 * sigmaR ** 2 / (2 * T), 13 * sigmaR ** 2 / (2 * T * T), 6 * sigmaR ** 2 / (T * T * T), 0, 0, 0],
                  [sigmaR ** 2 / (T * T), 6 * sigmaR ** 2 / (T * T * T), 6 * sigmaR ** 2 / (T * T * T * T), 0, 0, 0],
                  [0, 0, 0, sigmaR ** 2, 3 * sigmaR ** 2 / (2 * T), sigmaR ** 2 / (T * T)],
                  [0, 0, 0, 3 * sigmaR ** 2 / (2 * T), 13 * sigmaR ** 2 / (2 * T * T), 6 * sigmaR ** 2 / (T * T * T)],
                  [0, 0, 0, sigmaR ** 2 / (T * T), 6 * sigmaR ** 2 / (T * T * T), 6 * sigmaR ** 2 / (T * T * T * T)]])
    R = np.array([[sigmaR ** 2, 0],
                  [0, sigmaR ** 2]])
    RF = np.array([[sigmaX, sigmaXY],
                   [sigmaXY, sigmaY]])
    F = np.array([[1, T, T * T / 2, 0, 0, 0],
                  [0, 1, T, 0, 0, 0],
                  [0, 0, 1, 0, 0],
                  [0, 0, 0, 1, T, T * T / 2],
                  [0, 0, 0, 0, 1, T],
                  [0, 0, 1, 0, 0, 1]])
    H = np.array([[1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0]])
    Q = np.zeros(6)
    q = np.array([0, 0])
    I = np.array([[1, 0, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1]])
    # K = np.zeros( (6,2) )

    def __init__(self, name):
        self.name = name
        print type(self.F)
        print type(self.H)

    def filt(self):
        for i in range(100):
            if i == 1 or i == 2 or i == 3:
                pass
            else:
                # prediction
                print self.F.shape
                self.P = np.dot(np.dot(self.F, self.P), self.F) + self.Q
                self.PF = np.dot(np.dot(self.F, self.PF), self.F) + self.Q
                # correction
                k = np.dot(self.P, self.H) / (np.dot(np.dot(self.H, self.P), self.H) + self.R)
                kf = np.dot(self.PF, self.H) / (np.dot(np.dot(self.H, self.PF), self.H) + self.R)
                # x_kF = x_kF + KF*([z(1,i); z(4,i)] - H*x_kF)
                # x_k = x_k + K*([z(1,i); z(4,i)] - H*x_k)
                self.P = (self.I - k * self.H) * self.P
                self.PF = (self.I - kf * self.H) * self.PF


kalman = Kalman("KIrill")
kalman.filt()
