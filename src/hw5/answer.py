import pdb

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from torchvision import transforms


# %%
class NN(nn.Module):
    def __init__(self, arr=[]):
        super(NN, self).__init__()
        self.relu = nn.ReLU()
        self.fc1 = nn.Linear(30 * 30 * 3, 128)
        self.fc2 = nn.Linear(128, 5)

    def forward(self, x):
        batch_size = x.shape[0]
        x = x.view(batch_size, -1)
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x


# %%
class SimpleCNN(nn.Module):
    def __init__(self, arr=[]):
        super(SimpleCNN, self).__init__()
        self.conv_layer = nn.Conv2d(3, 8, 3)
        self.pool = nn.MaxPool2d(2)
        self.fc1 = nn.Linear(1568, 5)

    def forward(self, x):
        """
        Question 2

        TODO: fill this forward function for data flow
        """
        x=self.pool(F.relu(self.conv_layer(x)))
        x=torch.flatten(x, 1) # flatten all dimensions except batch
        x = self.fc1(x)
        return x


# %%
basic_transformer = transforms.Compose([transforms.ToTensor()])

"""
Question 3

TODO: Add color normalization to the transformer. For simplicity, let us use 0.5 for mean
      and 0.5 for standard deviation for each color channel.
"""
norm_transformer = transforms.Compose([transforms.ToTensor(),transforms.Normalize(0.5,0.5)])


# %%
class DeepCNN(nn.Module):
    def __init__(self, arr=[]):
        super(DeepCNN, self).__init__()
        """
        Question 4

        TODO: setup the structure of the network
        """
        sequence=[]
        prev=3
        size=30
        for step in arr:
            if type(step) == int:
                sequence.append(nn.Conv2d(prev, step, 3))
                prev=step
                size-=2
            else:
                #It must be a pool layer
                sequence.append(nn.MaxPool2d(2))
                size=size//2
        self.layers=nn.ModuleList(sequence)

        self.fc1 = nn.Linear(4608, 5) #size * size * prev

    def forward(self, x):
        """
        Question 4

        TODO: setup the flow of data (tensor)
        """
        batch_size = x.shape[0]
        for layer in self.layers:
            if str(layer)[0:3] == "Con":
                x=F.relu(layer(x))
            else:
                x=layer(x)
        
        x=x.reshape([batch_size,4608]) # flatten all dimensions except batch
        x = self.fc1(x)
        return x


# %%
"""
Question 5

TODO:
    change the train_transformer to a tranformer with random horizontal flip
    and random affine transformation

    1. It should randomly flip the image horizontally with probability 50%
    2. It should apply random affine transformation to the image, which randomly rotate the image 
        within 5 degrees, and shear the image within 10 degrees.
    3. It should include color normalization after data augmentation. Similar to question 3.
"""

"""Add random data augmentation to the transformer"""
aug_transformer = transforms.Compose([transforms.ToTensor(),transforms.RandomHorizontalFlip(),transforms.RandomAffine(5,shear=10),transforms.Normalize(0.5,0.5)])


