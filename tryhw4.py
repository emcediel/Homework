import sys
import pickle
sys.dont_write_bytecode = True
from uwimg import *

def softmax_model(inputs, outputs):
    l = [make_layer(inputs, outputs, SOFTMAX)]
    return make_model(l)

def neural_net(inputs, outputs):
    print(inputs)
    l = [   make_layer(inputs, 64, LRELU), #LOGISTIC
            make_layer(64,32,   LRELU),
            make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

print("loading data...")
train = load_classification_data(b"cifar.train", b"cifar/labels.txt", 1)#load_classification_data(b"mnist.train", b"mnist.labels", 1)
test  = load_classification_data(b"cifar.test", b"cifar/labels.txt", 1)#load_classification_data(b"mnist.test", b"mnist.labels", 1)
print("done")
print

print("training model...")
batch = 128
iters = 10000
rate = .05 #.01
momentum = .9
decay =pow(5,-4) #.0 #next do pow(5,-4)



#m = neural_net(train.X.cols, train.y.cols) #m = softmax_model(train.X.cols, train.y.cols) uncomment later
#train_model(m, train, batch, iters, rate, momentum, decay) uncomment later
#print("done")
#print
#print("evaluating model...")
#print("training accuracy: %0.2f %%"%(100*accuracy_model(m, train)))
#print("test accuracy:     %0.2f %%"%(100*accuracy_model(m, test)))
itera=[10000]
learning_rate=[pow(5,-2),pow(5,-3)]
weight_decay=[pow(10,-6),pow(10,-7)]#pow(10,-1),pow(10,-2),pow(10,-3),pow(10,-4),
res_train=[]
res_test=[]
for lr in learning_rate:
    m = neural_net(train.X.cols, train.y.cols)
    train_model(m, train, batch, iters, lr, momentum, decay) 
    res_train.append(100*accuracy_model(m, train))
    res_test.append(100*accuracy_model(m, test))

for i in range(len(res_train)):
    print("iter: %d %% training accuracy: %0.2f %%"%(learning_rate[i],res_train[i]))
    print("iter: %d %% test accuracy:     %0.2f %%"%(learning_rate[i],res_test[i]))
#m = softmax_model(train.X.cols, train.y.cols)
#train_model(m, train, batch, iters, rate, momentum, w_decay)
#print("evaluating model...")
#print("decay: %0.2f %% training accuracy: %0.2f %%"%(w_decay,100*accuracy_model(m, train)))
#print("decay: %0.2f %% test accuracy:     %0.2f %%"%(w_decay,100*accuracy_model(m, test)))