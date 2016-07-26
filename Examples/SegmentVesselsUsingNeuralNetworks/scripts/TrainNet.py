#!/usr/bin/python

###########################################################################
# TrainNet.py
###########################################################################

import subprocess
import matplotlib.pyplot as plt
import time

from pylab import *

import sys
sys.path.insert(0, caffe_root + 'python') #Append pycaffe to sys.path
import caffe
from caffe import layers as L, params as P

# Path variables
caffe_root = "./"  # this file should be run from {caffe_root} (otherwise change this line)

train_net_path = caffe_root + 'data/SegmentVesselsUsingNeuralNetworks/NetProto/custom_auto_train.prototxt'
test_net_path = caffe_root + 'data/SegmentVesselsUsingNeuralNetworks/NetProto/custom_auto_test.prototxt'
solver_config_path = caffe_root + 'data/SegmentVesselsUsingNeuralNetworks/NetProto/custom_auto_solver.prototxt'

# Features square plot :
#   Plot any layer's output features with the data parameter corresponding to
#   solver.net.params['conv1'][0].data, where 'conv1' is the layer name.
def squarePlot(data):
    """Take an array of shape (n, height, width) or (n, height, width, 3)
       and visualize each (height, width) thing in a grid of size approx. sqrt(n) by sqrt(n)"""

    # normalize data for display
    data = (data - data.min()) / (data.max() - data.min())

    # force the number of filters to be square
    n = int(np.ceil(np.sqrt(data.shape[0])))
    padding = (((0, n ** 2 - data.shape[0]),
               (0, 1), (0, 1))                 # add some space between filters
               + ((0, 0),) * (data.ndim - 3))  # don't pad the last dimension (if there is one)
    data = np.pad(data, padding, mode='constant', constant_values=1)  # pad with ones (white)

    # tile the filters into an image
    data = data.reshape((n, n) + data.shape[1:]).transpose((0, 2, 1, 3) + tuple(range(4, data.ndim + 1)))
    data = data.reshape((n * data.shape[1], n * data.shape[3]) + data.shape[4:])
    filters=plt.imshow(data[:,:,0], cmap='gray', animated=True); plt.axis('off')
    filters.set_data(data[:,:,0])
    plt.draw()
    plt.pause(0.05)


# Define net architecture
def custom_net(lmdb, batch_size):
    # define your own net!
    n = caffe.NetSpec()

    # keep this data layer for all networks
    n.data, n.label = L.Data(batch_size=batch_size, backend=P.Data.LMDB, source=lmdb,
                             ntop=2, transform_param=dict(scale=1./255))

    # EDIT HERE to try different networks
    # this single layer defines a simple linear classifier
    # (in particular this defines a multiway logistic regression)
    #n.score =   L.InnerProduct(n.data, num_output=10, weight_filler=dict(type='xavier'))

    # EDIT HERE this is the LeNet variant we have already tried
    n.conv1 = L.Convolution(n.data, kernel_size=6, num_output=48, weight_filler=dict(type='xavier'))
    n.pool1 = L.Pooling(n.conv1, kernel_size=2, stride=2, pool=P.Pooling.MAX)
    n.conv2 = L.Convolution(n.pool1, kernel_size=5, num_output=48, weight_filler=dict(type='xavier'))
    n.pool2 = L.Pooling(n.conv2, kernel_size=2, stride=2, pool=P.Pooling.MAX)
    n.conv3 = L.Convolution(n.pool2, kernel_size=4, num_output=48, weight_filler=dict(type='xavier'))
    n.pool3 = L.Pooling(n.conv3, kernel_size=2, stride=2, pool=P.Pooling.MAX)
    n.conv4 = L.Convolution(n.pool3, kernel_size=2, num_output=48, weight_filler=dict(type='xavier'))
    n.pool4 = L.Pooling(n.conv4, kernel_size=2, stride=2, pool=P.Pooling.MAX)
    n.fc1 =   L.InnerProduct(n.pool4, num_output=50, weight_filler=dict(type='xavier'))
    #EDIT HERE consider L.ELU or L.Sigmoid for the nonlinearity
    #n.relu1 = L.ReLU(n.fc1, in_place=True)
    #n.drop1 = L.Dropout(n.fc1, dropout_param=dict(dropout_ratio=0.5))
    n.score = L.InnerProduct(n.fc1, num_output=2, weight_filler=dict(type='xavier'))

    # keep this loss layer for all networks
    n.loss =  L.SoftmaxWithLoss(n.score, n.label)

    return n.to_proto()

########
# Main #
########
# Create testing and training nets
with open(train_net_path, 'w') as f:
    f.write(str(custom_net('data/SegmentVesselsUsingNeuralNetworks/Net_TrainData', 128)))
with open(test_net_path, 'w') as f:
    f.write(str(custom_net('data/SegmentVesselsUsingNeuralNetworks/Net_ValData', 100)))

# Define solver
from caffe.proto import caffe_pb2
s = caffe_pb2.SolverParameter()

# Set a seed for reproducible experiments:
# this controls for randomization in training.
s.random_seed = 0xCAFFE

# Specify locations of the train and (maybe) test networks.
s.train_net = train_net_path
s.test_net.append(test_net_path)
s.test_interval = 10000  # Test after every 500 training iterations.
s.test_iter.append(1) # Test on 100 batches each time we test.

s.max_iter = 500000  # no. of times to update the net (training iterations)

# EDIT HERE to try different solvers
# solver types include "SGD", "Adam", and "Nesterov" among others.
s.type = "SGD"

# Set the initial learning rate for SGD.
s.base_lr = 1e-2  # EDIT HERE to try different learning rates
# Set momentum to accelerate learning by
# taking weighted average of current and previous updates.
s.momentum = 0.9
# Set weight decay to regularize and prevent overfitting
s.weight_decay = 5e-4

# Set `lr_policy` to define how the learning rate changes during training.
# This is the same policy as our default LeNet.
s.lr_policy = 'inv'
s.gamma = 0.0001
s.power = 0.75
# EDIT HERE to try the fixed rate (and compare with adaptive solvers)
# `fixed` is the simplest policy that keeps the learning rate constant.
s.lr_policy = 'fixed'

# Display the current training loss and accuracy every 1000 iterations.
s.display = 50000

# Snapshots are files used to store networks we've trained.
# We'll snapshot every 5K iterations -- twice during training.
s.snapshot = 20000
s.snapshot_prefix = 'data/SegmentVesselsUsingNeuralNetworks/NetProto/net'

# Train on the GPU
s.solver_mode = caffe_pb2.SolverParameter.GPU

# Write the solver to a temporary file and return its filename.
with open(solver_config_path, 'w') as f:
    f.write(str(s))

### load the solver and create train and test nets
solver = None  # ignore this workaround for lmdb data (can't instantiate two solvers on the same data)
solver = caffe.get_solver(solver_config_path)

### solve
niter = 20000  # EDIT HERE increase to train for longer
test_interval = niter / 10
# losses will also be stored in the log
train_loss = zeros(niter)
loss = 0
test_acc = zeros(int(np.ceil(niter / test_interval)))

caffe.set_mode_gpu() # Use GPU

# the main solver loop
for it in range(niter):
    ## Uncomment time lines to get the step time
    #t=time.time()
    solver.step(1)  # SGD by Caffe
    print solver.net.blobs['data'].data[50, 0][25:35,27]
    #print it," : ",(time.time()-t)
    ## store the mean loss
    loss += solver.net.blobs['loss'].data
    train_loss[it] = loss/(it+1)
    ## run a full test every so often
    ## (Caffe can also do this for us and write to a log, but we show here
    ## how to do it directly in Python, where more complicated things are easier.)
    if it % test_interval == 0:
        print 'Iteration', it, 'testing...'
        correct = 0
        for test_it in range(100):
            solver.test_nets[0].forward()
            correct += sum(solver.test_nets[0].blobs['score'].data.argmax(1)
                           == solver.test_nets[0].blobs['label'].data)
        test_acc[it // test_interval] = correct / 1e4
        print "Test acc : ", test_acc[it // test_interval]
        print "Train Loss ", it, train_loss[it]

        #Plot hidden layer features
        features=solver.net.params['conv1'][0].data
        squarePlot(features.transpose(0, 2, 3, 1))

## Show 10 first testing results
plt.imshow(solver.test_nets[0].blobs['data'].data[:10, 0].transpose(1, 0, 2).reshape(65, 10*65), cmap='gray'); axis('off')
print 'test labels:', solver.test_nets[0].blobs['label'].data[:10]
print 'score:', solver.test_nets[0].blobs['score'].data[:10]
print 'Argmax : ', solver.test_nets[0].blobs['score'].data.argmax(1)
plt.show()

## Plot output loss and test accuracy
_, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(arange(niter), train_loss)
ax2.plot(test_interval * arange(len(test_acc)), test_acc, 'r')
ax1.set_xlabel('iteration')
ax1.set_ylabel('train loss')
ax2.set_ylabel('test accuracy')
ax2.set_title('Custom Test Accuracy: {:.2f}'.format(test_acc[-1]))
plt.show()