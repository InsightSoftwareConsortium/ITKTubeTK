#!/usr/bin/python

###########################################################################
# TrainNet.py
###########################################################################

import os
import sys
import json
import time
import shutil

import numpy as np

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import utils

# Define paths
script_params = json.load(open('params.json'))
caffe_root = str(script_params['CAFFE_SRC_ROOT'])
output_data_root = str(script_params['OUTPUT_DATA_ROOT'])
input_data_root = script_params['INPUT_DATA_ROOT']

# import caffe
sys.path.insert(0, os.path.join(caffe_root, 'python'))  # Add pycaffe
import caffe
from caffe import layers as L, params as P
from caffe.proto import caffe_pb2
import lmdb

# Define file paths
net_proto_path = os.path.join(output_data_root, 'NetProto')
snapshot_prefix = os.path.join(net_proto_path, 'net')

train_net_path = os.path.join(net_proto_path, 'net_train.prototxt')
test_net_path = os.path.join(net_proto_path, 'net_test.prototxt')
deploy_net_path = os.path.join(net_proto_path, 'net_deploy.prototxt')
solver_config_path = os.path.join(net_proto_path, 'net_solver.prototxt')

train_results_dir = os.path.join(net_proto_path, 'train_results')
utils.ensureDirectoryExists(train_results_dir)

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
    data = np.pad(data, padding, mode='constant',
                  constant_values=1)  # pad with ones (white)

    # tile the filters into an image
    data = data.reshape(
        (n, n) + data.shape[1:]).transpose((0, 2, 1, 3) + tuple(range(4, data.ndim + 1)))
    data = data.reshape(
        (n * data.shape[1], n * data.shape[3]) + data.shape[4:])

    if len(data.shape) == 3:
        data = data[:, :, 0]

    vis = plt.imshow(data, cmap='gray', animated=True)
    plt.axis('off')
    vis.set_data(data)
    plt.draw()
    plt.pause(0.05)


# Define net architecture
def custom_net(batch_size, lmdb=None):

    # define your own net!
    n = caffe.NetSpec()

    # keep this data layer for all networks
    if lmdb:

        n.data, n.label = L.Data(batch_size=batch_size, backend=P.Data.LMDB, source=lmdb,
                                 ntop=2, transform_param=dict(scale=1. / 255))

    else:

        # deploy mode
        patch_size = 2 * script_params['PATCH_RADIUS'] + 1

        n.data = L.Input(shape=[dict(dim=[batch_size, 1, patch_size, patch_size])],
                         transform_param=dict(scale=1.0 / 255.0))

    n.conv1 = L.Convolution(n.data, kernel_size=6,
                            num_output=48, weight_filler=dict(type='xavier'))
    n.relu1 = L.ReLU(n.conv1, negative_slope=0.1)
    n.pool1 = L.Pooling(n.relu1, kernel_size=2, stride=2, pool=P.Pooling.MAX)

    n.conv2 = L.Convolution(n.pool1, kernel_size=5,
                            num_output=48, weight_filler=dict(type='xavier'))
    n.relu2 = L.ReLU(n.conv2, negative_slope=0.1)
    n.pool2 = L.Pooling(n.relu2, kernel_size=2, stride=2, pool=P.Pooling.MAX)

    n.conv3 = L.Convolution(n.pool2, kernel_size=4,
                            num_output=48, weight_filler=dict(type='xavier'))
    n.relu3 = L.ReLU(n.conv3, negative_slope=0.1)
    n.pool3 = L.Pooling(n.relu3, kernel_size=2, stride=2, pool=P.Pooling.MAX)

    n.conv4 = L.Convolution(n.pool3, kernel_size=2,
                            num_output=48, weight_filler=dict(type='xavier'))
    n.relu4 = L.ReLU(n.conv4, negative_slope=0.1)
    n.pool4 = L.Pooling(n.relu4, kernel_size=2, stride=2, pool=P.Pooling.MAX)

    n.fc1 = L.InnerProduct(n.pool4, num_output=50,
                           weight_filler=dict(type='xavier'))
    n.relu5 = L.ReLU(n.fc1, negative_slope=0.1)

    n.drop1 = L.Dropout(n.relu5, dropout_param=dict(dropout_ratio=0.5))

    n.score = L.InnerProduct(n.drop1, num_output=2,
                             weight_filler=dict(type='xavier'))

    # keep this loss layer for all networks
    if lmdb:

        n.loss = L.SoftmaxWithLoss(n.score, n.label)

    else:

        # deploy mode
        n.loss = L.Softmax(n.score)

    return n.to_proto()


def custom_solver(num_train_samples, num_test_samples):

    # Define and save solver params
    solver_params = caffe_pb2.SolverParameter()

    # set solver_params type - "SGD", "Adam", and "Nesterov" etc
    solver_params.type = script_params['SOLVER_TYPE']

    # set solver params parameters
    solver_params.base_lr = script_params['BASE_LR']
    solver_params.momentum = script_params['MOMENTUM']
    solver_params.weight_decay = script_params['WEIGHT_DECAY']

    solver_params.lr_policy = script_params['LR_DECAY_POLICY']
    solver_params.gamma = script_params['GAMMA']
    solver_params.power = script_params['POWER']

    # Set a seed for reproducible experiments
    solver_params.random_seed = 0xCAFFE

    # Specify locations of the train and (maybe) test networks.
    solver_params.train_net = train_net_path
    solver_params.test_net.append(test_net_path)

    # specify train iterations
    # num passes through train dataset
    num_train_epochs = script_params['NUM_TRAIN_EPOCHS']
    train_batch_size = script_params['TRAIN_BATCH_SIZE']
    num_train_iters_per_epoch = num_train_samples / train_batch_size

    solver_params.max_iter = num_train_epochs * num_train_iters_per_epoch

    # specify test interval and iterations
    # num passes through test dataset
    num_test_epochs = script_params['NUM_TEST_EPOCHS']
    test_batch_size = script_params['TEST_BATCH_SIZE']
    num_test_iters_per_epoch = num_test_samples / test_batch_size

    # one pass through train dataset
    solver_params.test_interval = num_train_iters_per_epoch
    solver_params.test_iter.append(num_test_epochs * num_test_iters_per_epoch)

    solver_params.display = solver_params.test_interval
    solver_params.snapshot = solver_params.test_interval
    solver_params.snapshot_prefix = snapshot_prefix

    solver_params.solver_mode = caffe_pb2.SolverParameter.GPU

    return solver_params


def run():

    sys.stdout = utils.Logger(os.path.join(net_proto_path, 'net_train.log'))

    # Create testing and training net
    train_batch_size = script_params['TRAIN_BATCH_SIZE']
    train_lmdb_path = os.path.join(output_data_root, 'Net_TrainData')
    with open(train_net_path, 'w') as f:
        f.write(str(custom_net(train_batch_size, train_lmdb_path)))

    test_batch_size = script_params['TEST_BATCH_SIZE']
    test_lmdb_path = os.path.join(output_data_root, 'Net_ValData')
    with open(test_net_path, 'w') as f:
        f.write(str(custom_net(test_batch_size, test_lmdb_path)))

    with open(deploy_net_path, 'w') as f:
        f.write(str(custom_net(test_batch_size, None)))

    # get number of train/test samples
    num_train_samples = lmdb.open(
        train_lmdb_path, readonly=True).stat()['entries']
    num_test_samples = lmdb.open(
        test_lmdb_path, readonly=True).stat()['entries']

    # Define and save solver params
    solver_params = custom_solver(num_train_samples, num_test_samples)

    with open(solver_config_path, 'w') as f:
        f.write(str(solver_params))

    # Must be called before instantiating solver due to potential
    # Caffe bug
    caffe.set_mode_gpu()  # Use GPU

    # load the solver and create train and test nets
    # ignore this workaround for lmdb data (can't instantiate two solvers on
    # the same data)
    solver = None
    solver = caffe.get_solver(str(solver_config_path))

    # print the structure of the network
    print '\nNumber of training samples = ', num_train_samples
    print 'Number of testing samples = ', num_test_samples

    print('\nNetwork structure ...')
    for layer_name, blob in solver.net.blobs.iteritems():
        print(layer_name + '\t' + str(blob.data.shape))

    # print parameters at each network layer
    print('\nNetwork parameters ...')

    num_params = 0
    for layer_name, param in solver.net.params.iteritems():
        print layer_name + '\t' + str(param[0].data.shape), str(param[1].data.shape),

        cur_num_params = np.prod(param[0].data.shape) + param[1].data.shape
        print '\t%d parameters (%.2f MB)' % (cur_num_params, cur_num_params * 8.0 / 10.0**6)
        num_params += cur_num_params

    print 'Total number of params : %d (%.3f MB)' % (num_params, num_params * 8.0 / 10**6)

    # ask if user wants to start training
    flag_train = raw_input('start training (y/n)?')

    if flag_train != 'y':
        sys.exit(0)

    # solve

    fig = {}
    for layer_name in solver.net.params.keys():
        if layer_name.startswith('conv'):
            fig[layer_name] = plt.figure()

    num_train_epochs = script_params['NUM_TRAIN_EPOCHS']
    train_batch_size = script_params['TRAIN_BATCH_SIZE']
    num_train_iters_per_epoch = num_train_samples / train_batch_size
    num_train_iters = num_train_epochs * num_train_iters_per_epoch

    num_test_epochs = script_params['NUM_TEST_EPOCHS']
    test_batch_size = script_params['TEST_BATCH_SIZE']
    num_test_iters_per_epoch = num_test_samples / test_batch_size
    num_test_iters = num_test_epochs * num_test_iters_per_epoch

    assert(num_train_iters == solver_params.max_iter)
    assert(num_test_iters == solver_params.test_iter[0])
    assert(num_train_iters_per_epoch == solver_params.test_interval)

    train_loss = np.zeros(num_train_iters)
    test_acc = np.zeros(num_train_epochs)
    test_loss = np.zeros_like(test_acc)

    it_train = 0

    best_iter = -1
    best_epoch = -1
    best_acc = 0

    for ep_train in range(num_train_epochs):

        # make one pass through training data
        t = time.time()

        mean_train_loss = 0

        for i in range(num_train_iters_per_epoch):

            # make one training step
            solver.step(1)  # SGD by Caffe

            # get and store the loss
            cur_loss = solver.net.blobs['loss'].data

            train_loss[it_train] = cur_loss

            mean_train_loss += cur_loss

            it_train += 1

        epoch_train_time = time.time() - t

        mean_train_loss /= num_train_iters_per_epoch

        # validate model
        correct = 0.0
        num_test_cases = 0.0
        mean_test_loss = 0.0

        t = time.time()

        for ep_test in range(num_test_epochs):

            for i in range(num_test_iters_per_epoch):

                # make one validation step
                solver.test_nets[0].forward()

                # compute and store test accuracy
                pred_labels = solver.test_nets[0].blobs['score'].data.argmax(1)
                true_labels = solver.test_nets[0].blobs['label'].data
                correct += np.sum(pred_labels == true_labels)
                num_test_cases += len(pred_labels)

                # compute and store test loss
                mean_test_loss += solver.test_nets[0].blobs['loss'].data

        mean_test_loss /= num_test_iters
        cur_acc = correct / num_test_cases

        test_acc[ep_train] = cur_acc
        test_loss[ep_train] = mean_test_loss

        if cur_acc > best_acc:

            best_iter = it_train
            best_epoch = ep_train
            best_acc = cur_acc

        test_time = time.time() - t

        # print status
        print 'Training epoch %s / %s ...' % (ep_train + 1, num_train_epochs)

        print '\tTrain loss mean = %s' % mean_train_loss
        print '\tNo of training iterations finished = %s' % it_train
        print '\tEpoch training time = %s' % epoch_train_time

        print "\tTest loss mean = %s, accuracy = %s" % (mean_test_loss, cur_acc)
        print "\tTesting time = %s" % test_time

        # Visualize convolutional filters at each layer
        for layer_name in solver.net.params.keys():

            if layer_name.startswith('conv'):

                plt.figure(fig[layer_name].number)

                features = solver.net.params[layer_name][0].data

                n_filters, n_channels, h, w = features.shape
                n_images = n_filters * n_channels

                squarePlot(features.reshape(n_images, h, w))

                plt.title('%s - %s - Iter %s - Epoch %s' %
                          (layer_name, str(features.shape),
                           it_train, ep_train))
                plt.savefig(
                    os.path.join(
                        train_results_dir,
                        'filters_%s_ep_%.3d_it_%.7d.png' % (layer_name,
                                                            ep_train,
                                                            it_train)
                    )
                )

    print 'Test accuracy : ', test_acc
    print 'Test loss     : ', test_loss

    # save best model
    print 'Best model: accuracy = %s, iter = %s, epoch = %s' % (best_acc,
                                                                best_iter,
                                                                best_epoch)

    shutil.copyfile(
        snapshot_prefix + '_iter_%d.caffemodel' % best_iter,
        snapshot_prefix + '_best.caffemodel'
    )
    shutil.copyfile(
        snapshot_prefix + '_iter_%d.solverstate' % best_iter,
        snapshot_prefix + '_best.solverstate'
    )

    # Plot learning curve
    fig_learning_curve = plt.figure()

    ax1 = fig_learning_curve.add_subplot(111)

    train_iter_ticks = np.arange(
        num_train_iters, dtype=np.double) / num_train_iters_per_epoch
    ln1 = ax1.plot(train_iter_ticks, train_loss, 'g', label='train loss')

    ax1.plot(best_epoch, best_acc, 'k*', markersize=10)

    ax1.set_xlabel('Epochs')
    ax1.set_ylabel('train/test loss')
    ax1.set_xlim([0, num_train_epochs])

    ax2 = ax1.twinx()

    test_iter_ticks = np.arange(num_train_epochs, dtype=np.double)
    ln2 = ax1.plot(test_iter_ticks, test_loss, 'b', label='test loss')
    ln3 = ax2.plot(test_iter_ticks, test_acc, 'r', label='test accuracy')

    ax2.grid()
    ax2.set_ylabel('test accuracy')
    ax2.set_title('Learning curve : best epoch = %d, iter = %d' %
                  (best_epoch, best_iter))
    ax2.set_xlim([0, num_train_epochs])
    ax2.set_ylim([0.0, 1.0])

    lns = ln1 + ln2 + ln3
    labels = [l.get_label() for l in lns]
    ax2.legend(lns, labels, loc=0)

    fig_learning_curve.savefig(os.path.join(
        train_results_dir, 'learning_curve.png'))

    # Show sample images from each cell of the confusion matrix
    pred_labels = solver.test_nets[0].blobs['score'].data.argmax(1)
    true_labels = solver.test_nets[0].blobs['label'].data

    def generate_confusion_file(true_value, pred_value):
        adj = 'true' if true_value == pred_value else 'false'
        noun = 'positive' if pred_value == 1 else 'negative'

        c_data = solver.test_nets[0].blobs['data'].data[
            (true_labels == true_value) & (pred_labels == pred_value),
            0]
        c_count = c_data.shape[0]
        c_percent = c_count * 100. / test_batch_size
        print '%s %ss : %d / %d (%.2f%%)' % (adj, noun, c_count, test_batch_size, c_percent)
        plt.figure()
        squarePlot(c_data)
        plt.title('%s %s examples - %.2f%%' % (adj, noun, c_percent))
        plt.savefig(os.path.join(train_results_dir, 'sample_%s_%ss.png' % (adj, noun)))

    for truth in [True, False]:
        for val in [1, 0]:
            generate_confusion_file(true_value=(not truth) ^ val, pred_value=val)


if __name__ == "__main__":
    run()
