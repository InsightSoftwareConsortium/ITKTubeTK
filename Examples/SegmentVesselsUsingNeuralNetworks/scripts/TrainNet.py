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
output_data_root = str(script_params['OUTPUT_DATA_ROOT'])
input_data_root = script_params['INPUT_DATA_ROOT']

import keras
import keras.callbacks as C
import keras.layers as L
import keras.models as M
import keras.optimizers as O
import keras.utils as U

# Define file paths
net_proto_path = os.path.join(output_data_root, 'NetProto')
snapshot_format = os.path.join(net_proto_path, 'net').replace('{','{{').replace('}','}}') + \
                  '_{epoch:02d}.hdf5'

train_results_dir = os.path.join(net_proto_path, 'train_results')
utils.ensureDirectoryExists(train_results_dir)

patch_size = 2 * script_params['PATCH_RADIUS'] + 1

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


# Define net architecture, returning an uncompiled model
def create_uncompiled_model():
    # Channels go last
    inputs = L.Input(shape=(patch_size, patch_size, 1))

    # First layer set
    x = L.Conv2D(filters=48, kernel_size=6)(inputs)
    x = L.LeakyReLU(0.1)(x)
    x = L.MaxPooling2D(2)(x)

    # Second layer set
    x = L.Conv2D(filters=48, kernel_size=5)(x)
    x = L.LeakyReLU(0.1)(x)
    x = L.MaxPooling2D(2)(x)

    # Second layer set
    x = L.Conv2D(filters=48, kernel_size=4)(x)
    x = L.LeakyReLU(0.1)(x)
    x = L.MaxPooling2D(2)(x)

    # Second layer set
    x = L.Conv2D(filters=48, kernel_size=2)(x)
    x = L.LeakyReLU(0.1)(x)
    x = L.MaxPooling2D(2)(x)

    # Fully connected layer set
    x = L.Flatten()(x)
    x = L.Dense(50)(x)
    x = L.LeakyReLU(0.1)(x)

    x = L.Dropout(0.5)(x)

    # Classify
    x = L.Dense(2)(x)
    predictions = L.Activation('softmax')(x)

    return M.Model(inputs=inputs, outputs=predictions)


# Configure and compile model
def compile_model(model):

    model.compile(O.SGD(lr=script_params['BASE_LR'],
                        momentum=script_params['MOMENTUM'],
                        decay=script_params['GAMMA']),
                  'categorical_crossentropy',
                  metrics=['accuracy'])


def queryResultToModelArguments(result):
    """Convert the list of (image, label) pairs from a query into a pair
    of NumPy arrays to pass into various model functions.

    """
    image_data = (np.array([np.frombuffer(im, dtype=np.uint8) for im, _ in result])
                  .reshape((len(result), patch_size, patch_size, 1)) / 255.)
    labels = np.array([l for _, l in result])
    return image_data, labels


def run():

    sys.stdout = utils.Logger(os.path.join(net_proto_path, 'net_train.log'))

    # Create testing and training net
    train_batch_size = script_params['TRAIN_BATCH_SIZE']
    train_db_path = os.path.join(output_data_root, 'Net_TrainData')

    test_batch_size = script_params['TEST_BATCH_SIZE']
    test_db_path = os.path.join(output_data_root, 'Net_ValData')

    model = create_uncompiled_model()
    compile_model(model)

    def get_sample_count(path):
        db = utils.open_sqlite3_db(path)
        retval, = db.execute('''select count(*) from "Patches"''').fetchone()
        db.close()
        return retval

    # get number of train/test samples
    num_train_samples = get_sample_count(train_db_path)
    num_test_samples = get_sample_count(test_db_path)

    # print the structure of the network
    print '\nNumber of training samples = ', num_train_samples
    print 'Number of testing samples = ', num_test_samples

    # TODO translate
    '''
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
    '''

    # ask if user wants to start training
    flag_train = raw_input('start training (y/n)?')

    if flag_train != 'y':
        sys.exit(0)

    # solve
    # TODO figure out how to translate
    '''
    fig = {}
    for layer_name in solver.net.params.keys():
        if layer_name.startswith('conv'):
            fig[layer_name] = plt.figure()
    '''

    num_train_epochs = script_params['NUM_TRAIN_EPOCHS']
    train_batch_size = script_params['TRAIN_BATCH_SIZE']
    num_train_iters_per_epoch = num_train_samples / train_batch_size
    num_train_iters = num_train_epochs * num_train_iters_per_epoch

    num_test_epochs = script_params['NUM_TEST_EPOCHS']
    test_batch_size = script_params['TEST_BATCH_SIZE']
    num_test_iters_per_epoch = num_test_samples / test_batch_size
    num_test_iters = num_test_epochs * num_test_iters_per_epoch

    def data_generator(db_path, batch_size):
        """Generate batches of batch_size samples from the database at
        db_path, looping through the database repeatedly, and
        ignoring remaining samples when the total sample count is
        divided by batch_size.

        """
        db = utils.open_sqlite3_db(db_path)
        cursor = db.cursor()
        while True:
            cursor.execute('''select "image_data", "patch_index"
                              from "Patches"''')
            while True:
                result = cursor.fetchmany(batch_size)
                if len(result) < batch_size:
                    break
                image_data, labels = queryResultToModelArguments(result)
                yield image_data, U.to_categorical(labels, 2)

    history = model.fit_generator(data_generator(train_db_path, train_batch_size),
                                  steps_per_epoch=num_train_iters_per_epoch,
                                  epochs=num_train_epochs,
                                  callbacks=[C.ProgbarLogger('steps'),
                                             C.ModelCheckpoint(snapshot_format)],
                                  validation_data=data_generator(test_db_path, test_batch_size),
                                  validation_steps=num_test_iters_per_epoch).history

    print 'Test accuracy : ', history['val_acc']
    print 'Test loss     : ', history['val_loss']

    best_epoch = np.argmax(history['val_acc'])
    best_acc = history['val_acc'][best_epoch]

    # save best model
    print 'Best model: accuracy = %s, epoch = %s' % (best_acc, best_epoch)

    shutil.copyfile(
        snapshot_format.format(epoch=best_epoch),
        os.path.join(net_proto_path, 'net_best.hdf5')
    )

    '''
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
    '''

    # Plot learning curve
    fig_learning_curve = plt.figure()

    ax1 = fig_learning_curve.add_subplot(111)

    ln1 = ax1.plot(history['loss'], 'g', label='train loss')

    ax1.plot(best_epoch, best_acc, 'k*', markersize=10)

    ax1.set_xlabel('Epochs')
    ax1.set_ylabel('train/test loss')
    ax1.set_xlim([0, num_train_epochs])

    ax2 = ax1.twinx()

    ln2 = ax1.plot(history['val_loss'], 'b', label='test loss')
    ln3 = ax2.plot(history['val_acc'], 'r', label='test accuracy')

    ax2.grid()
    ax2.set_ylabel('test accuracy')
    ax2.set_title('Learning curve : best epoch = %d' % best_epoch)
    ax2.set_xlim([0, num_train_epochs])
    ax2.set_ylim([0.0, 1.0])

    lns = ln1 + ln2 + ln3
    labels = [l.get_label() for l in lns]
    ax2.legend(lns, labels, loc=0)

    fig_learning_curve.savefig(os.path.join(
        train_results_dir, 'learning_curve.png'))

    # Show sample images from each cell of the confusion matrix
    image_data, true_labels = queryResultToModelArguments(
        utils.open_sqlite3_db(test_db_path)
        .execute('''select "image_data", "patch_index" from "Patches"
                    order by random() limit ?''', (test_batch_size,))
        .fetchall())
    pred_labels = model.predict(image_data, test_batch_size).argmax(1)

    def generate_confusion_file(true_value, pred_value):
        adj = 'true' if true_value == pred_value else 'false'
        noun = 'positive' if pred_value == 1 else 'negative'

        c_data = image_data[(true_labels == true_value) & (pred_labels == pred_value), ..., 0]
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
