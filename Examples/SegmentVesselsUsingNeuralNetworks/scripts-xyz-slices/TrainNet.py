#!/usr/bin/python

###########################################################################
# TrainNet.py
###########################################################################

import itertools
import os
import sys
import time
import shutil

import numpy as np

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import utils
from utils import script_params

# Define paths
output_data_root = str(script_params['OUTPUT_DATA_ROOT'])
input_data_root = script_params['INPUT_DATA_ROOT']

import keras
import keras.callbacks as C
import keras.layers as L
import keras.models as M
import keras.optimizers as O
import keras.regularizers as R
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
    weight_decay = script_params['WEIGHT_DECAY']

    def wrap_regularizer(layer):
        def wrapper(*args, **kwargs):
            return layer(*args,
                         kernel_regularizer=R.l2(weight_decay),
                         bias_regularizer=R.l2(weight_decay),
                         **kwargs)
        return wrapper

    design = script_params['NETWORK_DESIGN']

    if design == 'xyz':
        L_ConvND = L.Conv2D
        MaxPoolingND = L.MaxPooling2D
        ndim = 2
        ninputs = 3
    elif design == 'full3d':
        L_ConvND = L.Conv3D
        MaxPoolingND = L.MaxPooling3D
        ndim = 3
        ninputs = 1
    else:
        raise ValueError('Unknown NETWORK_DESIGN')

    ConvND = wrap_regularizer(L_ConvND)
    Dense = wrap_regularizer(L.Dense)
    def BatchNormalization():
        return L.BatchNormalization(scale=False, beta_regularizer=R.l2(weight_decay))

    def convLayer(f=32, k=3):
        c = ConvND(filters=f, kernel_size=k, use_bias=False)
        n = BatchNormalization()
        r = L.LeakyReLU(0.1)
        return lambda x: r(n(c(x)))

    def denseLayer(u):
        d = Dense(u, use_bias=False)
        n = BatchNormalization()
        r = L.LeakyReLU(0.1)
        return lambda x: r(n(d(x)))

    # Channels go last
    input_shape = (patch_size,) * ndim + (1,)

    inputs = [L.Input(shape=input_shape) for _ in range(ninputs)]

    sharedInput = L.Input(shape=input_shape)

    out_size = patch_size
    x = sharedInput

    # Convolutional layer sets
    while out_size >= 4:
        x = convLayer()(x)
        x = MaxPoolingND(2)(x)
        out_size = int((out_size - 2) / 2)

    # Fully connected layer set
    x = L.Flatten()(x)
    x = denseLayer(32)(x)

    x = L.Dropout(0.5)(x)

    sharedModel = M.Model(inputs=sharedInput, outputs=x)

    x = L.Concatenate()([sharedModel(i) for i in inputs])

    x = denseLayer(20)(x)

    # Classify
    x = Dense(2)(x)
    # Use gamma because softmax is scale-sensitive
    x = L.BatchNormalization(beta_regularizer=R.l2(weight_decay),
                             gamma_regularizer=R.l2(weight_decay))(x)
    predictions = L.Activation('softmax')(x)

    return M.Model(inputs=inputs, outputs=predictions)


# Configure and compile model
def compile_model(model):
    optimizer = getattr(O, script_params['SOLVER_TYPE'])
    kwargs = {k.lower(): v for k, v in script_params['SOLVER_PARAMS'].items()}

    model.compile(optimizer(lr=script_params['BASE_LR'],
                            decay=script_params['GAMMA'],
                            **kwargs),
                  'categorical_crossentropy',
                  metrics=['accuracy'])


def queryResultToModelArguments(result):
    """Convert the list of (image, label) pairs from a query into a pair
    of a list of inputs and a label array to pass into various model
    functions.

    """
    design = script_params['NETWORK_DESIGN']
    if design == 'xyz':
        patch_shape = patch_size, patch_size, 3
    elif design == 'full3d':
        patch_shape = (patch_size,) * 3
    else:
        raise ValueError('Unknown NETWORK_DESIGN')
    image_data = np.stack(
        np.frombuffer(im, dtype=np.float16).reshape(patch_shape)
        for im, _ in result
    )
    labels = np.array([l for _, l in result])
    return utils.prepareInputArray(image_data), labels


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

    def with_estimated_mem(num):
        return "{} ({:.3f} MB)".format(num, num*8/1e6)

    def flatten_model(m):
        if isinstance(m, M.Model):
            return itertools.chain.from_iterable(flatten_model(l) for l in m.layers)
        else:
            return [m]

    layers = list(flatten_model(model))

    # Warning: this summation assumes the memory used by a layer for
    # its tensors is determined only by its output.  This estimate
    # will be an overestimate to the extent that any operation can be
    # done "in-place", and an underestimate to the extent a layer is
    # more complicated than a single step of input to output.
    total_elements = 0
    print '\nNetwork structure ...'
    for layer in layers:
        print layer.name + '\t' + str(layer.output_shape),
        elements = train_batch_size * np.prod(layer.output_shape[1:])
        print '\tElements:', with_estimated_mem(elements)
        total_elements += elements
    print 'Total elements:', with_estimated_mem(total_elements)

    print '\nNetwork weights ...'

    total_params = 0
    for layer in layers:
        shapes = [w.shape for w in layer.get_weights()]
        if not shapes:
            continue
        print layer.name + '\t' + ' '.join(map(str, shapes)),
        num_params = sum(map(np.prod, shapes))
        print '\tParameters:', with_estimated_mem(num_params)
        total_params += num_params
    print 'Total parameters:', with_estimated_mem(total_params)

    # ask if user wants to start training
    flag_train = raw_input('start training (y/n)?')

    if flag_train != 'y':
        sys.exit(0)

    # solve
    conv_layers = [l for l in layers if isinstance(l, L.Conv2D)]
    # Auto-generated names are supposed to be unique
    fig = {l.name: plt.figure() for l in conv_layers}

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
        utils.best_model_path(),
    )

    design = script_params['NETWORK_DESIGN']

    # Visualize convolutional filters at each layer
    ep_train = num_train_epochs - 1
    for l in conv_layers:

        plt.figure(fig[l.name].number)

        features = l.get_weights()[0]

        if design == 'xyz':
            h, w, n_channels, n_filters = features.shape
            n_images = n_filters * n_channels

            squarePlot(features.transpose(2, 3, 0, 1).reshape(n_images, h, w))
        elif design == 'full3d':
            # TODO how to visualize 3D data?
            pass
        else:
            raise ValueError('Unknown NETWORK_DESIGN')

        plt.title('%s - %s - Epoch %s' %
                  (l.name, str(features.shape), ep_train))
        plt.savefig(
            os.path.join(
                train_results_dir,
                'filters_%s_ep_%.3d.png' % (l.name, ep_train)
            )
        )

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

    if design == 'xyz':
        # TODO this generates color plots.  Not helpful.
        image_data = np.stack(image_data, -2)
    elif design == 'full3d':
        pass
    else:
        raise ValueError('Unknown NETWORK_DESIGN')
    image_data.shape = image_data.shape[:-1]

    def generate_confusion_file(true_value, pred_value):
        adj = 'true' if true_value == pred_value else 'false'
        noun = 'positive' if pred_value == 1 else 'negative'

        c_data = image_data[(true_labels == true_value) & (pred_labels == pred_value)]
        c_count = c_data.shape[0]
        c_percent = c_count * 100. / test_batch_size
        print '%s %ss : %d / %d (%.2f%%)' % (adj, noun, c_count, test_batch_size, c_percent)
        plt.figure()
        if design == 'xyz':
            squarePlot(c_data)
        elif design == 'full3d':
            # TODO how to visualize 3D data?
            pass
        else:
            raise ValueError('Unknown NETWORK_DESIGN')
        plt.title('%s %s examples - %.2f%%' % (adj, noun, c_percent))
        plt.savefig(os.path.join(train_results_dir, 'sample_%s_%ss.png' % (adj, noun)))

    for truth in [True, False]:
        for val in [1, 0]:
            generate_confusion_file(true_value=(not truth) ^ val, pred_value=val)


if __name__ == "__main__":
    t = time.time()
    try:
        run()
    finally:
        print("Running time: " + str(time.time() - t))
