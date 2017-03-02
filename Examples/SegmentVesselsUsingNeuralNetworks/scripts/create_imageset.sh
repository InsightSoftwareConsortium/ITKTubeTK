#!/usr/bin/env sh
# Create the imagenet lmdb inputs
# N.B. set the path to the imagenet train + val data dirs

CAFFE_SRC_ROOT=`jq '.CAFFE_SRC_ROOT' params.json`
CNN_DATA_ROOT=`jq '.CNN_DATA_ROOT' params.json`
PROJECT_REL_PATH=`jq '.PROJECT_REL_PATH' params.json`

CAFFE_PROJ_ROOT="${CAFFE_SRC_ROOT}/data/${PROJECT_REL_PATH}"
CAFFE_TOOLS_DIR="${CAFFE_SRC_ROOT}/build/tools"

TRAIN_DATA_ROOT="${CNN_DATA_ROOT}/${PROJECT_REL_PATH}/training/patches/"
VAL_DATA_ROOT="${CNN_DATA_ROOT}/${PROJECT_REL_PATH}/testing/patches/"

# Set RESIZE=true to resize the images to 256x256. Leave as false if images have
# already been resized using another tool.
RESIZE=false
if $RESIZE; then
  RESIZE_HEIGHT=256
  RESIZE_WIDTH=256
else
  RESIZE_HEIGHT=0
  RESIZE_WIDTH=0
fi

if [ ! -d "$TRAIN_DATA_ROOT" ]; then
  echo "Error: TRAIN_DATA_ROOT is not a path to a directory: $TRAIN_DATA_ROOT"
  echo "Set the TRAIN_DATA_ROOT variable in create_imagenet.sh to the path" \
       "where the training patches data is stored."
  exit 1
fi

if [ ! -d "$VAL_DATA_ROOT" ]; then
  echo "Error: VAL_DATA_ROOT is not a path to a directory: $VAL_DATA_ROOT"
  echo "Set the VAL_DATA_ROOT variable in create_imagenet.sh to the path" \
       "where the validation patches data is stored."
  exit 1
fi

echo "Creating train lmdb..."

GLOG_logtostderr=1 $CAFFE_TOOLS_DIR/convert_imageset \
    --resize_height=$RESIZE_HEIGHT \
    --resize_width=$RESIZE_WIDTH \
    --shuffle \
    --gray \
    $TRAIN_DATA_ROOT \
    $TRAIN_DATA_ROOT/train.txt \
    $CAFFE_PROJ_ROOT/Net_TrainData

echo "Creating val lmdb..."

GLOG_logtostderr=1 $CAFFE_TOOLS_DIR/convert_imageset \
    --resize_height=$RESIZE_HEIGHT \
    --resize_width=$RESIZE_WIDTH \
    --shuffle \
    --gray \
    $VAL_DATA_ROOT \
    $VAL_DATA_ROOT/val.txt \
    $CAFFE_PROJ_ROOT/Net_ValData



# echo "Creating train mean ... "
# $CAFFE_TOOLS_DIR/compute_image_mean $DATA/Net_TrainData/ $DATA/Net_TrainData/mean.binaryproto
#
# echo "Creating val mean ... "
# $CAFFE_TOOLS_DIR/compute_image_mean $DATA/Net_ValData/ $DATA/Net_ValData/mean.binaryproto

echo "Done."
