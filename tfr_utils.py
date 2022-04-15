import tensorflow as tf
import numpy as np

#-------------------------------------------------------------------------------
# TF Records
#-------------------------------------------------------------------------------

def write_TFRecord_dataset(filepath, seqs, targets, coords=None):
  """Write x y data pairs to TF record -- with optional coordinates"""

  tf_opts = tf.io.TFRecordOptions(compression_type='ZLIB')
  with tf.io.TFRecordWriter(filepath) as writer:
    if coords:
      for seq, target, coord in zip(seqs, targets, coords):
        example_proto = serialize_example_with_coord(seq, target, coord.encode())
        writer.write(example_proto)
    else:
      for seq, target, coord in zip(seqs, targets):
        example_proto = serialize_example(seq, target)
        writer.write(example_proto)
  

def write_TFRecord_example(writer, seq, targets, coords=None):
  """Write x y data pairs to TFrecord -- with optional coordinates"""

  if coords:
    for seq, target, coord in zip(seqs, targets, coords):
      example_proto = serialize_example_with_coord(seq, target, coord.encode())
      writer.write(example_proto)
  else:
    for seq, target, coord in zip(seqs, targets):
      example_proto = serialize_example(seq, target)
      writer.write(example_proto)



def read_TFRecord_dataset(filepath, shuffle=False, batch_size=None):
  """ Load a TFRecord dataset"""
  dataset = tf.data.TFRecordDataset(filepath)
  dataset = dataset.map(read_tfrecord_example)
  dataset = dataset.prefetch(tf.data.experimental.AUTOTUNE)
  if shuffle:
    dataset = dataset.shuffle(True)
  if batch_size:
    dataset = dataset.batch(10)
  return dataset


def load_TFRecord_dataset(filepath, shuffle=False, batch_size=None, coord=None):
  """ Load a TFRecord dataset"""
  dataset = tf.data.TFRecordDataset(filepath)
  dataset = dataset.map(read_tfrecord_example)
  dataset = dataset.prefetch(tf.data.experimental.AUTOTUNE)
  inputs = []
  targets = []
  if coord:
    coords = []

  for features in dataset:
    inputs.append(features['inputs'].numpy())
    targets.append(features['targets'].numpy())
    if coord:
      coords.append(features['coord'])

  return np.array(inputs), np.array(targets), np.array(coords)


def read_tfrecord_example(serialized_example):

  if coord:
    feature_description={
      'input': tf.io.FixedLenFeature((), tf.float32),
      'target':tf.io.FixedLenFeature((), tf.float32),                
      'coord': tf.io.FixedLenFeature((), tf.str)
    }
    example = tf.io.parse_single_example(serialized_example, feature_description)
    return example['input'], example['target'], example['coord']
  else:
    feature_description={
      'input': tf.io.FixedLenFeature((), tf.float32),
      'target':tf.io.FixedLenFeature((), tf.float32)                
    }
    example = tf.io.parse_single_example(serialized_example, feature_description)
    return example['input'], example['target']




def serialize_example(seq, target):
  """
  Creates a tf.train.Example message ready to be written to a file.
  """

  feature_dict = {
      'input': feature_bytes(seq),
      'target': feature_bytes(target),
  }
  example_proto = tf.train.Example(features=tf.train.Features(feature=feature_dict))
  return example_proto.SerializeToString()

def serialize_example_with_coord(seq, target, coord):
  """
  Creates a tf.train.Example message ready to be written to a file.
  """
  
  feature_dict = {
      'coord': feature_str(coord),
      'input': feature_bytes(seq),
      'target': feature_bytes(target),
  }
  example_proto = tf.train.Example(features=tf.train.Features(feature=feature_dict))
  return example_proto.SerializeToString()


def feature_bytes(values):
  """Convert numpy arrays to bytes features."""
  
  values = values.flatten().tobytes()
  return tf.train.Feature(bytes_list=tf.train.BytesList(value=[values]))


def feature_str(values):
  """Convert str to bytes features."""
  
  return tf.train.Feature(bytes_list=tf.train.BytesList(value=[values]))

