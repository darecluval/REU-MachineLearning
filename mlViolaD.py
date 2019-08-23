# Machine Learning - Viola D - D Major
# Clare DuVal - MSU/Clemson
# Updated: August 23, 2019

#IMPORTS
import numpy as np
import tensorflow as tf
from tensorflow.python.training import checkpoint_state_pb2
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import random
import uuid
import datetime
from PIL import Image
from sklearn.model_selection import train_test_split
import helper as helper

tf.logging.set_verbosity(tf.logging.INFO)

#NETWORK PARAMETERS
epochs = 10
learning_rate = 0.01
batch_size = 156
batch_size2 = batch_size / 2
batch_size4 = batch_size / 4
print("Number of Epochs: " + str(epochs))
print("Learning Rate: " + str(learning_rate))
print("Batch Size: " + str(batch_size))

#Variables
train_num = 56
test_num = 70 - train_num
num_nets = 9
num_classes = 13
dim = 200 #height and width

#-------Input Train and Test Data-------#

#Pre-Allocating Images
train_images = np.zeros((num_sets * train_num * num_classes, dim*dim))
train_labels = np.zeros((num_sets * train_num * num_classes, num_classes))
test_images = np.zeros((num_sets * test_num * num_classes, dim*dim))
test_labels = np.zeros((num_sets * test_num * num_classes, num_classes))

#Randomly shuffles notes to be collected
notes = list(range(1,71))
np.random.shuffle(notes)
#print(notes)

#Change working directory
PATH = "/Users/clareduval/Documents/Python/Spectrograms"
os.chdir(PATH)

#Files to be read are in format "set_note_j.jpeg"
#Training Data
count = 0
for j in range(train_num):
    for set in range(num_sets):
        for note in range(num_classes):
            filename = str(set) + '_' + str(note+1) + '_' + str(notes[j]) + '.jpeg'
            temp_image = Image.open(filename)
            train_images[count,:] = np.reshape(temp_image, -1)
            train_labels[count, note] = 1
            count += 1

train_images = np.array(train_images, dtype = np.uint8)

#Test Data
count = 0
for k in range(test_num):
    for set in range(num_sets):
        for note in range(num_classes):
            filename = str(set) + '_' + str(note+1) + '_' + str(notes[j]) + '.jpeg'
            temp_image = Image.open(filename)
            test_images[count,:] = np.reshape(temp_image, -1)
            test_labels[count, note] = 1
            count += 1

#Shuffles images and labels simultaneously
train_random = list(zip(train_images, train_labels))
random.shuffle(train_random)
train_images, train_labels = zip(*train_random)

test_random = list(zip(test_images, test_labels))
random.shuffle(test_random)
test_images, test_labels = zip(*test_random)

#Label Dictionary
labels = {
    1: 'Detache Down',
    2: 'Detache Up',
    3: 'Slur 2',
    4: 'Slur 4',
    5: 'Slur 8',
    6: 'Staccato',
    7: 'Staccato Hooked',
    8: 'Portato Hooked',
    9: 'Portato Separate',
    10: 'Spiccato',
    11: 'Spiccato Fast',
    12: 'Ricochet Fast',
    13: 'Ricochet'
}

#Convert Images to Tensors
train = np.asarray(train_images, np.float32)
test = np.asarray(test_images, np.float32)
train_x = train.reshape(-1, dim, dim, 1)
test_x = test.reshape(-1, dim, dim, 1)

print("Images imported successfully!")

#Placeholder
x = tf.placeholder("float", [None, dim, dim, 1])
y = tf.placeholder("float", [None, num_classes])

#Weights
weights = {
    'wc1': tf.get_variable('W0', shape = (3, 3,  1,  batch_size), initializer = tf.contrib.layers.xavier_initializer()),
    'wc2': tf.get_variable('W1', shape = (3, 3, batch_size,  batch_size2), initializer = tf.contrib.layers.xavier_initializer()),
    'wc3': tf.get_variable('W2', shape = (3, 3, batch_size2, batch_size4), initializer = tf.contrib.layers.xavier_initializer()),
    'wd1': tf.get_variable('W3', shape = (5*5*batch_size4, batch_size4), initializer = tf.contrib.layers.xavier_initializer()),
    'out': tf.get_variable('W6', shape = (batch_size4, num_classes), initializer = tf.contrib.layers.xavier_initializer())
}

#Biases
biases = {
    'bc1': tf.get_variable('B0', shape = (batch_size), initializer = tf.contrib.layers.xavier_initializer()),
    'bc2': tf.get_variable('B1', shape = (batch_size2), initializer = tf.contrib.layers.xavier_initializer()),
    'bc3': tf.get_variable('B2', shape = (batch_size4), initializer = tf.contrib.layers.xavier_initializer()),
    'bd1': tf.get_variable('B3', shape = (batch_size4), initializer = tf.contrib.layers.xavier_initializer()),
    'out': tf.get_variable('B4', shape = (num_classes), initializer = tf.contrib.layers.xavier_initializer())
}

#Loss and Optimizer Nodes
pred = helper.conv_net(x, weights, biases)
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits = pred, labels = y))
optimizer = tf.train.AdamOptimizer(learning_rate = learning_rate).minimize(cost)

#Evaluate Model Node
correct_prediction = tf.equal(tf.argmax(pred,1), tf.argmax(y, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

#Initialize the variables
init = tf.global_variables_initializer()

#---TRAIN & TEST THE MODEL---#
i = 1
with tf.Session() as sess:

    sess.run(init)
    train_loss = []
    test_loss = []
    train_accuracy = []
    test_accuracy = []
    summary_writer = tf.summary.FileWriter('./Output', sess.graph)

    for ep in range(epochs):
        for batch in range(len(train_x)//batch_size):
            print("Batch " + str(batch+1))
            batch_x = train_x[batch*batch_size:min((batch+1)*batch_size, len(train_x))]
            batch_y = train_labels[batch*batch_size:min((batch+1)*batch_size, len(train_labels))]

            #Calculate batch loss and accuracy
            opt = sess.run(optimizer, feed_dict = {x: batch_x, y: batch_y})
            loss, acc = sess.run([cost, accuracy], feed_dict = {x: batch_x, y: batch_y})

        #Calculate loss and training accuracy
        print("Epoch: " + str(ep+1) + " out of " + str(epochs))
        print("Loss = " + \
                        "{:.6f}".format(loss) + ", Training Accuracy = " + \
                        "{:.5f}".format(acc))
        print("Optimization finished!")

        #Calculate and append testing accuracy
        print("Calculating testing accuracy...")
        test_acc, valid_loss = sess.run([accuracy, cost], feed_dict = {x: test_X, y: test_labels})
        train_loss.append(loss)
        test_loss.append(valid_loss)
        train_accuracy.append(acc)
        test_accuracy.append(test_acc)

        print("Testing Accuracy: ", "{:.5f}".format(test_acc))
        learning_rate = learning_rate * 0.85

    summary_writer.close()
