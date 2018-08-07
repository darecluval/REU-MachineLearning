#Machine Learning - Viola D - D Major
#Clare DuVal - MSU/Clemson
#Updated: July 24, 2018

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

tf.logging.set_verbosity(tf.logging.INFO)

#NETWORK PARAMETERS
epochs = 10
learning_Rate = 0.01
batch_Size = 156
batch_Size2 = batch_Size / 2
batch_Size4 = batch_Size / 4
print("Number of Epochs: " + str(epochs))
print("Learning Rate: " + str(learning_Rate))
print("Batch Size: " + str(batch_Size))

#Variables
train_Num = 56
test_Num = 70 - train_Num
num_Sets = 9
num_Classes = 13
dim = 200 #height and width

#-------Input Train and Test Data-------#

#Pre-Allocating Images
train_Images = np.zeros((num_Sets * train_Num * num_Classes, dim*dim))
train_Labels = np.zeros((num_Sets * train_Num * num_Classes, num_Classes))
test_Images = np.zeros((num_Sets * test_Num * num_Classes, dim*dim))
test_Labels = np.zeros((num_Sets * test_Num * num_Classes, num_Classes))

#Randomly shuffles notes to be collected
notes = list(range(1,71))
np.random.shuffle(notes)
#print(notes)

#Change working directory
PATH = "/Users/clareduval/Documents/Python/Spectrograms"
os.chdir(PATH)

#Files to be read are in format "set_note_j.jpeg"

#print("Train Images:")

#Training Data
count = 0
j = 0
while j < train_Num:
    #num_Sets
    set = 0
    while set < num_Sets:
        #13 styles
        note = 1
        while note <= 13:
            FILENAME = str(set) + '_' + str(note) + '_' + str(notes[j]) + '.jpeg'
            temp_Image = Image.open(FILENAME)
            train_Images[count,:] = np.reshape(temp_Image, -1)
            train_Labels[count, note-1] = 1
            #print(str(count+1))
            note += 1
            count += 1
        set += 1
    j += 1
    #print(j)

train_Images = np.array(train_Images, dtype = np.uint8)

#Test Data
count = 0
k = 0
while k < test_Num:
    #num_Sets
    set = 0
    while set < num_Sets:
        #13 styles
        note = 1
        while note <= 13:
            FILENAME = str(set) + '_' + str(note) + '_' + str(notes[j]) + '.jpeg'
            temp_Image = Image.open(FILENAME)
            test_Images[count,:] = np.reshape(temp_Image, -1)
            test_Labels[count, note-1] = 1
            #print(str(count+1))
            note += 1
            count += 1
        set += 1
    j += 1
    k += 1
    #print(j)

#Shuffles images and labels simultaneously
trainRANDOM = list(zip(train_Images, train_Labels))
random.shuffle(trainRANDOM)
train_Images, train_Labels = zip(*trainRANDOM)
testRANDOM = list(zip(test_Images, test_Labels))
random.shuffle(testRANDOM)
test_Images, test_Labels = zip(*testRANDOM)

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
#print(labels)

#Convert Images to Tensors
train = np.asarray(train_Images, np.float32)
test = np.asarray(test_Images, np.float32)
train_X = train.reshape(-1, dim, dim, 1)
test_X = test.reshape(-1, dim, dim, 1)

print("Images imported successfully!")

#Placeholder
x = tf.placeholder("float", [None, dim, dim, 1])
y = tf.placeholder("float", [None, num_Classes])

#Define Convolutional layer
def conv2d(x, W, b, strides = 1):
    x = tf.nn.conv2d(x, W, strides = [1, strides, strides, 1], padding = "SAME")
    x = tf.nn.bias_add(x, b)
    return tf.nn.relu(x)

#Define Max-Pooling Operation
def maxpool2d(x, k):
    return tf.nn.max_pool(x, ksize = [1, k, k, 1], strides = [1, k, k , 1], padding = "SAME")

#Weights
weights = {
    'wc1': tf.get_variable('W0', shape = (3, 3,  1,  batch_Size), initializer = tf.contrib.layers.xavier_initializer()),
    'wc2': tf.get_variable('W1', shape = (3, 3, batch_Size,  batch_Size2), initializer = tf.contrib.layers.xavier_initializer()),
    'wc3': tf.get_variable('W2', shape = (3, 3, batch_Size2, batch_Size4), initializer = tf.contrib.layers.xavier_initializer()),
    'wd1': tf.get_variable('W3', shape = (5*5*batch_Size4, batch_Size4), initializer = tf.contrib.layers.xavier_initializer()),
    'out': tf.get_variable('W6', shape = (batch_Size4, num_Classes), initializer = tf.contrib.layers.xavier_initializer())
}

#Biases
biases = {
    'bc1': tf.get_variable('B0', shape = (batch_Size), initializer = tf.contrib.layers.xavier_initializer()),
    'bc2': tf.get_variable('B1', shape = (batch_Size2), initializer = tf.contrib.layers.xavier_initializer()),
    'bc3': tf.get_variable('B2', shape = (batch_Size4), initializer = tf.contrib.layers.xavier_initializer()),
    'bd1': tf.get_variable('B3', shape = (batch_Size4), initializer = tf.contrib.layers.xavier_initializer()),
    'out': tf.get_variable('B4', shape = (num_Classes), initializer = tf.contrib.layers.xavier_initializer())
}

#Network Architecture
def conv_net(x, weights, biases):

    #Convolutional Layers
    conv1 = conv2d(x, weights['wc1'], biases['bc1'])
    conv1 = maxpool2d(conv1, k = 5)

    conv2 = conv2d(conv1, weights['wc2'], biases['bc2'])
    conv2 = maxpool2d(conv2, k = 4)

    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'])
    conv3 = maxpool2d(conv3, k = 2)

    #Fully Connected Layer
    fcl = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    fcl = tf.add(tf.matmul(fcl, weights['wd1']), biases['bd1'])
    fcl = tf.nn.relu(fcl)

    #Output, class prediction
    out = tf.add(tf.matmul(fcl, weights['out']), biases['out'])
    return out

#Loss and Optimizer Nodes
pred = conv_net(x, weights, biases)
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits = pred, labels = y))
optimizer = tf.train.AdamOptimizer(learning_rate = learning_Rate).minimize(cost)

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

    while i <= epochs:
        for batch in range(len(train_X)//batch_Size):
            print("Batch " + str(batch+1))
            batch_x = train_X[batch*batch_Size:min((batch+1)*batch_Size, len(train_X))]
            batch_y = train_Labels[batch*batch_Size:min((batch+1)*batch_Size, len(train_Labels))]

            #Calculate batch loss and accuracy
            opt = sess.run(optimizer, feed_dict = {x: batch_x, y: batch_y})
            loss, acc = sess.run([cost, accuracy], feed_dict = {x: batch_x, y: batch_y})

        #Calculate loss and training accuracy
        print("Epoch: " + str(i) + " out of " + str(epochs))
        print("Loss = " + \
                        "{:.6f}".format(loss) + ", Training Accuracy = " + \
                        "{:.5f}".format(acc))
        print("Optimization finished!")

        #Calculate and append testing accuracy
        print("Calculating testing accuracy...")
        test_acc, valid_loss = sess.run([accuracy, cost], feed_dict = {x: test_X, y: test_Labels})
        train_loss.append(loss)
        test_loss.append(valid_loss)
        train_accuracy.append(acc)
        test_accuracy.append(test_acc)

        print("Testing Accuracy: ", "{:.5f}".format(test_acc))
        learning_Rate = learning_Rate * 0.85

        #a = 1
        #while a == 1:
        #    another = input("Would you like to run another epoch? (yes/no) ")
        #    if another == "yes":
        #        a = 0
        #        epochs += 1
        #        i =+ 1
        #    elif another == "no":
        #        a = 0
        #    else:
        #        print("Please answer 'yes' or 'no'.")

        i =+ 1

    summary_writer.close()
