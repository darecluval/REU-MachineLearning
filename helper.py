# Machine Learning - Viola D - D Major
# Clare DuVal - MSU/Clemson
# Updated: August 23, 2019
import tensorflow as tf

class helper:

    #Define Convolutional layer
    def conv2d(x, w, b, strides = 1):
        x = tf.nn.conv2d(x, w, strides = [1, strides, strides, 1], padding = "SAME")
        x = tf.nn.bias_add(x, b)
        return tf.nn.relu(x)

    #Define Max-Pooling Operation
    def maxpool2d(x, k):
        return tf.nn.max_pool(x, ksize = [1, k, k, 1], strides = [1, k, k , 1], padding = "SAME")

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
